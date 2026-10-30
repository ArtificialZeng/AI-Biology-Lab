"""SBML (Systems Biology Markup Language) import/export for CRNetwork.

This is an *optional* interoperability layer.  It requires ``python-libsbml``;
if that package is not installed, importing the helpers raises a clear
:class:`ImportError` with an install hint, but the rest of ``mantis`` keeps
working without it.

The mapping between mantis and SBML is mass-action-centric:

* Each directed :class:`~mantis.parsing.Reaction` becomes one irreversible SBML
  reaction.  A reversible ``A <-> B`` round-trips as two irreversible reactions
  (matching mantis' internal representation).
* Each reaction carries a mass-action kinetic law ``k_i * Π reactant^coeff``
  with the rate constant stored as a global parameter ``k_i``.
* Chemostatted species are exported as SBML *boundary* species
  (``boundaryCondition=True``) holding their fixed concentration; on import,
  boundary / constant species are recovered as chemostatted.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .network import CRNetwork


_INSTALL_HINT = (
    "SBML support requires the optional 'python-libsbml' package. "
    "Install it with:  pip install mantis-delta[sbml]   (or  pip install python-libsbml)"
)


def _require_libsbml():
    """Import and return the libsbml module, or raise a helpful ImportError."""
    try:
        import libsbml  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without the dep
        raise ImportError(_INSTALL_HINT) from exc
    return libsbml


def _check(libsbml, code, what: str) -> None:
    """Raise if a libsbml call returned a failure status code."""
    if code is None:
        return
    if isinstance(code, int) and code != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise RuntimeError(f"libsbml error while {what}: status {code}")


# ── Export ────────────────────────────────────────────────────────────────────

def network_to_sbml(net: "CRNetwork", *, level: int = 3, version: int = 2) -> str:
    """Serialize a :class:`CRNetwork` to an SBML document string (Level 3 v2)."""
    libsbml = _require_libsbml()

    doc = libsbml.SBMLDocument(level, version)
    model = doc.createModel()
    model.setId("mantis_model")

    # A single well-mixed compartment of unit size.
    comp = model.createCompartment()
    comp.setId("c")
    comp.setConstant(True)
    comp.setSize(1.0)
    comp.setSpatialDimensions(3)

    chem = net.chemostatted  # name -> fixed concentration

    # All species names appearing anywhere in the reactions (dynamic + chemostatted).
    names: list[str] = []
    seen: set[str] = set()
    for rxn in net._reactions:
        for name, _ in list(rxn.reactants) + list(rxn.products):
            if name not in seen:
                seen.add(name)
                names.append(name)
    names.sort()

    for name in names:
        sp = model.createSpecies()
        sp.setId(name)
        sp.setCompartment("c")
        sp.setHasOnlySubstanceUnits(False)
        if name in chem:
            sp.setInitialConcentration(float(chem[name]))
            sp.setBoundaryCondition(True)
            sp.setConstant(True)
        else:
            sp.setInitialConcentration(0.0)
            sp.setBoundaryCondition(False)
            sp.setConstant(False)

    # One global rate parameter per reaction, plus one irreversible reaction.
    for i, rxn in enumerate(net._reactions, start=1):
        kid = f"k_{i}"
        param = model.createParameter()
        param.setId(kid)
        param.setConstant(True)
        # Only assign a value for rates that were actually supplied; leaving it
        # unset lets ``from_sbml`` distinguish "unknown" from a literal 1.0.
        if rxn.rate_key in net._rates:
            param.setValue(float(net._rates[rxn.rate_key]))

        r = model.createReaction()
        r.setId(f"r_{i}")
        r.setReversible(False)
        # The reaction's human-readable rate key is preserved as the SBML name.
        r.setName(rxn.rate_key)

        for name, coeff in sorted(rxn.reactants, key=lambda x: x[0]):
            sr = r.createReactant()
            sr.setSpecies(name)
            sr.setStoichiometry(float(coeff))
            sr.setConstant(True)
        for name, coeff in sorted(rxn.products, key=lambda x: x[0]):
            sp_ref = r.createProduct()
            sp_ref.setSpecies(name)
            sp_ref.setStoichiometry(float(coeff))
            sp_ref.setConstant(True)

        # Mass-action kinetic law: k_i * Π species^coeff (chemostatted included).
        factors = [kid]
        for name, coeff in sorted(rxn.reactants, key=lambda x: x[0]):
            factors.extend([name] * int(coeff))
        formula = " * ".join(factors)
        kl = r.createKineticLaw()
        math_ast = libsbml.parseL3Formula(formula)
        if math_ast is None:  # pragma: no cover - formula is machine-generated
            raise RuntimeError(f"Failed to parse kinetic-law formula: {formula!r}")
        kl.setMath(math_ast)

    writer = libsbml.SBMLWriter()
    return writer.writeSBMLToString(doc)


def write_sbml(net: "CRNetwork", path: str) -> None:
    """Write a :class:`CRNetwork` to an SBML file at ``path``."""
    xml = network_to_sbml(net)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(xml)


# ── Import ────────────────────────────────────────────────────────────────────

def _stoich_terms(species_refs) -> list[tuple[str, int]]:
    """Extract (species_id, integer_coeff) pairs from an SBML reactant/product list."""
    terms: list[tuple[str, int]] = []
    for ref in species_refs:
        coeff = ref.getStoichiometry()
        if coeff is None or (isinstance(coeff, float) and coeff != coeff):  # NaN
            coeff = 1.0
        terms.append((ref.getSpecies(), int(round(coeff))))
    return terms


def _complex_str(terms: list[tuple[str, int]]) -> str:
    """Render [(name, coeff), ...] as 'a + 2 b' for a reaction string."""
    parts = []
    for name, coeff in terms:
        parts.append(f"{coeff} {name}" if coeff != 1 else name)
    return " + ".join(parts) if parts else ""


def _collect_param_names(ast, libsbml) -> list[str]:
    """Depth-first collect AST_NAME identifiers from a libsbml math tree."""
    if ast is None:
        return []
    out: list[str] = []
    if ast.getType() == libsbml.AST_NAME:
        out.append(ast.getName())
    for i in range(ast.getNumChildren()):
        out.extend(_collect_param_names(ast.getChild(i), libsbml))
    return out


def _param_value(model, kinetic_law, name: str):
    """Look up a parameter value by id, preferring the reaction-local scope."""
    if kinetic_law is not None:
        lp = kinetic_law.getLocalParameter(name) if hasattr(kinetic_law, "getLocalParameter") else None
        if lp is not None and lp.isSetValue():
            return float(lp.getValue())
        p = kinetic_law.getParameter(name) if hasattr(kinetic_law, "getParameter") else None
        if p is not None and p.isSetValue():
            return float(p.getValue())
    gp = model.getParameter(name)
    if gp is not None and gp.isSetValue():
        return float(gp.getValue())
    return None


def _rates_for_reaction(model, reaction, species_ids: set[str], libsbml):
    """Best-effort extraction of (forward_value, reverse_value) rate constants.

    Returns a tuple; either element may be None when it cannot be determined.
    Reversible reactions whose math is ``forward - reverse`` are split on the
    top-level minus node so each side's lone parameter is attributed correctly.
    """
    kl = reaction.getKineticLaw()
    if kl is None:
        return None, None
    math_ast = kl.getMath()
    if math_ast is None:
        return None, None

    def lone_param(node):
        names = [n for n in _collect_param_names(node, libsbml) if n not in species_ids]
        # Keep only names that resolve to an actual parameter value.
        vals = [(n, _param_value(model, kl, n)) for n in names]
        vals = [(n, v) for n, v in vals if v is not None]
        if len(vals) == 1:
            return vals[0][1]
        return None

    if reaction.getReversible() and math_ast.getType() == libsbml.AST_MINUS \
            and math_ast.getNumChildren() == 2:
        fwd = lone_param(math_ast.getChild(0))
        rev = lone_param(math_ast.getChild(1))
        return fwd, rev

    return lone_param(math_ast), None


def network_from_sbml(source: str) -> "CRNetwork":
    """Construct a :class:`CRNetwork` from an SBML file path or XML string.

    ``source`` is treated as a file path if it points to an existing file,
    otherwise as a raw SBML/XML document string.
    """
    import os

    from .network import CRNetwork
    from .parsing import parse_reactions, normalize_rate_key

    libsbml = _require_libsbml()

    if os.path.exists(source):
        doc = libsbml.readSBMLFromFile(source)
    else:
        doc = libsbml.readSBMLFromString(source)
    if doc is None or doc.getModel() is None:
        raise ValueError("Could not parse SBML: no model found in document.")
    if doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) > 0:
        msg = doc.getErrorLog().toString()
        raise ValueError(f"SBML document has fatal errors:\n{msg}")

    model = doc.getModel()
    species_ids = {model.getSpecies(i).getId() for i in range(model.getNumSpecies())}

    # Chemostatted species: SBML boundary or constant species, with their fixed value.
    chemostatted: dict[str, float] = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        if sp.getBoundaryCondition() or sp.getConstant():
            if sp.isSetInitialConcentration():
                val = float(sp.getInitialConcentration())
            elif sp.isSetInitialAmount():
                val = float(sp.getInitialAmount())
            else:
                val = 0.0
            chemostatted[sp.getId()] = val

    reaction_strings: list[str] = []
    rates: dict[str, float] = {}
    for i in range(model.getNumReactions()):
        rxn = model.getReaction(i)
        reactants = _stoich_terms(
            [rxn.getReactant(j) for j in range(rxn.getNumReactants())]
        )
        products = _stoich_terms(
            [rxn.getProduct(j) for j in range(rxn.getNumProducts())]
        )
        lhs = _complex_str(reactants)
        rhs = _complex_str(products)
        arrow = "<->" if rxn.getReversible() else "->"
        rxn_str = f"{lhs} {arrow} {rhs}"
        reaction_strings.append(rxn_str)

        fwd, rev = _rates_for_reaction(model, rxn, species_ids, libsbml)
        if fwd is not None:
            rates[normalize_rate_key(f"{lhs} -> {rhs}")] = fwd
        if rev is not None:
            rates[normalize_rate_key(f"{rhs} -> {lhs}")] = rev

    reactions = parse_reactions(reaction_strings)
    return CRNetwork(reactions, rates=rates, chemostatted=chemostatted or None)
