"""Reaction string parser: Complex and Reaction data types, canonical rate-key normalization."""
import re
from dataclasses import dataclass

# frozenset of (species_name, stoichiometric_coefficient) pairs
Complex = frozenset[tuple[str, int]]

# A valid term is an optional integer coefficient followed immediately by a
# species name that starts with a letter: "2B", "B", "miR21_H1".
# Coeff and name may optionally be separated by whitespace: "2 B" is allowed.
_TERM_RE = re.compile(r'^(\d+)?\s*([A-Za-z][A-Za-z0-9_]*)$')


def parse_complex(text: str) -> Complex:
    """Parse "A + 2B" → frozenset({('A',1), ('B',2)})."""
    text = text.strip()
    if not text:
        return frozenset()
    terms: dict[str, int] = {}
    for raw_term in text.split("+"):
        token = raw_term.strip()
        if not token:
            continue
        m = _TERM_RE.match(token)
        if not m:
            raise ValueError(
                f"Invalid species name '{token}': must be [coeff]Name where "
                f"Name starts with a letter and contains only letters, digits, underscores"
            )
        coeff = int(m.group(1)) if m.group(1) else 1
        name = m.group(2)
        terms[name] = terms.get(name, 0) + coeff
    if not terms:
        raise ValueError(f"Could not parse complex: '{text}'")
    return frozenset(terms.items())


def canonical_rate_key(reactants: Complex, products: Complex) -> str:
    """Produce a deterministic string 'A + 2B -> C' for rate dict lookup."""
    def fmt(c: Complex) -> str:
        parts = []
        for name, coeff in sorted(c, key=lambda x: x[0]):
            parts.append(f"{coeff}{name}" if coeff > 1 else name)
        return " + ".join(parts)
    return f"{fmt(reactants)} -> {fmt(products)}"


@dataclass(frozen=True)
class Reaction:
    reactants: Complex
    products: Complex
    rate_key: str


def parse_reaction_string(text: str) -> list[Reaction]:
    """
    Parse one reaction string. Returns 1 Reaction for '->' or 2 for '<->'.
    Spaces around arrows are optional.
    """
    text = text.strip()
    if "<->" in text:
        left, right = text.split("<->", 1)
        r = parse_complex(left)
        p = parse_complex(right)
        return [
            Reaction(r, p, canonical_rate_key(r, p)),
            Reaction(p, r, canonical_rate_key(p, r)),
        ]
    elif "->" in text:
        left, right = text.split("->", 1)
        r = parse_complex(left)
        p = parse_complex(right)
        return [Reaction(r, p, canonical_rate_key(r, p))]
    else:
        raise ValueError(f"Reaction string must contain '->' or '<->': '{text}'")


def parse_reactions(strings: list[str]) -> list[Reaction]:
    """Parse a list of reaction strings into directed Reaction objects."""
    reactions = []
    for s in strings:
        reactions.extend(parse_reaction_string(s))
    return reactions


def normalize_rate_key(key: str) -> str:
    """
    Normalize a user-provided rate key string to canonical form.
    Handles arbitrary spacing and species order on each side.
    """
    if "->" not in key:
        raise ValueError(f"Rate key must contain '->': '{key}'")
    left, right = key.split("->", 1)
    r = parse_complex(left)
    p = parse_complex(right)
    return canonical_rate_key(r, p)


def reduce_complex(c: Complex, chemostatted: set[str]) -> Complex:
    """Remove chemostatted species from a complex for the reduced CRNT graph."""
    return frozenset((name, coeff) for name, coeff in c if name not in chemostatted)
