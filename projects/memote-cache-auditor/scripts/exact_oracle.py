"""Elementary exact mass constraints, independent of MEMOTE and LP solvers."""

from fractions import Fraction


def exact_oracle():
    rows = []
    for multiplier in (2, 3):
        stoich = ((-1, 1), (-1, multiplier))
        determinant = stoich[0][0] * stoich[1][1] - stoich[0][1] * stoich[1][0]
        witness = (Fraction(1), Fraction(1))
        assert sum(c * m for c, m in zip(stoich[0], witness)) == 0
        assert all(m > 0 for m in witness)
        assert determinant != 0  # The pair has only the zero conserved mass.
        rows.append(
            {
                "multiplier": multiplier,
                "transposed_stoichiometry": [list(row) for row in stoich],
                "determinant": determinant,
                "single_reaction_positive_mass": [str(m) for m in witness],
                "both_reactions_positive_mass_exists": False,
            }
        )
    return {"arithmetic": "integer/Fraction", "cases": rows}
