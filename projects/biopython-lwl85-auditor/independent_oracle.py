"""Independent LWL85 oracle for the bounded single-difference fixtures.

No Biopython imports/helpers. Standard genetic-code map uses NCBI table 1 in
TCAG order. Threefold or cross-degeneracy changes are outside this fixture oracle.
"""
from fractions import Fraction
from math import log1p
from itertools import product

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = dict(zip(("".join(c) for c in product(BASES, repeat=3)), AMINO_ACIDS))


def folds(codon):
    aa = GENETIC_CODE[codon]
    if aa == "*":
        raise ValueError("Stop codon excluded")
    result = []
    for position in range(3):
        synonymous = sum(GENETIC_CODE[codon[:position] + b + codon[position+1:]] == aa
                         for b in BASES)
        if synonymous not in (1, 2, 4):
            raise ValueError("Threefold-site conventions outside this oracle")
        result.append({1: 0, 2: 2, 4: 4}[synonymous])
    return result


def estimate(sequence1, sequence2):
    if len(sequence1) != len(sequence2) or len(sequence1) % 3:
        raise ValueError("Require equally long codon-aligned strings")
    sites = {0: Fraction(0), 2: Fraction(0), 4: Fraction(0)}
    transition = {0: 0, 2: 0, 4: 0}
    transversion = {0: 0, 2: 0, 4: 0}
    for offset in range(0, len(sequence1), 3):
        c1, c2 = sequence1[offset:offset+3], sequence2[offset:offset+3]
        f1, f2 = folds(c1), folds(c2)
        for cls in f1 + f2:
            sites[cls] += Fraction(1, 2)
        differences = [j for j in range(3) if c1[j] != c2[j]]
        if len(differences) > 1:
            raise ValueError("Only zero/single-nucleotide codon differences supported")
        for j in differences:
            if f1[j] != f2[j]:
                raise ValueError("Cross-degeneracy convention outside this oracle")
            cls = f1[j]
            if {c1[j], c2[j]} in ({"A", "G"}, {"C", "T"}):
                transition[cls] += 1
            else:
                transversion[cls] += 1
    A, B, fractions = {}, {}, {}
    for cls in (0, 2, 4):
        if sites[cls] <= 0:
            raise ValueError("All site classes must be populated")
        P, Q = Fraction(transition[cls], sites[cls]), Fraction(transversion[cls], sites[cls])
        if 2*P+Q >= 1 or 2*Q >= 1:
            raise ValueError("Saturated fixture")
        A[cls] = -.5*log1p(-float(2*P+Q)) + .25*log1p(-float(2*Q))
        B[cls] = -.5*log1p(-float(2*Q))
        fractions[str(cls)] = {"P": str(P), "Q": str(Q)}
    L0, L2, L4 = map(float, (sites[0], sites[2], sites[4]))
    dS = 3*(L2*A[2]+L4*(A[4]+B[4]))/(L2+3*L4)
    dN = 3*(L2*B[2]+L0*(A[0]+B[0]))/(2*L2+3*L0)
    return {"dN": dN, "dS": dS, "omega": dN/dS if dS else None,
            "sites": {str(k): str(v) for k, v in sites.items()},
            "transitions": transition, "transversions": transversion,
            "fractions": fractions, "A": A, "B": B}


def fixtures():
    cases = []
    for name, phe, ala, changes in [
        ("asymmetric_synonymous", 20, 5, {0:"TTC",1:"TTC"}),
        ("asymmetric_omega",20,5,{0:"TTC",1:"TTC",2:"TAT"}),
        ("asymmetric_twofold_transversion",20,5,{0:"TTC",1:"TTC",2:"TTG",3:"TAT"}),
        ("balanced_negative_control",10,10,{0:"TTC",1:"TTC",2:"TAT"}),
    ]:
        before = ["TTT"]*phe+["GCT"]*ala
        after = before.copy()
        for j,codon in changes.items():
            after[j] = codon
        cases.append({"name":name,"sequence1":"".join(before),"sequence2":"".join(after)})
    return cases

if __name__ == "__main__":
    import json
    print(json.dumps([{**c,"oracle":estimate(c["sequence1"],c["sequence2"])} for c in fixtures()],sort_keys=True,indent=2))
