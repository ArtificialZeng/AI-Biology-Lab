# Exact formula scope

Index the site counts as L=[L0,L2,L4] and corrections as A=[A0,A2,A4], B=[B0,B2,B4]. Original LWL85 aggregates

    dS = 3*(L2*A2 + L4*(A4+B4))/(L2+3*L4),
    dN = 3*(L2*B2 + L0*(A0+B0))/(2*L2+3*L0),
    Ai = -ln(1-2*Pi-Qi)/2 + ln(1-2*Qi)/4,
    Bi = -ln(1-2*Qi)/2.

The audited Biopython implementations used L4 in the first numerator term of both formulas. The patch replaces these two factors with L2 in each API. This does not audit saturation handling, ambiguous codons, pathway averaging, genetic-code conventions, other estimators or all possible alignments.

For the simple Phe/Ala fixture with no fourfold changes, A4=B4=B2=0, so the erroneous dS is L4/L2 times the formula dS. In the stated L2/L4=4 fixture dN is unchanged and dN/dS is four times the corrected ratio. This is an exact coefficient observation; numerical logarithms only evaluate the displayed values.
