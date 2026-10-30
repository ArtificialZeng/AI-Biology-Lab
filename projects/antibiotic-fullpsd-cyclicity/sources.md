# Sources and exact relationship

Checked28October2026UTC. This is a bounded source screen, not exhaustive priority clearance.

1. Athreya GS, Gokhale CS, Verma P (2025), Coexistence rules for small, antagonistically interacting microbial communities. PLoS Computational Biology21(12):e1013763. https://doi.org/10.1371/journal.pcbi.1013763
   Publisher: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013763
   Native Eqs(1),(2),(4)-(6) define normalized discrete abundance dynamics, constant metabolic-cost rates, Poisson inhibition and independent-antibiotic inclusion-exclusion survival. Results use interior fixed points with Jacobian spectral radius below1. The Summary conjectures broader channel-cyclicity necessity after finite topology searches at two cost settings and selected sweeps. This note supports a restricted natural subclass; it neither corrects the model nor resolves the full conjecture.
2. Author code: https://github.com/tecoevo/antibiotics_biodiversity
   Inspected current head97e3d7d50fbe5428c5de72f51bd7512b7d5f1618, with native model.m/symmodel.m unchanged from release e3d9f192af4175908d4b01f1cf3b0dfda5342aea.
   Native map: https://github.com/tecoevo/antibiotics_biodiversity/blob/97e3d7d50fbe5428c5de72f51bd7512b7d5f1618/SimulationCode/mixed_inhibition_zone/functions/model.m
   Current SearchSpaceOfCommunities_mod.m filters PSD_fraction=1, motivating the every-channel-PSD class. Its cyclicity_degree is a motif-linked antibiotic-coverage metric, distinct from the paper's circulant-column pattern. No correction or theorem about that newer metric is claimed.
3. Kelsic ED etal(2015), Counteraction of antibiotic production and degradation stabilizes microbial communities. Nature521:516-519. https://doi.org/10.1038/nature14485
   Supplies the mixed-inhibition mechanism and cyclic coexistence baseline. The inspected primary Methods do not supply this arbitrary-channel fixed-role exclusion proof.
4. Wright ES, Vetsigian KH(2016), Inhibitory interactions promote frequent bistability among competing bacteria. Nature Communications7:11274. https://doi.org/10.1038/ncomms11274
   Its necessary bistable-pair condition is a distinct result. For columns PSD/PDS, g=1,c_P=1/5,c_D=1/10,c_R=1/20,K_P=log2 and any positive K_D, the PP/SD pair has invasion multipliers3/4 and2/3. Thus a no-cyclic-pair full-PSD topology can satisfy that prior necessary condition, which alone cannot exclude its coexistence. This comparison is not a new coexistence construction.

In this note, a cyclic pair consists of two distinct column words related by one or two rotations. Column ordering and strain relabeling do not change existence of such a pair. Duplicate channels are permitted but do not themselves count. With exactly one P,S,D in each three-entry column, there is no ambiguity from the source's near-circulant D/R examples.

Effective growth rates may be any positive constants for the impossibility proof. Original common-basal-growth, strictly ordered-cost parameterizations with positive rates are a subset. Artificially fitted positive rates in the60 numerical fixtures are not presented as physical parameter estimates; the two separate exact examples do satisfy the original cost order.

No inspected later primary resolution was found in bounded exact-title/DOI/cyclicity searches or current author issues/PRs. Absence of a hit is not proof that this theorem is new. No source catalogue-membership or all-literature claim is made.
