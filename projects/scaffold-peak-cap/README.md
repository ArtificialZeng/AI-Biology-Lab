# Sharp peak cap for exclusive tri-ligand scaffold binding

This original code and proof packet concerns the newly posed question for three independent ligand sites and **exclusive** occupancy of sites 1 and 2. For all positive binding parameters and total ligand concentrations, the proposed theorem gives at most two strict local maxima as scaffold concentration varies. It does not address arbitrary numbers of sites and is not attributed to the source author as a conjecture.

The proof text reduces the stationary equation to a third-derivative sign. The checker reconstructs the resulting exact polynomial and verifies all 322 integer coefficients against the supplied certificate. Separate Fraction-only checks certify the known two-peak parameter example, establishing sharpness when combined with the cap proof. This finite algebra check supports the proof; the analytic reductions and Rolle argument must also be read. Non-exhaustive negative floating-point searches are not proofs and are not required to reproduce this packet.

A separate AI context independently reconstructed the derivative identity and all 322 coefficients, checked the root branches, complete positive domain and Rolle argument, and verified sharpness using fresh rational bounds. A separate source/value review accepted the bounded contribution while retaining literature-priority uncertainty. These are computational and mathematical checks, not human peer review; source-author verification is pending.

## Reproduce

Requires Python 3.12 or 3.13 and SymPy 1.13.3. Install the dependency with `python -m pip install -r requirements.txt`, then run:

```sh
python reproduce.py --work-dir /absolute/path/outside-this-project
python -O reproduce.py --work-dir /another/absolute/external/path
```

The wrapper verifies the file manifest and launches fresh scientific child interpreters with assertions enabled, even if the parent uses `-O` or `PYTHONOPTIMIZE` is set. It reconstructs the coefficient identity and checks both rational sharpness certificates. All generated logs and the result JSON go to the supplied external work directory; the project can be read-only. The manifest is an integrity check, not mathematical evidence by itself.

## Model and credit

The model and the previously known two-peak example come from Eduardo D. Sontag, [Dynamics and dose response in scaffold ligand binding](https://arxiv.org/html/2508.06599v6), version 6 (29 June 2026), Theorem 2, section 6.1, and Remark 3/Figure 8. The original contribution proposed here is the full positive-parameter cap for exclusive doubly occupied three-site scaffolds, not construction of the known sharpness example.

An earlier primary model is J. Yang and W. S. Hlavacek, [Scaffold-mediated nucleation of protein signaling complexes: elementary principles](https://arxiv.org/abs/1105.4688). Complete occupancy and marginal pair occupancy differ from the exclusive partial occupancy studied here. The bounded source audit found no equivalent exclusive cap, but does not establish exhaustive literature priority. No source-author endorsement is claimed. Third-party papers are linked rather than redistributed.

The MIT license applies only to the original text, code, and certificates in this packet. It does not relicense the linked sources or SymPy.
