# Licensing and attribution

GLPK 5.0 is GNU software under GPL-3.0-or-later; the unchanged release license is in LICENSE.glpk. The context-bearing native-guard.patch is a derivative distributed under those upstream terms. The lab MIT license does not relicense it.

The two C enumeration harnesses, reproduction wrapper, explanations and computed receipts are original lab contributions under the repository MIT license. Linking the harnesses to GPL GLPK creates temporary local executable tests; no solver binary or complete third-party source tree is distributed here. The official source must be supplied separately.

AI assisted source isolation, implementation and a separate skeptical review. The independently written second harness found both additional validation evidence and an unresolved tolerance family. These checks are not external human peer review or a solver-wide correctness proof.

Optlang is Apache-2.0 software; its unchanged license is in LICENSE.optlang. The context-bearing wrapper/optlang-mip-validation.patch and wrapper/test_glpk_mip_validation.py are distributed under Apache-2.0, matching that upstream project. The lab MIT license does not relicense upstream context. The original reproduce_wrapper.py and wrapper explanations remain MIT. No Optlang source tree, Python dependency wheel or modified GLPK library is bundled in this project.
