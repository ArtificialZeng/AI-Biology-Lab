# Validated limiting oscillator — candidate study

This study attempts to certify a locally attracting periodic orbit of the published model's five-dimensional high-kinase limit at P=8, S=40. It uses fixed-source CAPD interval integration, a four-dimensional Poincare section, Krawczyk inclusion, and a transverse contraction bound. These are standard verification methods.

**Status:** development inequalities passed on Apple Silicon; Linux x86/FILIB reproduction and fresh mathematical review are pending. No novelty or full-Conjecture-6.2 claim is made. The original paper already numerically observed oscillations, and directly overlapping 2023/2025 conference work needs priority clarification. See [candidate argument](candidate_lemma.md).

Run `bash build_and_run.sh` in this folder. Dependencies build locally and are not installed globally. The script downloads CAPD revision `37057e5899bccaa6cd5701057b7e6ea5e232ef8b`; CMake, a C++17 compiler, and Git are required. Mac native-interval output is development evidence only. The Linux workflow uses x86 and FILIB. No binary dependency is redistributed.

The original checker source is MIT-licensed. CAPD is GPL-3.0; consult its pinned source's COPYING before redistributing linked binaries. The numerical seed and preconditioner are only proposed inputs; every decisive inequality must be recomputed by the interval checker.
