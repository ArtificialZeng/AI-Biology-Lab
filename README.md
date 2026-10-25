# AI Biology Lab

**AI-assisted discovery in mathematical and computational biology, with inspectable proofs, code, and evidence.**

[中文](README.zh-CN.md) · [First result](projects/mixed-phosphorylation/README.md) · [Research roadmap](docs/RESEARCH_ROADMAP.md) · [For reviewers](docs/VALIDATION.md) · [媒体资料](docs/MEDIA_BRIEF.zh-CN.md)

Led by [Zijian Zeng / ArtificialZeng](https://github.com/ArtificialZeng). The lab pursues questions a small computational team can investigate and others can independently check: biochemical oscillations, molecular switches, reaction networks, and rigorous scientific computing.

The ambition is to produce discoveries worth reporting. Each public claim is accompanied by its scope, evidence, and unresolved limitations.

## First project: can a biochemical oscillator resist changes in enzyme concentration?

The mixed-mechanism phosphorylation model studied by Conradi, Mincheva, and Shiu has a unique positive equilibrium. Their 2019 paper asks whether it can remain unstable as total kinase grows arbitrarily large.

**Our exact result:** at the paper's Table 1 reaction rates, with total phosphatase **8** and total substrate **40**, the equilibrium has exactly **two unstable and four stable eigenvalues** on its six-dimensional compatibility class for every finite

$$K_{\mathrm{tot}}\ge \frac{4660166}{209569}\approx22.2369.$$

In particular, the result holds for every **K ≥ 23**. The bound is sufficient and is not claimed to be sharp. Rational equilibrium parametrization and **125 strictly positive rational Bernstein coefficients** give a certificate for the entire unbounded range, rather than a finite parameter scan.

| Evidence level | Current status |
| --- | --- |
| Equilibrium instability over the stated range | Exact proof and executable certificate |
| Periodic trajectories and concentration corrections | Floating-point numerical evidence |
| Stable periodic orbits for every sufficiently large K | Open in this project |
| Original Conjecture 6.2 | First assertion addressed; full conjecture unresolved |
| Novelty and priority | Unconfirmed; directly related 2023/2025 conference work exists |
| External human peer review / biological experiments | None reported |

**Unstable equilibrium does not by itself prove a periodic orbit.** The current result concerns a published mathematical model; it does not establish new behavior in living cells.

[Read the five-page note](projects/mixed-phosphorylation/manuscript/main.pdf) · [Read the proof](projects/mixed-phosphorylation/evidence/finite_kinase_proof.md) · [Inspect the certificate](projects/mixed-phosphorylation/evidence/finite_kinase_sign_certificate.json)

![Numerical concentration effects: waveforms, period and amplitude](projects/mixed-phosphorylation/evidence/concentration_effects.png)

*Numerical illustration, not a rigorous periodic-orbit certificate. The curves show finite-concentration effects and convergence toward a limiting model.*

## Verify the result

Python 3.12 or 3.13. The exact verification requires only SymPy and the Python standard library; no GPU, API key, or model weights are needed.

```sh
git clone https://github.com/ArtificialZeng/AI-Biology-Lab.git
cd AI-Biology-Lab
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-exact.txt
python scripts/reproduce.py verify
```

This runs three checks: reconstruction of the manuscript matrix, a direct nine-equation formulation, and a second reduction with an exact tangent-space identity and Routh root count. Each uses rational arithmetic. Runs take place in a temporary copy; results are saved under ignored `artifacts/` and the archived evidence stays unchanged. See [reproduction instructions](docs/REPRODUCIBILITY.md).

For numerical reproduction and figures:

```sh
python -m pip install -r requirements.txt
python scripts/reproduce.py numerics
python scripts/reproduce.py figure
```

## What we are working toward

1. Certify an attracting periodic orbit, including rigorous stability bounds.
2. Prove persistence for an unbounded range of finite kinase totals.
3. Bound how period and amplitude change with concentration.
4. Identify an interpretable parameter region and where robustness fails.

These are research targets, not completed discoveries. [The roadmap](docs/RESEARCH_ROADMAP.md) defines their acceptance criteria. [External validation](docs/VALIDATION.md) comes before claims of priority or resolution of the full conjecture.

## Repository map

```text
projects/mixed-phosphorylation/  First research project, proof, data, and note
scripts/reproduce.py            Reproduction entry point
docs/                          Research plan, validation, and media information
.github/workflows/             Continuous exact verification
```

The project uses AI for derivation, coding, search, and review. Separate AI review passes and independent implementations are useful checks; they are not external human peer review. See [AI use and provenance](docs/AI_USE.md).

## Sources and credit

- Conradi, Mincheva, and Shiu (2019), [*Emergence of oscillations in a mixed-mechanism phosphorylation system*](https://doi.org/10.1007/s11538-019-00580-6), [arXiv version used](https://arxiv.org/html/1809.02886v2).
- Suwanmajo and Krishnan (2015), [*Mixed mechanisms of multi-site phosphorylation*](https://doi.org/10.1098/rsif.2014.1405).
- Related work requiring priority clarification: [Grasse, SBHD 2023, poster 48](https://www.sbhdberlin.org/pdf/SBHD2023_Poster_Session_Program.pdf); [Grasse and Conradi, BIOMATH 2025 abstract](https://biomath.math.bas.bg/biomath/index.php/bmcs/article/download/1787/1594/7828).

Third-party papers are linked, not redistributed. Original repository material is available under the [MIT license](LICENSE). Citation metadata is in [CITATION.cff](CITATION.cff).
