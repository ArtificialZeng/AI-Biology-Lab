# ETE RF polytomy correction on the common-leaf induced tree

Current ETE can return -1 for two presentations of the same unrooted quartet:
`((A,B),(C,D));` and `(A,B,(C,D));`, with `unrooted_trees=True` and
`correct_by_polytomy_size=True`. Their ordinary RF distance is zero.
[Issue745](https://github.com/etetoolkit/ete/issues/745) already reports the
related treatment of a degree-three master root as a polytomy and is credited.

The three-file patch replaces the raw child-count correction by the number of
missing resolutions on the common-label induced topology. It handles rooted
and unrooted comparisons, different leaf sets, degree-two/three master roots,
unary paths and custom leaf properties. Two genuinely unresolved induced trees
still raise the existing exception. No RF result is clipped. The new regression
module is added to the official default fast runner.

This is a scientific-software correction, not a new distance or biological
discovery. Nonzero support filtering already fails in the current baseline and
remains outside this repair; 24 unchanged interactions are recorded as known
failures, not successful functionality. No new input prohibition is introduced.
The corrected score is not asserted to minimize RF over all binary refinements.

## Offline verification

Python3.12 or newer, standard library only:

```sh
python -B verify.py verify
python -B -O verify.py verify
python -B verify.py regenerate --output-dir /absolute/new/external-packet
```

This verifies the complete file/hash inventory and reruns 534 exact graph
identities. It does not execute ETE. The offline certificate labels that limit.
Verification writes only to external temporary directories. Extra or missing
files, symlinks and altered bytes fail. Regeneration requires an exclusive new
external directory and creates a complete replayable packet.

## Actual native replay

Use a fresh virtual environment and a C compiler, then:

```sh
python -m pip install -r requirements-native.txt
python -B verify.py native --output-dir /absolute/new/external-native-run
```

The command fetches the exact official source archive at commit
`9562dfb6a02795dfda3975be5025f83b8dc884b1`, verifies its SHA256, applies the
byte-bound patch, compiles both genuine Cython versions, and runs the original
frozen audit workers. For an already downloaded exact archive add
`--source-archive /absolute/archive.tar.gz`; the same checksum remains mandatory.
Network is needed for dependencies and source unless already available. Build
logs, source trees, compiled runtimes and outputs stay in the external directory.
The project itself may be read-only. No model, API key, GPU or paid service is used.

Native replay requires the baseline's known 10 failed /8 errored regression
subcases, followed by six successful focused methods on the candidate. It runs
189,204 actual API comparisons against an independent graph pruning/degree
oracle, including 14,991 expected two-polytomy exceptions, and verifies 24
unchanged inherited support observations. Worker assertions run without Python
optimization even when the outer wrapper uses `-O`.

The repository-required original default fast suite was also run locally:
both native baseline and candidate reported 71 passed/10 failed with identical
failed test IDs. Nine failures lacked optional lapjv; one differed from the
downloaded ortholog-annotation fixture. The candidate populated a shared ETE
fixture cache which the baseline reused. These are not clean full-suite results,
and the updated default runner's newly added module was tested separately.
Interactive GUI and slow suites were not run. Local Python was3.13.5; hosted
Python3.12 validation must be read from the actual PR checks.

GPL-3.0-or-later licensing and original ETE attribution are retained. The work
used substantial AI assistance and fresh-context AI review, which is not human
peer review, maintainer acceptance or formal proof-assistant certification.
