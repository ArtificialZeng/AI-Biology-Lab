# Source and licensing

- [COBRApy frozen source](https://github.com/opencobra/cobrapy/tree/a9c5ced29670bed5035fb19e87e7b0b588adb62d) is offered under LGPL-2.0-or-later OR GPL-2.0-or-later. [LICENSE.upstream](LICENSE.upstream) is the unchanged upstream license file containing both license texts.
- `upstream.patch` contains modifications and context from the upstream source. It is distributed under the upstream GPL-2.0-or-later option, preserving those terms. The lab's MIT license does not relicense this derivative patch.
- The independent review tests, reproduction wrapper, exact finite-family verifier, documentation and receipts were written with AI assistance. Their original material is covered by the lab's MIT license. No SPECTRA fixture code is included. Prior negative-only coverage from [PR1495](https://github.com/opencobra/cobrapy/pull/1495) is credited in `SOURCE.md`.
- Dependencies are installed separately and retain their respective licenses. No solver binaries, complete upstream source tree or third-party paper PDFs are bundled. The original FASTCC paper and related work are linked in `SOURCE.md`.
