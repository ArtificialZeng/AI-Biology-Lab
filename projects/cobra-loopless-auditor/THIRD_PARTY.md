# Source and licensing

- [COBRApy frozen source](https://github.com/opencobra/cobrapy/tree/a9c5ced29670bed5035fb19e87e7b0b588adb62d), public API documentation and tests supply the baseline. COBRApy offers LGPL-2.0-or-later OR GPL-2.0-or-later; the supplied license file includes both texts.
- `upstream.patch` modifies GPL-covered upstream files and is distributed under the upstream GPL-2.0-or-later option; [LICENSE.upstream](LICENSE.upstream) preserves the license. The lab's MIT license does not relicense those derivative patch contents.
- The independent contract test and reproduction wrapper were newly written with AI assistance. Their original code is MIT licensed. Dependencies are installed separately under their own licenses; no solver binaries or complete upstream source tree are bundled.
- [PR #1472](https://github.com/opencobra/cobrapy/pull/1472) was inspected for overlap. It proposes related formulation work but retained the same empty-index paths when checked. No claim is made about later commits or universal absence of another fix.
