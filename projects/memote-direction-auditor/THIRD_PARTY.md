# Source and license notices

The pinned MEMOTE helper and consistency module are from `opencobra/memote` commit `bf76477b4087feefc8ad429ebd492e30ef10eebc` (release 0.17.0). Their original copyright notices and Apache-2.0 license are retained in `fixtures/`. The proposed helper changes are offered under the same upstream Apache-2.0 terms. The lab's original test code is MIT-licensed.

COBRApy, optlang and GLPK are runtime dependencies; they are installed through their distributions and are not vendored here. Their own license terms apply. No upstream package is modified in place: the verifier imports original and patched copies from temporary directories.
