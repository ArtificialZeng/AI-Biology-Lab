# Example 8 Jacobian correction and stability verification

In Example 8 of Eilertsen, Romanovski, Schnell and Walcher's [*Lumping of reaction networks: Generic and critical parameters*, arXiv:2606.28895v2](https://arxiv.org/abs/2606.28895v2), the displayed Jacobian entry `(2,1)` on printed p.38 is `-k2*c2`. Differentiating the preceding equation gives `-k3*c2`.

The printed expression has an exact positive eigenvalue at an admissible positive rate choice. The corrected derivative is Hurwitz for all eight fast reaction rates strictly positive and `c1,c2 >= 0`, including zero concentration boundaries. [The proof](PROOF.md) supplies the full argument; the scripts check its algebra exactly.

This is a source-specific correction and verification of the intended local stability statement. The source's broad future question about useful all-time approximation criteria remains unresolved. The companion repository represents the reactions correctly; no production-code defect or changed biological behavior is established here.

From the repository root, with Python 3.12 or 3.13:

```sh
python -m pip install -r projects/lumping-example8-auditor/requirements.txt
python -B projects/lumping-example8-auditor/replay.py
```

The replay verifies archived hashes, runs two separate exact implementations in ordinary and optimized Python modes, and compares their output with four archived JSON reports. It writes no reports or archive files. To run either checker directly:

```sh
python -B projects/lumping-example8-auditor/verify.py
python -B projects/lumping-example8-auditor/independent_verify.py
```

Both default to one JSON document on stdout. An explicit `--output NEW_FILE` requests exclusive creation of a new report; existing files are refused.

Read [the source record](SOURCE.md), [scope](SCOPE.md), [verification limits](VERIFICATION.md), and [credit and licensing](THIRD_PARTY.md). The self-replication model is due to Gijima and Peacock-López and the cited source authors. Prepared by Zijian Zeng with AI assistance; separate AI review and symbolic verification are not external human peer review or proof-assistant formalization.
