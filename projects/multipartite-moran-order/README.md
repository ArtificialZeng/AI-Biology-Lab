# Mixed-Moran ordering on complete multipartite graphs

An optional exact count-chain oracle and a proposed subsidiary proof note for the fixed-mixture Moran model. Prepared by Shawn Zeng (Zijian Zeng), AI Biology Lab, UCSI University.

The native process chooses birth-death (Bd) with fixed probability lambda and death-Birth (dB) with probability1-lambda. On a complete multipartite graph, count mutants in each part. We prove two connected statements:

* For a fixed mixture, the native count generator preserves coordinatewise stochastic order at EVERY positive mutant fitness exactly when the process is pure dB, the graph has two parts, or all part sizes are equal.
* On complete bipartite and equal-part complete multipartite graphs, fixation probability is nondecreasing in positive fitness and initial part counts, for every fixed mixture. The known pure-dB comparison is included and credited.

The criterion is universal in fitness. It does not say an excluded graph fails at every individual fitness. Attractivity here means stochastic-order preservation for the native generator L=P-I, without a state-dependent clock change.

The proof classifies a count-coupling method, not the original general-graph fixation question. On unequal multipartite graphs, failure of this particular generator order does not refute fixation monotonicity or exclude other orders/time changes. The accompanying exact example has reversed instantaneous birth rates while its true fixation probabilities remain correctly ordered.

## Use the exact oracle

No third-party Python dependency is required.

~~~sh
python count_moran.py --parts 2,3 --fitness 2 --bd 1/2 --initial 1,0
python verify.py
python -O verify.py
~~~

The oracle reduces binary vertex configurations to the product of (part_size+1) count states. This can be much smaller for nontrivial parts. It uses exact rational arithmetic and a dense harmonic solve, with a default 256-count-state resource guard; this finite implementation limit does not restrict the analytic theorem. Larger guards should be chosen deliberately.

See proof.md for the complete argument and sources.md for the native contract and prior comparison. The verifier independently enumerates parent-to-target events in the fine binary-state model, checks exact lumpability, compares fine/coarse harmonic solutions, and checks finite order/witness cases. It is read-only; normal and optimized Python retain the same explicit checks. Ordinary pre-existing bytecode caches are tolerated, but verification executes current source directly and does not rewrite the packet. These finite checks do not mechanically certify the infinite theorem.

## Status and credit

This is a modest proposed proof/code contribution. The historical broad monotonicity belief and its present universal status remain unresolved/uncertain here. The later journal manuscript is distinct from the ITCS manuscript. Pure-rule comparisons, earlier bipartite count reduction, and generic coupling machinery are credited; no new empirical biology, fitness paradox, first-proof or major-discovery claim is made.

Substantial AI assistance was used for source discovery, derivation, implementation and separate fresh-context checking. No source-author review, external human peer review, proof-assistant certificate or priority confirmation has been obtained.

All prose/code in this directory is original and MIT-licensed (LICENSE). No author source code, paper PDF, private correspondence or internal workflow material is redistributed. The original sources retain their own licenses.
