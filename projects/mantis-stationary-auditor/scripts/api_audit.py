"""Actual pinned CRNetwork API comparisons, with deterministic boolean output."""
import json,sys
from mantis import CRNetwork
from exact_oracle import require
NA=6.02214076e23
variant=sys.argv[1]
rows=[]
for reaction,rates in [("A <-> B",{"A -> B":1.,"B -> A":1.}),("2A <-> 2B",{"2A -> 2B":1.,"2B -> 2A":1.})]:
 law=CRNetwork.from_string([reaction],rates=rates).stationary_distribution({"A":2,"B":0},volume_L=1/NA,initial_as="count")
 values,probs=law.marginal("A");mean=float(values@probs);variance=float((values-mean)**2@probs)
 is_double=reaction.startswith("2")
 expected_states={(0,2),(2,0)} if variant=="candidate" and is_double else {(0,2),(1,1),(2,0)}
 expected_outside=0 if variant=="candidate" else .5
 expected_middle=0 if variant=="candidate" and is_double else .5
 expected_var=1 if variant=="candidate" and is_double else .5
 require(set(map(tuple,law._states))==expected_states,"unexpected actual API support")
 require(abs(law.probability({"A":0,"B":0})-expected_outside)<1e-7,"outside probability")
 require(abs(law.probability({"A":1,"B":1})-expected_middle)<1e-7,"middle probability")
 require(abs(variance-expected_var)<1e-7,"variance")
 if variant=="candidate":require(law.probability({"A":0,"B":0})==0,"outside must be exactly zero")
 rows.append({"reaction":reaction,"expected_support":[list(v) for v in sorted(expected_states)],"expected_variance_A":"1" if expected_var==1 else "1/2","expected_probability_00":"0" if expected_outside==0 else "1/2","expected_probability_11":"0" if expected_middle==0 else "1/2","actual_library_comparisons_passed":True,"tolerance":1e-7})
print(json.dumps({"variant":variant,"cases":rows,"legacy_law_is_stationary_mixture":True},indent=2,sort_keys=True))
