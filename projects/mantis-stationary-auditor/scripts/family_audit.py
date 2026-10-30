import json
import sys
from fractions import Fraction
import numpy as np
from mantis import CRNetwork
from mantis.parsing import parse_reactions
from mantis.stochastic_stationary import StationaryDistribution,_reachable_finite_class
from exact_oracle import stationary,require
NA=6.02214076e23
records=[];max_error=0.0;float_diagnostics=[]
for d in range(1,5):
 for forward,backward in [(1,1),(2,3)]:
  expression=f"{d}A <-> {d}B" if d>1 else "A <-> B"
  reactions=parse_reactions([expression]);rates={reactions[0].rate_key:forward,reactions[1].rate_key:backward}
  for total in range(1,9):
   states=np.array([(a,total-a) for a in range(total+1)],dtype=np.int64)
   ratio=(backward/forward)**(1/d);means=np.array([total*ratio/(1+ratio),total/(1+ratio)])
   for a0 in range(total+1):
    initial=(a0,total-a0);qreactions=[((d,0),(0,d),forward),((0,d),(d,0),backward)]
    exact=stationary(initial,qreactions)
    support=_reachable_finite_class(states,np.array(initial),reactions,["A","B"],rates,{})
    require(set(map(tuple,support))==set(exact),"integer reachable support mismatch")
    law=StationaryDistribution(["A","B"],means,1/NA,True,True,_states=support)
    error=max(abs(law.probability({"A":a,"B":b})-float(exact.get((a,b),0))) for a,b in states)
    require(error<1e-10,"finite CTMC/product-form probability mismatch")
    require(law.probability({"A":0,"B":0})==0,"outside total support")
    max_error=max(max_error,error)
    records.append({"step":d,"total":total,"initial":list(initial),"rates":[forward,backward],"support":[list(x) for x in exact],"exact_pi":[str(v) for v in exact.values()],"float_comparison_tolerance":1e-10,"float_comparison_passed":True})
    float_diagnostics.append(error)
# End-to-end source API controls, including an under-populated catalyst.
end_to_end=[]
for expr,rates,initial,qreact in [
 ("A <-> B",{"A -> B":1.,"B -> A":1.},(2,0),[((1,0),(0,1),1),((0,1),(1,0),1)]),
 ("2A <-> 2B",{"2A -> 2B":1.,"2B -> 2A":1.},(2,0),[((2,0),(0,2),1),((0,2),(2,0),1)]),
 ("2A <-> 2B",{"2A -> 2B":2.,"2B -> 2A":3.},(1,1),[((2,0),(0,2),2),((0,2),(2,0),3)]),
 ("3A <-> 3B",{"3A -> 3B":2.,"3B -> 3A":3.},(8,0),[((3,0),(0,3),2),((0,3),(3,0),3)]),
 ("A + 2C <-> B + 2C",{"A + 2C -> B + 2C":1.,"B + 2C -> A + 2C":1.},(2,0,1),[((1,0,2),(0,1,2),1),((0,1,2),(1,0,2),1)])]:
 names=["A","B"] if len(initial)==2 else ["A","B","C"];net=CRNetwork.from_string([expr],rates=rates)
 law=net.stationary_distribution(dict(zip(names,initial)),volume_L=1/NA,initial_as="count");exact=stationary(initial,qreact)
 require(set(map(tuple,law._states))==set(exact),"end-to-end support mismatch")
 error=max(abs(law.probability(dict(zip(names,x)))-float(v)) for x,v in exact.items())
 require(error<1e-7,"end-to-end numeric equilibrium tolerance")
 end_to_end.append({"reaction":expr,"initial":initial,"support":[list(x) for x in exact],"exact_pi":[str(v) for v in exact.values()],"float_comparison_tolerance":1e-7,"float_comparison_passed":True})
 float_diagnostics.append(error)
# Absent catalyst: exact support helper only; initial-class-positive-equilibrium restriction unchanged.
rx=parse_reactions(["2A + C <-> 2B + C"]);rv={r.rate_key:1 for r in rx};allowed=np.array([(0,2,0),(1,1,0),(2,0,0)])
support=_reachable_finite_class(allowed,np.array([2,0,0]),rx,["A","B","C"],rv,{})
require(set(map(tuple,support))=={(2,0,0)},"absent catalyst control")
legacy_q=[[Fraction(-2),Fraction(0),Fraction(2)],[Fraction(0),Fraction(0),Fraction(0)],[Fraction(2),Fraction(0),Fraction(-2)]]
legacy_pi=[Fraction(1,4),Fraction(1,2),Fraction(1,4)]
require(all(sum(legacy_pi[i]*legacy_q[i][j] for i in range(3))==0 for j in range(3)),"legacy mixture is stationary negative control")
legacy_control={"states":[[0,2],[1,1],[2,0]],"legacy_pi":list(map(str,legacy_pi)),"legacy_full_plane_stationarity_residuals":["0","0","0"],"legacy_variance_A":"1/2","initial_component_variance_A":"1","lesson":"Stationarity alone cannot identify the initial-state-conditioned law."}
result={"ok":True,"oracle":"independent integer BFS and exact Fraction CTMC stationary equations","complete_bounded_family":{"steps":[1,2,3,4],"totals":"1..8","initial_a":"every0..total","rate_pairs":[[1,1],[2,3]],"cases":len(records)},"helper_and_product_law_cases":records,"end_to_end_cases":end_to_end,"legacy_stationarity_negative_control":legacy_control,"helper_only_absent_catalyst":True,"helper_float_comparison_tolerance":1e-10,"end_to_end_float_comparison_tolerance":1e-7,"float_comparisons_all_passed":True,"scope":"finite complete conditioned class, valid integer queries; no infinite-component or new ACK theorem claim"}
if "--diagnostics" in sys.argv:
 result["observed_float_errors"] = float_diagnostics
 result["observed_max_helper_error"] = max_error
print(json.dumps(result,indent=2,sort_keys=True))
