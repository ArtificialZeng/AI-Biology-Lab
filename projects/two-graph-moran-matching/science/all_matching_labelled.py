"""Separately implemented full-labelled parent-target event check.
Does not import certificate.py or any R14/reduced-chain checker.
"""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
cert=json.loads((ROOT/'evidence'/'uniform-certificate.json').read_text())

def orbit(mask,n,m):
 a=b=0
 for pair in range(m):
  v=((mask>>(2*pair))&1)+((mask>>(2*pair+1))&1)
  a+=(v==2); b+=(v==1)
 c=sum((mask>>i)&1 for i in range(2*m,n))
 return a,b,m-a-b,c,n-2*m-c

def w_scaled(mask,n,m):
 # Exactly n(n-1)W; compute directly from labelled configuration.
 count=mask.bit_count()
 full=sum(((mask>>(2*i))&3)==3 for i in range(m))
 return ((n-2)*n+m-2)*(n-1)*count+n*(n-1)*full+m*count*(n-count)

def p_value(args):
 return sum(t['coefficient']*__import__('math').prod(x**j for x,j in zip(args,t['powers'])) for t in cert['P_terms'])

cases=[];configs=0;events=0
for n,m in [(3,1),(4,1),(4,2),(5,2),(6,3),(7,3),(8,4),(9,3),(10,5),(11,4),(11,5)]:
 nconfigs=0;nevents=0;zeros=0;singletons=[]
 for mask in range(1<<n):
  old=w_scaled(mask,n,m); scaled_q=0
  for parent in range(n):
   mutant=(mask>>parent)&1
   paired_mutant=mutant and parent<2*m
   degree=n-2 if paired_mutant else n-1
   multiplier=(n-1) if paired_mutant else (n-2)
   neighbors=0
   for target in range(n):
    if target==parent:continue
    if paired_mutant and target==(parent^1):continue
    neighbors+=1
    new=(mask|(1<<target)) if mutant else (mask&~(1<<target))
    scaled_q+=(w_scaled(new,n,m)-old)*multiplier
    nevents+=1
   assert neighbors==degree
  args=orbit(mask,n,m)
  expected=p_value(args)
  assert -scaled_q==expected,(n,m,mask,args,scaled_q,expected)
  assert scaled_q<=0
  zeros+=(scaled_q==0)
  if mask.bit_count()==1:
   assert scaled_q==-2*m*(n-2)*(n-1)
   singletons.append(mask.bit_length()-1)
  nconfigs+=1
 cases.append({'N':n,'m':m,'labelled_configurations':nconfigs,'parent_target_events':nevents,'zero_drift_configurations':zeros,'singletons_checked':len(singletons)})
 configs+=nconfigs;events+=nevents
out={'status':'exact integer labelled-event checks passed','independence':'separate implementation; no import of quotient rates, certificate generator or R14 checkers; checks the supplied symbolic polynomial against direct labelled parent-target sums','scaled_potential':'n(n-1)W, integer','scaled_generator':'D QW where D=n(n-2)(n-1)^2, integer','labelled_configurations':configs,'parent_target_events':events,'cases':cases,'limit':'finite supplemental check; uniform proof depends on exact symbolic certificate and general lumpability argument'}
(ROOT/'evidence'/'labelled-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
