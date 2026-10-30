"""Exact algebra for the regular-opponent theorem; no graph enumeration."""
import argparse
import json
from pathlib import Path
import sympy as sp

def checked_zero(expr, name):
    value=sp.factor(expr)
    if value != 0:
        raise ArithmeticError(f'{name}: {value}')
    return name

def certificate():
    N,d,k,e,T=sp.symbols('N d k e T')
    r=N-k;B=d*k-2*e;U=B/d;L=k*r/(N-1)
    Qk=U-L;Qk2=2*k*Qk+U+L;Qe=T/d-2*e*r/(N-1)
    QV=d*((2*N-1)*Qk+Qk2)/(2*(N-1))-Qe
    D=T-2*d*k*B/(N-1)+r*(d*k/(N-1))**2
    checks=[checked_zero(QV+D/d,'generator_sum_of_squares')]
    V=lambda kval,eval: d*kval*(2*N-1+kval)/(2*(N-1))-eval
    C=d*N**2/(N-1);q=N-1-d
    checks += [checked_zero(V(0,0),'extinction'),checked_zero(V(N,N*d/2)-C,'fixation'),checked_zero(V(1,0)-C/N,'singleton_calibration')]
    Dsingle=d*(1-d/(N-1))**2+q*(d/(N-1))**2
    checks += [checked_zero(Dsingle-d*q/(N-1),'singleton_square_sum'),checked_zero((C/N-Dsingle/(d*N))/C-(1/N-q/(d*N**3)),'singleton_gap')]
    checks += [checked_zero(V(N-1,N*d/2-d)/C-(N-1)/N,'nearfix_calibration')]
    EDnext=2*(d*(1-d/(N-1))**2+q*(d/(N-1))**2)/(N*(N-1))
    checks += [checked_zero(EDnext-2*d*q/(N*(N-1)**2),'nearfix_next_square_sum'),checked_zero(EDnext/(d*N*C)-2*q/(d*N**4*(N-1)),'reverse_two_step_gap')]
    altA=N-3*q/2-q/(N-1)
    altV=altA*k+(k*(k-1)/2-e)+q*k*(N-k)/(2*(N-1))
    checks += [checked_zero(altV-V(k,e),'deleted_graph_potential_equivalence')]
    return {'schema':'regular-opponent-symbolic-v1','checks':checks,'identity':'N(P-I)V=-sum_resident(g_j-d*k/(N-1))^2/d','potential':'V=d*k*(2*N-1+k)/(2*(N-1))-e','singleton_deficit':'(N-1-d)/(d*N^3)','reverse_singleton_surplus':'2*(N-1-d)/(d*N^4*(N-1))','status':'exact-symbolic-identities; graph-domain and optional-stopping proof in proof.md'}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path);ns=a.parse_args();result=certificate();s=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if ns.output:
        with ns.output.open('x') as f:f.write(s)
    else:print(s,end='')
