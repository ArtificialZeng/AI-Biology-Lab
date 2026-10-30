"""Exact finite-state algebra for proper uniform modulo-2 RNA colorings.
Relative residues: unpaired=0, gray=1. Bits X,Y,Z mean a non-gray
head at entry 0, a non-gray head at entry 1, and a gray head at entry 1.
No sequence energetics or real-valued computation is used.
"""
from itertools import product
import json
B,W,G=0,1,2
INV={B:W,W:B,G:G}
F,Q=7,5

def flag(p,c):
    return (1 if p==0 else 2) if c!=G else (4 if p==1 else 0)
def proper(parent,children):
    xs=([INV[parent]] if parent is not None else [])+list(children)
    return xs.count(B)<=1 and xs.count(W)<=1 and xs.count(G)<=2

def step(mask):
    out=0
    for p,c in product(range(2),range(3)):
        q=p^(c!=G)
        if c==G and q!=1: continue
        if any(proper(c,[d]) and mask&flag(q,d) for d in range(3)):
            out|=flag(p,c)
    return out

def terminal(child_masks,dots):
    out=0
    for p,c in product(range(2),range(3)):
        q=p^(c!=G)
        if c==G and q!=1: continue
        if dots and q!=0: continue
        for cs in product(range(3),repeat=len(child_masks)):
            if proper(c,cs) and all(m&flag(q,d) for m,d in zip(child_masks,cs)):
                out|=flag(p,c)
                break
    return out

def helix(h,children,dots):
    mask=terminal(children,dots)
    for _ in range(h-1): mask=step(mask)
    return mask

def root(children,dots):
    feasible=[]
    for p in range(2):
        if dots and p!=0: continue
        for cs in product(range(3),repeat=len(children)):
            if proper(None,cs) and all(m&flag(p,c) for m,c in zip(children,cs)):
                feasible.append([p,list(cs)])
    return feasible

def cases():
    rows=[]
    for h in [2,3,4,5,6,7,8]:
        for d in range(4):
            for dots in [False,True]:
                if dots and d>1: continue
                if not dots and d==1: continue # then chains merge; not a maximal-helix terminal
                for qchild in ([None] if d==0 else [None]+list(range(d))):
                    children=[Q if j==qchild else F for j in range(d)]
                    rows.append(dict(h=h,d=d,dots=dots,qchild=qchild,mask=helix(h,children,dots)))
    return rows

def run():
    rows=cases()
    checks={
      'all_long_support_F':all(r['mask']==F for r in rows if r['h']>=3 and r['qchild'] is None),
      'short_or_one_demand_support_Q':all(r['mask']&Q==Q for r in rows),
      'one_demand_root_exists':all(root([Q if j==qchild else F for j in range(d)],dots)
         for d in range(5) for dots in [False,True] if not dots or d<=2
         for qchild in ([None] if d==0 else [None]+list(range(d)))),
      'all_long_h_ge3_stable_in_admitted_cases':all(r['mask']==F for r in rows if r['h']>=3) and step(F)==F,
      'arbitrary_mask_stability_after_three_edges':all(step(step(step(step(m))))==step(step(step(m))) for m in range(8)),
    }
    # For the admitted cases, h=3 reaches F and step(F)=F,
    # proving the result for every h>=3, not merely h<=8.
    result={'checks':checks,'rows':rows,'step_table':{m:step(m) for m in range(8)}}
    if not all(checks.values()): raise ValueError('A finite-state invariant failed')
    return result

if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
