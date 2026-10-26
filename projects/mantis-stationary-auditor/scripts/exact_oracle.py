"""Independent finite CTMC oracle: integer enabling and Fraction linear algebra."""
from collections import deque
from fractions import Fraction as F
from math import prod

def require(ok,message):
    if not ok:raise RuntimeError(message)

def solve(matrix,rhs):
    n=len(rhs);rows=[[F(x) for x in row]+[F(v)] for row,v in zip(matrix,rhs)]
    for col in range(n):
        pivot=next((r for r in range(col,n) if rows[r][col]),None)
        require(pivot is not None,"singular finite CTMC normalization")
        rows[col],rows[pivot]=rows[pivot],rows[col]
        scale=rows[col][col];rows[col]=[x/scale for x in rows[col]]
        for r in range(n):
            if r!=col and rows[r][col]:
                a=rows[r][col];rows[r]=[x-a*y for x,y in zip(rows[r],rows[col])]
    return [row[-1] for row in rows]

def jumps(state,reactions):
    out={}
    for reactant,product,rate in reactions:
        if all(x>=a for x,a in zip(state,reactant)):
            propensity=F(rate)*prod(prod(range(x-a+1,x+1)) for x,a in zip(state,reactant))
            after=tuple(x-a+b for x,a,b in zip(state,reactant,product))
            if propensity and after!=state:out[after]=out.get(after,F(0))+propensity
    return out

def stationary(initial,reactions,budget=1000):
    initial=tuple(initial);seen={initial};queue=deque([initial])
    while queue:
        state=queue.popleft()
        for after in jumps(state,reactions):
            if after not in seen:
                require(len(seen)<budget,"finite oracle state budget")
                seen.add(after);queue.append(after)
    states=sorted(seen);index={x:i for i,x in enumerate(states)};q=[[F(0) for _ in states] for _ in states]
    for x in states:
        for y,rate in jumps(x,reactions).items():
            q[index[x]][index[y]]+=rate;q[index[x]][index[x]]-=rate
    matrix=[[q[j][i] for j in range(len(states))] for i in range(len(states))]
    matrix[-1]=[F(1)]*len(states);rhs=[F(0)]*(len(states)-1)+[F(1)];pi=solve(matrix,rhs)
    require(sum(pi)==1 and min(pi)>0,"positive normalized finite law")
    require(all(sum(pi[i]*q[i][j] for i in range(len(states)))==0 for j in range(len(states))),"exact CTMC stationary residual")
    return dict(zip(states,pi))
