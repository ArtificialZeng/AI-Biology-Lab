"""Independent direct dot-bracket checks; does not import the state algebra.
For the finite 3-stack classification, brute-force all 3^6 pair colorings.
For each positive witness, construct a sequence and count all optimal folds.
"""
from itertools import product
from functools import lru_cache
import json

@lru_cache(None)
def forests(n,d):
    if d==0:return ('',) if n==0 else ()
    out=[]
    for k in range(1,n+1):
        for a in helices(k):
            for b in forests(n-k,d-1):out.append(a+b)
    return tuple(out)
@lru_cache(None)
def helices(n):
    out=[]
    for d in range(min(3,n-1)+1):
        for inside in forests(n-1,d):
            # split the forest into d top-level child strings
            parts=split(inside)
            gaps=[()] if d>=2 else list(product(range(2),repeat=d+1))
            for bits in gaps:
                if d==1 and not any(bits):continue
                if d>=2:content=inside
                else:content=''.join(('.' if bits[i] else '')+(parts[i] if i<d else '') for i in range(d+1))
                out.append('(('+content+'))')
    return tuple(out)
def split(s):
    parts=[]; level=0; start=0
    for i,c in enumerate(s):
        level+=(c=='(')-(c==')')
        if level==0:parts.append(s[start:i+1]);start=i+1
    return parts

def parse(s):
    pairs=[]; children={-1:[]}; unpaired={-1:[]}; stack=[]
    for i,c in enumerate(s):
        if c=='(':
            k=len(pairs);pairs.append([i,None]);children[k]=[];unpaired[k]=[]
            children[stack[-1] if stack else -1].append(k);stack.append(k)
        elif c==')':pairs[stack.pop()][1]=i
        else:unpaired[stack[-1] if stack else -1].append(i)
    return pairs,children,unpaired

def brute(s):
    pairs,ch,up=parse(s);n=len(pairs);good=[]
    for colors in product(range(3),repeat=n):
        valid=True
        for p,kids in ch.items():
            cs=[colors[k] for k in kids]
            if cs.count(0)>1 or cs.count(1)>1 or cs.count(2)>2:valid=False;break
            if p>=0:
                if cs.count(colors[p])>1:valid=False;break
                if colors[p]==0 and 1 in cs or colors[p]==1 and 0 in cs:valid=False;break
        if not valid:continue
        levels={-1:0}
        for p in range(n):
            par=next(k for k,kids in ch.items() if p in kids)
            levels[p]=levels[par]+(1 if colors[p]==0 else -1 if colors[p]==1 else 0)
        ul={levels[p]%2 for p in up if up[p]}
        gl={levels[p]%2 for p in range(n) if colors[p]==2}
        if any(all(v==xi for v in ul) and all(v!=xi for v in gl) for xi in range(2)):
            good.append(colors)
    return good

def sequence(s,cols):
    pairs,ch,up=parse(s);w=['A']*len(s);right={-1:None}
    for p in [-1]+list(range(len(pairs))):
        used=set()
        for k in ch[p]:
            c=cols[k]
            if c==0:a,b='G','C'
            elif c==1:a,b='C','G'
            else:
                a=next(x for x in 'AU' if x not in used and x!=right[p]);b='U' if a=='A' else 'A'
            used.add(a);right[k]=b;i,j=pairs[k];w[i]=a;w[j]=b
    return ''.join(w)

def nussinov(w):
    comp={'A':'U','U':'A','G':'C','C':'G'}
    @lru_cache(None)
    def f(i,j):
        if i>j:return 0,1
        best,count=f(i+1,j)
        for k in range(i+1,j+1):
            if w[k]!=comp[w[i]]:continue
            a,ca=f(i+1,k-1);b,cb=f(k+1,j);v=1+a+b
            if v>best:best,count=v,ca*cb
            elif v==best:count+=ca*cb
        return best,count
    return f(0,len(w)-1)

def nested(s):
    pairs,ch,_=parse(s)
    heads=[p for p,(i,j) in enumerate(pairs) if not any(a==i-1 and b==j+1 for a,b in pairs)]
    return all((pairs[a][0]<pairs[b][0]<pairs[b][1]<pairs[a][1]) or (pairs[b][0]<pairs[a][0]<pairs[a][1]<pairs[b][1]) for a in heads for b in heads if a!=b)

def all_targets3():
    out=set()
    for d in range(1,4):
        for forest in forests(3,d):
            parts=split(forest)
            for gaps in ([()] if d>2 else product(range(2),repeat=d+1)):
                if d>2:s=forest
                else:s=''.join(('.' if gaps[i] else '')+(parts[i] if i<d else '') for i in range(d+1))
                out.add(s)
    return sorted(out,key=lambda s:(len(s),s))


# Separate direct interval-tree recursion for larger held-out examples.
def find_direct(s):
    pairs,ch,up=parse(s)
    @lru_cache(None)
    def node(p,level,c):
        here=(level+(c!=2))%2
        if c==2 and here!=1 or up[p] and here!=0:return None
        for cs in product(range(3),repeat=len(ch[p])):
            if cs.count(0)>1 or cs.count(1)>1 or cs.count(2)>2:continue
            if cs.count(c)>1 or c==0 and 1 in cs or c==1 and 0 in cs:continue
            subs=[node(k,here,d) for k,d in zip(ch[p],cs)]
            if all(z is not None for z in subs):
                out={p:c}
                for z in subs:out.update(z)
                return out
        return None
    for lev in range(2):
        if up[-1] and lev!=0:continue
        for cs in product(range(3),repeat=len(ch[-1])):
            if cs.count(0)>1 or cs.count(1)>1 or cs.count(2)>2:continue
            subs=[node(k,lev,c) for k,c in zip(ch[-1],cs)]
            if all(z is not None for z in subs):
                out={}
                for z in subs:out.update(z)
                return tuple(out[k] for k in range(len(pairs)))
    return None

def short_heads(s):
    ps,_,_=parse(s);st=set(map(tuple,ps));hs=[]
    for i,j in ps:
        if (i-1,j+1) in st:continue
        h=1
        while (i+h,j-h) in st:h+=1
        if h==2:hs.append((i,j))
    return hs


def run():
    results=[]
    for s in all_targets3():
        good=brute(s);r={'target':s,'colorings_checked':3**6,'uniform_colorings':len(good),'shorts_nested':nested(s)}
        if good:
            w=sequence(s,good[0]);r.update(sequence=w,colors=list(good[0]),optimum_and_count=list(nussinov(w)))
            if not (nussinov(w)==(6,1)): raise ValueError('Exact RNA certificate check failed')
        if not (good or not nested(s)): raise ValueError('Exact RNA certificate check failed')
        results.append(r)
    summary={'targets':len(results),'nested_targets':sum(r['shorts_nested'] for r in results),'nonseparable_uniform':sum(r['uniform_colorings']==0 for r in results),'nussinov_positive_witnesses':sum(r['uniform_colorings']>0 for r in results),'results':results}
    held=[]
    for k in [3,4,5,7,9,12]:
        t='((...))'
        for a in range(k-1):
            side='('* (3+a%2)+'...' +')'*(3+a%2)
            t='(('+ (t+side if a%2 else side+t) +'))'
        t='.'+'(((.'+t+'.)))' +'.'
        sh=short_heads(t)
        if not (len(sh)==k and all(a==b or a[0]<b[0]<b[1]<a[1] or b[0]<a[0]<a[1]<b[1] for a in sh for b in sh)): raise ValueError('Exact RNA certificate check failed')
        cols=find_direct(t)
        if not (cols is not None): raise ValueError('Exact RNA certificate check failed')
        w=sequence(t,cols);best,count=nussinov(w)
        if not ((best,count)==(t.count('('),1)): raise ValueError('Exact RNA certificate check failed')
        held.append({'target':t,'short_helices':k,'sequence':w,'optimum':best,'optimum_count':count,'length':len(t),'coloring':list(cols)})
    return {'classification':summary,'held_out':held}

if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
