"""Portable independent exact verifier; Python standard library only.

Generation starts with every 6-pair Dyck word and inserts binary dot gaps,
rather than generating contracted helix forests. Folding DP decomposes at the
rightmost nucleotide. Does not import either author's implementation, write
files, or read private audit/snapshot files. Successful stdout is one JSON value.
Every validation uses explicit exceptions and remains active under python -O.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVERSE = (1, 0, 2)
COMPLEMENT = {'A':'U', 'U':'A', 'G':'C', 'C':'G'}

class VerificationError(ValueError):
    """An exact check failed; never suppressed by Python optimization."""

def require(condition, message):
    if not condition:
        raise VerificationError(message)

def tree(s):
    require(isinstance(s,str), 'Target must be a dot-bracket string')
    pairs, parent, kids, dots, stack = [], [], [[]], [[]], [0]
    for i, ch in enumerate(s):
        if ch == '(':
            node = len(kids)
            parent.append(stack[-1])
            kids[stack[-1]].append(node)
            kids.append([]); dots.append([]); pairs.append([i, None])
            stack.append(node)
        elif ch == ')':
            require(len(stack)>1, 'Unmatched closing bracket')
            node = stack.pop()
            pairs[node-1][1] = i
        elif ch == '.':
            dots[stack[-1]].append(i)
        else:
            raise VerificationError('Invalid target character: '+repr(ch))
    require(stack == [0], 'Unmatched opening bracket')
    return pairs, parent, kids, dots

@lru_cache(None)
def dyck(n):
    if n == 0:
        return ('',)
    return tuple('('+a+')'+b for k in range(n) for a in dyck(k) for b in dyck(n-1-k))

def helix_lengths(kids, dots):
    # Consecutive nesting is a unary paired vertex with no singleton child.
    continues = {u: kids[u][0] for u in range(1,len(kids)) if len(kids[u])==1 and not dots[u]}
    heads = set(range(1,len(kids))) - set(continues.values())
    lengths = []
    for h in heads:
        length = 1
        while h in continues:
            h = continues[h]; length += 1
        lengths.append(length)
    return sorted(lengths)

def generate():
    targets = set()
    for skeleton in dyck(6):
        _, _, kids, _ = tree(skeleton)
        if len(kids[0]) > 4 or any(len(k)>3 for k in kids[1:]):
            continue
        # Owner of each insertion gap, before/after the skeleton's symbols.
        owners, stack = [0], [0]
        next_node = 1
        for ch in skeleton:
            if ch == '(':
                stack.append(next_node); next_node += 1
            else:
                stack.pop()
            owners.append(stack[-1])
        eligible = [i for i,owner in enumerate(owners) if len(kids[owner]) <= (2 if owner==0 else 1)]
        for selected in product((False,True), repeat=len(eligible)):
            occupied = {i for i,bit in zip(eligible,selected) if bit}
            text = ''.join(('.' if i in occupied else '') + (skeleton[i] if i<len(skeleton) else '') for i in range(len(skeleton)+1))
            _, _, children, singletons = tree(text)
            if helix_lengths(children,singletons) == [2,2,2]:
                targets.add(text)
    return targets

def count_colorings(s):
    pairs, parents, kids, dots = tree(s)
    counts = {'proper':0, 'mod2':0, 'uniform':0, 'mod2_vacuous':0}
    for colors in product(range(3), repeat=len(pairs)):
        colors = (None,) + colors
        if any(any(n>cap for n,cap in zip((Counter(([INVERSE[colors[u]]] if u else []) + [colors[v] for v in kids[u]])[c] for c in range(3)),(1,1,2))) for u in range(len(kids))):
            continue
        counts['proper'] += 1
        level = [0] * len(kids)
        for v,p in enumerate(parents,1):
            level[v] = level[p] + (1,-1,0)[colors[v]]
        gray = {level[v]%2 for v in range(1,len(kids)) if colors[v]==2}
        unpaired = {level[u]%2 for u in range(len(kids)) if dots[u]}
        if not gray & unpaired:
            counts['mod2'] += 1
            counts['mod2_vacuous'] += int(not gray or not unpaired)
        if any(unpaired <= {residue} and gray <= {1-residue} for residue in (0,1)):
            counts['uniform'] += 1
    return counts

def fold(w):
    """Exact optimal-score/count DP using the last position's unique partner."""
    n = len(w)
    score, ways = {}, {}
    def get(i,j):
        return (0,1) if i>j else (score[i,j],ways[i,j])
    for width in range(1,n+1):
        for i in range(n-width+1):
            j = i+width-1
            best,count = get(i,j-1)
            for k in range(i,j):
                if COMPLEMENT[w[k]] != w[j]:
                    continue
                left,cl = get(i,k-1); inside,ci = get(k+1,j-1)
                trial = left+inside+1
                if trial>best:
                    best,count = trial,cl*ci
                elif trial==best:
                    count += cl*ci
            score[i,j], ways[i,j] = best,count
    return get(0,n-1)

def check_witness(s,w,expected_score):
    pairs,_,_,_ = tree(s)
    require(isinstance(w,str) and set(w)<=set(COMPLEMENT), 'Invalid RNA sequence')
    require(len(s)==len(w), 'Target/sequence lengths differ')
    require(all(COMPLEMENT[w[i]]==w[j] for i,j in pairs), 'Target is incompatible with supplied sequence')
    answer = fold(w)
    require(answer == (expected_score,1), 'Nonunique or incorrect optimum: '+repr((s,answer)))
    return answer

def short_intervals(pairs):
    intervals=set(map(tuple,pairs))
    shorts=[]
    for i,j in pairs:
        if (i-1,j+1) in intervals:
            continue
        length=1
        while (i+length,j-length) in intervals:
            length+=1
        if length==2:
            shorts.append((i,j))
    return shorts

def comparable(shorts):
    return all(a==b or a[0]<b[0]<b[1]<a[1] or b[0]<a[0]<a[1]<b[1] for a in shorts for b in shorts)

def check_coloring_witness(s,provided):
    pairs,parents,kids,dots=tree(s)
    require(isinstance(provided,list) and len(provided)==len(pairs) and all(type(c) is int and c in range(3) for c in provided), 'Invalid serialized coloring')
    colors=[None]+provided
    for u in range(len(kids)):
        histogram=Counter(([INVERSE[colors[u]]] if u else [])+[colors[v] for v in kids[u]])
        require(histogram[0]<=1 and histogram[1]<=1 and histogram[2]<=2, 'Serialized coloring is not proper')
    levels=[0]*len(kids)
    for v,p in enumerate(parents,1):
        levels[v]=levels[p]+(1,-1,0)[colors[v]]
    gray={levels[v]%2 for v in range(1,len(kids)) if colors[v]==2}
    unpaired={levels[u]%2 for u in range(len(kids)) if dots[u]}
    require(any(unpaired<={r} and gray<={1-r} for r in (0,1)), 'Serialized coloring is not uniformly modulo-two separated')

def check_local_states():
    # Explicit states (parent parity, head color), retaining B/W separately.
    states = ((0,0),(0,1),(1,0),(1,1),(1,2))
    X,Y,Z = set(states[:2]),set(states[2:4]),set(states[4:])
    F,Q = X|Y|Z,X|Z
    def permitted(c,cs):
        hist = Counter([INVERSE[c]]+list(cs))
        return hist[0]<=1 and hist[1]<=1 and hist[2]<=2
    def terminal(children,has_dots):
        result = set()
        for p,c in states:
            q = (p+int(c!=2))%2
            if has_dots and q!=0:
                continue
            for cs in product(range(3),repeat=len(children)):
                if permitted(c,cs) and all((q,d) in child for d,child in zip(cs,children)):
                    result.add((p,c));break
        return result
    def outside(S):
        return {(p,c) for p,c in states if any(permitted(c,(d,)) and (((p+int(c!=2))%2),d) in S for d in range(3))}
    labels = [set(),X,Y,X|Y,Z,X|Z,Y|Z,F]
    table = [sum(bit for bit,part in ((1,X),(2,Y),(4,Z)) if part<=outside(mask)) for mask in labels]
    require(table == [0,2,5,7,5,7,5,7], 'Incorrect local transfer table')
    terminal_rows = []
    for degree,dotted,expected in ((0,False,F),(0,True,Y),(1,True,Y),(2,False,Q),(3,False,Q)):
        for exceptional in range(-1,degree):
            children=[Q if j==exceptional else F for j in range(degree)]
            observed=terminal(children,dotted)
            require(observed==expected, 'Incorrect terminal state set')
            require(outside(observed)>=Q, 'Length-two guarantee failed')
            require(outside(outside(observed))==F, 'Length-three guarantee failed')
            terminal_rows.append({'degree':degree,'dots':dotted,'Q_child':exceptional,'states':sorted(observed)})
    # Critical fork with TWO Q children, which is excluded from theorem 1's lemma.
    require(terminal([Q,Q],False)==X, 'Two-Q fork terminal differs from X')
    require(outside(X)==Y, 'Two-Q fork outer helix differs from Y')
    return {'D_table':table,'terminal_rows':terminal_rows,'two_Q_fork_terminal':sorted(X),'two_Q_fork_h2':sorted(Y)}

def check_held_out_shape(row):
    pairs,_,kids,dots = tree(row['target'])
    require(len(kids[0])<=4 and all(len(k)<=3 for k in kids[1:]), 'Held-out m5 motif')
    require(not dots[0] or len(kids[0])<=2, 'Held-out root m3dot motif')
    require(all(not dots[u] or len(kids[u])<=1 for u in range(1,len(kids))), 'Held-out nonroot m3dot motif')
    lengths=helix_lengths(kids,dots)
    require(1 not in lengths, 'Held-out length-one helix')
    require(lengths.count(2)==row['short_helices'], 'Held-out short-helix count differs')
    require(comparable(short_intervals(pairs)), 'Held-out short helices are not ancestry-comparable')
    require(row['length']==len(row['target']), 'Incorrect held-out length metadata')
    return len(pairs)

def verify(input_path):
    raw=input_path.read_bytes()
    data=json.loads(raw)
    require(isinstance(data,dict) and isinstance(data.get('classification'),dict) and isinstance(data.get('held_out'),list), 'Expected classification object and held_out list')
    classification=data['classification']
    supplied=classification.get('results')
    require(isinstance(supplied,list) and all(isinstance(row,dict) and isinstance(row.get('target'),str) for row in supplied), 'Expected classification.results target records')
    generated=generate()
    by_target={row['target']:row for row in supplied}
    require(len(by_target)==len(supplied), 'Duplicate target records')
    require(len(generated)==288, 'Independent generation did not produce 288 classes')
    require(generated==set(by_target), 'Supplied target domain differs from independent generation')
    rows=[]
    for target in sorted(generated):
        counts=count_colorings(target)
        pairs,_,_,_=tree(target)
        is_nested=comparable(short_intervals(pairs))
        require(counts['uniform']==by_target[target].get('uniform_colorings'), 'Uniform coloring count differs for '+target)
        require(by_target[target].get('colorings_checked')==729, 'Incorrect coloring-domain metadata')
        require(by_target[target].get('shorts_nested') is is_nested, 'Incorrect ancestry-comparability metadata')
        require(bool(counts['mod2'])==bool(counts['uniform']), 'Modulo-two/uniform existence differs')
        witness=None
        if counts['uniform']:
            witness=check_witness(target,by_target[target]['sequence'],6)
            check_coloring_witness(target,by_target[target].get('colors'))
            require(by_target[target].get('optimum_and_count')==list(witness), 'Incorrect classification optimum metadata')
        rows.append({'target':target,**counts,'shorts_nested':is_nested,'independent_optimum_and_count':witness})
    held=[]
    require([row.get('short_helices') for row in data['held_out']]==[3,4,5,7,9,12], 'Unexpected held-out family')
    for row in data['held_out']:
        target_pairs=check_held_out_shape(row)
        optimum=check_witness(row['target'],row['sequence'],target_pairs)
        check_coloring_witness(row['target'],row.get('coloring'))
        require((row.get('optimum'),row.get('optimum_count'))==optimum, 'Incorrect held-out optimum metadata')
        held.append({'short_helices':row['short_helices'],'length':len(row['target']),'score_count':optimum})
    # Negative controls: alternate compatible perfect folds, incompatible target.
    require(fold('AUAU')==(2,2), 'Nonunique-fold negative control failed')
    require(fold('AA')==(0,1), 'Incompatible-pair negative control failed')
    require(fold('')==(0,1), 'Empty-sequence boundary control failed')
    failures=[row['target'] for row in rows if not row['mod2']]
    require(set(failures)=={'((((.))((.)))).','.((((.))((.))))','.((((.))((.)))).','((.))((.))((.))'}, 'Failure classes differ from claimed geometries')
    positives=sum(bool(r['uniform']) for r in rows)
    require(classification.get('targets')==len(rows) and classification.get('nonseparable_uniform')==len(failures) and classification.get('nussinov_positive_witnesses')==positives, 'Incorrect classification summary metadata')
    nested=sum(r['shorts_nested'] for r in rows)
    require(classification.get('nested_targets')==nested==72, 'Incorrect ancestry-comparable class count')
    states=check_local_states()
    return {'schema':'independent-exact-check-v1','input_sha256':sha256(raw).hexdigest(),'generation':'all 132 six-pair Dyck skeletons with independently inserted motif-legal binary dot gaps','normalized_targets':len(rows),'nested_targets':nested,'colorings_checked_per_target':729,'mod2_failures':failures,'unique_positive_witnesses':positives,'ordinary_mod2_coloring_witnesses':sum(r['mod2'] for r in rows),'uniform_coloring_witnesses':sum(r['uniform'] for r in rows),'targets_with_vacuity_count_difference':sum(r['mod2']!=r['uniform'] for r in rows),'local_transfer_table':states['D_table'],'local_terminal_rows_checked':len(states['terminal_rows']),'two_Q_fork_checked':True,'held_out_checks':held,'negative_controls':{'AUAU':[2,2],'AA':[0,1],'empty':[0,1]}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=HERE/'evidence'/'designs.json',help='Public designs receipt (default: beside this script in evidence/designs.json)')
    args=parser.parse_args()
    print(json.dumps(verify(args.input),sort_keys=True,separators=(',',':')))

if __name__=='__main__':
    main()
