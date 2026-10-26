"""Deterministic exact rational phylogenies, independent of TreeSwift."""
from dataclasses import dataclass, field
from fractions import Fraction
from copy import deepcopy
import json
from pathlib import Path
from independent_oracle import oracle_newick, gamma_from_intervals

@dataclass
class Spec:
    time: Fraction
    children: list = field(default_factory=list)


def ranked_tree(intervals, mode):
    g = [Fraction(v) for v in intervals]
    root = Spec(Fraction(0), [Spec(Fraction(0)),Spec(Fraction(0))])
    active = list(root.children)
    time = Fraction(0)
    for event in range(1,len(g)):
        time += g[event-1]
        which = 0 if mode == 'comb' else (event*7+3) % len(active)
        chosen = active.pop(which)
        chosen.time = time
        chosen.children = [Spec(time),Spec(time)]
        active.extend(chosen.children)
    horizon = sum(g)
    for leaf in active:
        leaf.time = horizon
    return root


def serialize(root, labels='all', stem=Fraction(0), unary=False, scale=Fraction(1)):
    serial = [0]
    def render(n,parent_time=None):
        serial[0] += 1
        label = 'N'+str(serial[0])
        if not n.children:
            label = 'T'+str(serial[0]) if labels != 'none' else ''
        elif labels == 'none' or labels == 'tips' or (labels == 'mixed' and serial[0] % 2):
            label = ''
        elif labels == 'duplicate':
            label = 'X'
        body = '('+','.join(render(c,n.time) for c in n.children)+')' if n.children else ''
        edge = stem if parent_time is None else (n.time-parent_time)*scale
        # All generated denominators are powers of 2, so the decimal float is exact.
        if unary and parent_time is not None:
            return '('+body+label+':'+str(float(edge/2))+')U:'+str(float(edge/2))
        return body+label+':'+str(float(edge))
    result = render(root)
    if unary:
        result = '('+result+')PREFIX:1e100'
    return result+';'


def cases():
    data=[]
    for n in range(3,9):
        patterns = {'uniform':[Fraction(1)]*(n-1),
                    'late':[Fraction(9)]+[Fraction(1)]*(n-2),
                    'early':[Fraction(1)]*(n-2)+[Fraction(9)],
                    'rational':[Fraction(k+1,8) for k in range(n-1)]}
        # Exact zero gamma with strictly positive integer intervals: solve for
        # the final interval, then clear its denominator before serialization.
        zero_g = [Fraction(n*(n-2))] + [Fraction(1)]*(n-3)
        weighted = sum(k*x for k,x in enumerate(zero_g,2))
        cumulative = sum(k*(n-k)*x for k,x in enumerate(zero_g,2))/Fraction(n-2)
        final = (2*cumulative-weighted)/n
        assert final > 0
        patterns['zero-gamma'] = [x*final.denominator for x in zero_g] + [Fraction(final.numerator)]
        for zero in range(n-1):
            values = [Fraction(1)]*(n-1); values[zero] = Fraction(0)
            patterns['zero'+str(zero)] = values
        for name,g in patterns.items():
            for mode in ['comb','alternating']:
                root=ranked_tree(g,mode)
                for variant,labels,stem,unary,scale in [
                    ('labeled','all',Fraction(0),False,Fraction(1)),
                    ('unlabeled','tips',Fraction(0),False,Fraction(1)),
                    ('none','none',Fraction(0),False,Fraction(1)),
                    ('mixed-scaled','mixed',Fraction(17,8),False,Fraction(1,8)),
                    ('duplicate-large-stem','duplicate',Fraction(10**100),False,Fraction(8)),
                    ('unary-prefix','tips',Fraction(123),True,Fraction(1))]:
                    nw=serialize(root,labels,stem,unary,scale)
                    expected=oracle_newick(nw)
                    schedule=gamma_from_intervals(g)
                    assert abs(expected['value']-schedule['value']) < 1e-14
                    data.append({'id':f'{n}-{name}-{mode}-{variant}','newick':nw,'oracle':expected})
        # True multifurcation, equivalent to zero-length sequential binary resolutions.
        star=Spec(Fraction(0),[Spec(Fraction(3)) for _ in range(n)])
        for labels in ['all','tips','none']:
            nw=serialize(star,labels)
            data.append({'id':f'{n}-star-{labels}','newick':nw,'oracle':oracle_newick(nw)})
        # A mixed polytomy: root has n-1 children, one of which splits later.
        mixed=Spec(Fraction(0),[Spec(Fraction(3)) for _ in range(n-2)]+[Spec(Fraction(2),[Spec(Fraction(3)),Spec(Fraction(3))])])
        for labels in ['all','tips']:
            nw=serialize(mixed,labels,stem=Fraction(17))
            data.append({'id':f'{n}-mixed-polytomy-{labels}','newick':nw,'oracle':oracle_newick(nw)})
    return data

if __name__ == '__main__':
    data=cases()
    Path(__file__).with_name('exact-cases.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'case_count':len(data),'tip_counts':sorted({d['oracle']['tip_count'] for d in data})}))
