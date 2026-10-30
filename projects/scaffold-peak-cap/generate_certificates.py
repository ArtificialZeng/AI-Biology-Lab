"""Rational construction and exact finite radical certificates; no floats."""
from fractions import Fraction as Q
from math import isqrt
import argparse
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def centers(p):
    require(type(p) is int and p >= 2, 'integer p >= 2 required')
    return [Q(1, 10)] + [Q(2*j-3, 2*j-1) for j in range(2, p+1)]


def witnesses(p):
    return [(Q(4*j-7, 4*j-3), 1, Q(1, 4))
            if side == 0 else (Q(4*j-5, 4*j-1), -1, Q(-1, 4))
            for j in range(2, p+1) for side in (0, 1)] + [
                (Q(p, p+1), 1, Q(1)), (Q(3), -1, Q(3, 2)-p)]


def error(c, x):
    return sum((x*t/(x-t)**2 for t in c), Q(0)) + x*(1+x)/abs(x-1)**3


def ceil_q(x):
    return -((-x.numerator)//x.denominator)


def bracket_sqrt(z, bits):
    require(z > 0 and type(bits) is int and bits >= 1, 'sqrt inputs')
    scale = 1 << bits
    n = isqrt((z.numerator*scale*scale)//z.denominator)
    lo, hi = Q(n, scale), Q(n+1, scale)
    if lo*lo == z:
        lo -= Q(1, scale)
    require(0 <= lo and lo*lo < z < hi*hi, 'sqrt isolation')
    return lo, hi


def polynomial(u, q, alpha):
    return u*u + (q-alpha-1)*u-q


def make(p, bits=160):
    c, w = centers(p), witnesses(p)
    xs = [x for x, _, _ in w]
    require(all(a < b for a,b in zip(xs,xs[1:])), 'witness order')
    E = [error(c, x) for x in xs]
    emax = max(E)
    M = 1 + ceil_q(8*emax)
    K = [M/t for t in c] + [Q(M)]
    alpha = [Q(M-1)]*(p+1)
    evidence = []
    for (x, sign, target), e in zip(w, E):
        roots, terms = [], []
        for k, a in zip(K, alpha):
            q = k*x
            rad = (q-a-1)**2+4*q
            sl, sh = bracket_sqrt(rad, bits)
            lo, hi = (a+1-q+sl)/2, (a+1-q+sh)/2
            require(1 < lo < hi and polynomial(lo,q,a) < 0 < polynomial(hi,q,a),
                    'positive quadratic branch')
            roots.append([str(lo), str(hi)])
            terms.append((q/(hi*hi+q), q/(lo*lo+q)))
        ulo, uhi = map(Q, roots[-1])
        q = K[-1]*x
        mlo = (ulo-1)*q/(uhi*uhi+q)
        mhi = (uhi-1)*q/(ulo*ulo+q)
        flo = 1-sum((t[1] for t in terms[:-1]), Q(0))+mlo
        fhi = 1-sum((t[0] for t in terms[:-1]), Q(0))+mhi
        require((sign == 1 and flo > 0) or (sign == -1 and fhi < 0), 'exact sign')
        require(e/M < Q(1,8), 'finite-M error budget')
        evidence.append({'X':str(x), 'expected_sign':sign, 'F_infty':str(target),
                         'E':str(e), 'error_bound':str(e/M),
                         'U_brackets':roots, 'F_bounds':[str(flo),str(fhi)]})
    return {'p':p,'m':p+1,'occupied':list(range(1,p+1)), 'vacant':[p+1],
            'centers':list(map(str,c)),'M':M,'Emax':str(emax),
            'K':list(map(str,K)),'alpha':list(map(str,alpha)),
            'sqrt_bits':bits,'method':'strict Fraction dyadic radical isolation of original U-quadratic',
            'witnesses':evidence}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--minimum-p',type=int,default=2)
    parser.add_argument('--maximum-p',type=int,default=6)
    parser.add_argument('--bits',type=int,default=160)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('cases.json'))
    args = parser.parse_args()
    require(2 <= args.minimum_p <= args.maximum_p, 'finite range')
    docs = [make(p,args.bits) for p in range(args.minimum_p,args.maximum_p+1)]
    doc = {'schema':'one-vacancy-rational-radical-v1','cases':docs,
           'scope':'Finite checks supplement the all-integer-p proof; no floating-point evidence.'}
    args.output.write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps({'written':str(args.output),'p':[d['p'] for d in docs],
                      'M':[d['M'] for d in docs]}))


if __name__ == '__main__':
    main()
