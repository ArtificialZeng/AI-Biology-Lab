"""Fixed-rate mixed phosphorylation model; source arXiv:1809.02886v2.

Six coordinates are (S0, S0K, S1K, S2, S2P, S1P).
The eliminated species follow the three exact conservation laws.
"""
import numpy as np

RATES = np.array([1, 1, 1, 1, 100, 1, .9, 3, 1, 100], dtype=float)

def finite_rhs(t, y, K, P=8., S=40.):
    s0, u, v, w, a, b = y
    kinase = K-u-v
    phosphatase = P-a-b
    s1 = S-s0-u-v-w-a-b
    q = s0*kinase
    return np.array([-q+u+100*b, q-2*u, u-v,
                     v-100*w*phosphatase+a,
                     100*w*phosphatase-1.9*a,
                     3*s1*phosphatase-101*b])

def finite_jac(t, y, K, P=8., S=40.):
    s0,u,v,w,a,b=y
    z=K-u-v; p=P-a-b; r=S-s0-u-v-w-a-b
    return np.array([
        [-z,s0+1,s0,0,0,100],
        [z,-s0-2,-s0,0,0,0],
        [0,1,-1,0,0,0],
        [0,0,1,-100*p,100*w+1,100*w],
        [0,0,0,100*p,-100*w-1.9,-100*w],
        [-3*p,-3*p,-3*p,-3*p,-3*(r+p),-3*(r+p)-101]])

def limit_rhs(t, y, P=8., S=40.):
    u,v,w,a,b=y
    p=P-a-b; r=S-u-v-w-a-b
    return np.array([100*b-u,u-v,v-100*w*p+a,100*w*p-1.9*a,3*r*p-101*b])

def limit_jac(t, y, P=8., S=40.):
    u,v,w,a,b=y
    p=P-a-b; r=S-u-v-w-a-b
    return np.array([[-1,0,0,0,100],[1,-1,0,0,0],
                     [0,1,-100*p,100*w+1,100*w],
                     [0,0,100*p,-100*w-1.9,-100*w],
                     [-3*p,-3*p,-3*p,-3*(r+p),-3*(r+p)-101]])

def physical_species(y,K,P=8.,S=40.):
    s0,u,v,w,a,b=y
    return np.array([s0,K-u-v,u,v,w,P-a-b,a,S-s0-u-v-w-a-b,b])
