"""Full16parameter reducedmatrix transcribed from source §2, no symmetry imposed."""
NAMES=('a','b4','b10','c','e1','e4','f7','f10','u2','u3','v5','v6','x8','x9','y11','y12')
def matrix(p,constructor=list):
    a,b4,b10,c,e1,e4,f7,f10,u2,u3,v5,v6,x8,x9,y11,y12=p
    return constructor([[0,0,-u3,0,y12,0],[0,0,0,v6,0,-x9],
       [a,0,-a-e1-u2-u3,-e1,0,0],[-b4,-b4,-e4,-b4-e4-v5-v6,-b4,0],
       [-b10,-b10,0,-b10,-b10-f10-y11-y12,-f10],[0,c,0,0,-f7,-c-f7-x8-x9]])
