"""Independent group-law construction for the dihedral/quaternion ladder."""
from math import gcd,lcm
from pathlib import Path
import json
from sympy import factorint
from runners.experiments.finite_group_adams_star_census_20260724.run_census import analyze_group
ROOT=Path(__file__).resolve().parent

def metacyclic_row(n,twist=0):
    # a^n=1, b a b^-1=a^-1, b^2=a^twist (twist=0 or n/2).
    elements=[(i,j) for j in range(2) for i in range(n)]
    def mul(x,y):return ((x[0]+(-1 if x[1] else 1)*y[0]+(twist if x[1] and y[1] else 0))%n,(x[1]+y[1])%2)
    def power(x,k):
        out=(0,0)
        while k:
            if k&1:out=mul(out,x)
            x=mul(x,x);k//=2
        return out
    inverses={x:next(y for y in elements if mul(x,y)==(0,0)) for x in elements}
    cells=[];seen=set()
    for x in elements:
        if x in seen:continue
        cell=sorted({mul(mul(h,x),inverses[h]) for h in elements});seen.update(cell);cells.append(cell)
    index={x:i for i,c in enumerate(cells) for x in c}
    orders=[next(k for k in range(1,2*n+1) if power(c[0],k)==(0,0)) for c in cells]
    exponent=lcm(*orders);directions=[p**a for p,e in factorint(exponent).items() for a in range(1,e+1)]
    return dict(label=('Q' if twist else 'D')+str(2*n),family='quaternion' if twist else 'dihedral',
                order=2*n,class_count=len(cells),exponent=exponent,directions=directions,
                class_sizes=list(map(len,cells)),class_orders=orders,
                power_maps=[[index[power(c[0],p)]+1 for c in cells] for p in directions])

def cyclic_row(p,a):
    n=p**a
    return dict(label=f'C{n}',family='cyclic_prime_power',order=n,class_count=n,exponent=n,
                directions=[p**i for i in range(1,a+1)],class_sizes=[1]*n,
                class_orders=[n//gcd(n,i) for i in range(n)],
                power_maps=[[(p**j*i)%n+1 for i in range(n)] for j in range(1,a+1)])

def main():
    out=[]
    for row in [metacyclic_row(n) for n in range(3,17)]+[metacyclic_row(2*m,m) for m in range(2,9)]:
        r=analyze_group(row);out.append(r);print(r['label'],r['algebra_dimension'],r['center_dimension'],r['commutant_dimension'],flush=True)
    controls=[]
    for p,a in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(5,2)]:
        r=analyze_group(cyclic_row(p,a))
        expected=(a+1)**2+sum(j*j for j in range(1,a+(p!=2)))
        assert r['algebra_dimension']==expected
        controls.append(r)
    (ROOT/'data/families.json').write_text(json.dumps(dict(schema='eidetic-families-v1',ladder=out,cyclic_controls=controls),indent=2)+'\n')
if __name__=='__main__':main()
