"""Independent inertia-orbit and stabilizer-quotient Euler-factor verification."""
from pathlib import Path
from math import gcd,prod
from collections import Counter
import json
from sympy import factorint
from eidetic import cycle_lengths
ROOT=Path(__file__).resolve().parent

def units(m):return [u for u in range(1,m+1) if gcd(u,m)==1]
def orbit_partition(points,perms):
    remaining=set(points);out=[]
    while remaining:
        todo=[min(remaining)];orb=set(todo)
        for x in todo:
            for pm in perms:
                y=pm[x]
                if y not in orb:orb.add(y);todo.append(y)
        assert orb<=set(points)
        remaining-=orb;out.append(sorted(orb))
    return out

def group_factors(row):
    maps=[[x-1 for x in pm] for pm in row['local_power_maps']]
    def act(u,x):return maps[(u-1)%row['class_orders'][x]][x]
    records=[];rational=[]
    for m in sorted(set(row['class_orders'])):
        layer=[i for i,o in enumerate(row['class_orders']) if o==m]
        us=units(m);perms=[{i:act(u,i) for i in layer} for u in us]
        for orb in orbit_partition(layer,perms):
            stabilizer=[u for u in us if act(u,orb[0])==orb[0]]
            rational.append((m,orb,stabilizer))
    for p in sorted(map(int,factorint(row['exponent']))):
        factors=[];second=[];details=[];completion=[]
        for m,orb,stabilizer in rational:
            mp=m
            while mp%p==0 and mp>1:mp//=p
            pa=m//mp;us=units(m)
            inertia=[u for u in us if (u-1)%mp==0]
            kappa=next(u for u in us if (u-p)%mp==0 and (u-1)%pa==0)
            ps=[{i:act(u,i) for i in orb} for u in inertia]
            iorbs=orbit_partition(orb,ps);idx={x:i for i,o in enumerate(iorbs) for x in o}
            permutation=[idx[act(kappa,o[0])] for o in iorbs]
            assert sorted(permutation)==list(range(len(permutation)))
            first=cycle_lengths(permutation)
            si={(s*u)%m for s in stabilizer for u in inertia}
            v=kappa%m;f=1
            while v not in si:v=v*kappa%m;f+=1;assert f<=len(us)
            g=len(us)//len(si)//f
            quotient=[f]*g
            assert first==quotient,(row['label'],p,m,first,quotient)
            factors+=first;second+=quotient
            if m%p==0:completion+=first
            details.append(dict(element_order=m,classes=[row['class_names'][i] for i in orb],
                                stabilizer=stabilizer,inertia=inertia,kappa=kappa,
                                inertia_orbits=iorbs,frobenius_cycles=first,residue_degree=f,prime_count=g))
        # A separately computed determinant of P_p uses only p-regular cycles.
        pm=maps[p-1];uncompleted=cycle_lengths(pm)
        assert sorted(uncompleted+completion)==sorted(factors)
        records.append(dict(prime=p,cycle_factors=dict(sorted(Counter(factors).items())),
                            completed_degree=sum(factors),uncompleted_degree=sum(uncompleted),
                            correction_degree=sum(completion),routes_agree=True,rational_classes=details))
    return dict(label=row['label'],rational_class_count=len(rational),bad_primes=records)

def cyclic_row(moduli):
    from itertools import product
    points=list(product(*(range(m) for m in moduli)));index={x:i for i,x in enumerate(points)}
    from math import lcm
    orders=[lcm(*(m//gcd(m,a) for m,a in zip(moduli,x))) for x in points]
    return dict(label='x'.join('C'+str(m) for m in moduli),class_count=len(points),class_sizes=[1]*len(points),
                class_orders=orders,exponent=lcm(*moduli),class_names=[str(x) for x in points],
                local_power_maps=[[index[tuple(n*a%m for m,a in zip(moduli,x))]+1 for x in points] for n in range(1,max(orders)+1)])

def main():
    rows=json.loads((ROOT/'data/sporadic_tables.json').read_text())['groups'];out=[]
    for r in rows:
        x=group_factors(r);out.append(x);print(r['label'],len(x['bad_primes']),flush=True)
    controls=[]
    for n in range(1,73):
        row=cyclic_row([n]);r=group_factors(row)
        for slot in r['bad_primes']:
            p=slot['prime']
            # For Q[C_n], each order-d rational class is Q(zeta_d).
            expected=[]
            for d in range(1,n+1):
                if n%d:continue
                m=d
                while m%p==0 and m>1:m//=p
                if m==1:f=1
                else:
                    f=1;v=p%m
                    while v!=1:v=v*p%m;f+=1
                expected += [f]*(len(units(m))//f)
            assert dict(Counter(expected))==slot['cycle_factors']
        controls.append(dict(label=row['label'],verified=True))
    for moduli in ([2,2],[2,4],[3,3],[2,6],[4,4],[3,9],[2,2,2],[2,3,5]):
        row=cyclic_row(moduli);group_factors(row);controls.append(dict(label=row['label'],verified=True))
    result=dict(schema='eidetic-ramified-factors-v1',groups=out,slot_count=sum(len(x['bad_primes']) for x in out),all_routes_agree=True,controls=controls,
                note='All local unit maps re-exported; no slots omitted because of historical sparse-map availability.')
    (ROOT/'data/ramified_local_factors.json').write_text(json.dumps(result,indent=2)+'\n');print('Verified',result['slot_count'],'slots and',len(controls),'controls')
if __name__=='__main__':main()
