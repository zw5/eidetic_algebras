"""Exact block characters reconstructed from central idempotents and GAP unit maps."""
from pathlib import Path
from math import gcd,lcm
import json
from flint import fmpq,fmpq_mat
from sympy import factorint,kronecker_symbol
from symmetry_helpers import _local_unit_generators
from runners.experiments.sporadic_adams_star_census_20260724.run_census import residual_basis,stable_diagonal_partition,partition_cells,exact_pullback
ROOT=Path(__file__).resolve().parent

def matrix(rows):return fmpq_mat([[fmpq(x) for x in r] for r in rows])
def trace(a):return sum((a[i,i] for i in range(a.nrows())),fmpq(0))
def restriction_setup(row):
    labels,_=stable_diagonal_partition(row['class_sizes'],row['power_maps']);cells=partition_cells(labels)
    b,piv=residual_basis(row['class_sizes'],cells)
    if row['label']=='Fi23':
        v=json.loads((ROOT/'data/fi23_independent_check.json').read_text())['invariant_vector']
        b=fmpq_mat([[b[i,j] for j in range(b.ncols())]+[fmpq(v[i])] for i in range(b.nrows())])
    if b.ncols()==0:return b,[],None
    rr,rank=b.transpose().rref();piv=[next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(rank)]
    inv=fmpq_mat([[b[i,j] for j in range(rank)] for i in piv]).inv()
    return b,piv,inv

def restrict(pm,b,piv,inv):
    image=exact_pullback(pm)*b
    small=inv*fmpq_mat([[image[i,j] for j in range(image.ncols())] for i in piv])
    assert b*small==image
    return small

def center_coordinates(a,e,z,degree):
    bs=[];v=e
    for _ in range(degree):bs.append(v);v=v*z
    constraints=fmpq_mat([list(x.entries()) for x in bs]).transpose()
    _,rank=constraints.transpose().rref();assert rank==degree
    rr,_=constraints.transpose().rref();piv=[next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(degree)]
    square=fmpq_mat([[constraints[i,j] for j in range(degree)] for i in piv])
    values=list(a.entries());coeff=square.solve(fmpq_mat([[values[i]] for i in piv]))
    assert constraints*coeff==fmpq_mat([[x] for x in values]),'unit action is not a center scalar'
    return [str(coeff[i,0]) for i in range(degree)]

def action_order(a,e):
    x=e
    for n in range(1,121):
        x=x*a
        if x==e:return n
    raise ArithmeticError('unit image order exceeded bound')

def conductor(exponent,orders):
    at=0;answer=1
    for p,a in factorint(exponent).items():
        p=int(p);a=int(a);count=len(_local_unit_generators(p,a));os=orders[at:at+count];at+=count
        if all(o==1 for o in os):continue
        if p!=2:
            order=os[0];b=next(b for b in range(1,a+1) if ((p-1)*p**(b-1))%order==0)
        else:
            minus=os[0];five=os[1] if len(os)>1 else 1
            b=next(b for b in range(2,a+1) if (b>=2 or minus==1) and (2**max(b-2,0))%five==0)
        answer*=p**b
    assert at==len(orders)
    return answer

def group_characters(row,census):
    b,piv,inv=restriction_setup(row);unit=[restrict(pm,b,piv,inv) for pm in row['unit_power_maps']] if b.ncols() else []
    bad=[restrict(pm,b,piv,inv) for pm in row['power_maps']] if b.ncols() else []
    for u in unit:
        for a in bad:assert a*u==u*a
    blocks=[dict(center_field='Q',center_degree=1,block_dimension_over_Q=census['principal_block_size'],matrix_size_over_center=census['principal_block_size'],module_multiplicity=1,character_order=1,character_conductor=1,fundamental_discriminant=1,unit_values=['1']*len(row['unit_generators']),source='identity-class cyclic block')]
    cr=census['residual_algebra'].get('generic_center_certificate')
    if cr:
        z=matrix(cr['central_matrix'])
        for f,r in zip(census['residual_algebra']['rational_wedderburn_factors'],cr['factor_degrees_and_exponents']):
            e=matrix(r['idempotent']);assert e*e==e;dimension=int(trace(e))
            orders=[];coords=[];values=[]
            for u in unit:
                a=u*e;assert e*u==a
                coords.append(center_coordinates(a,e,z,f['center_degree']))
                orders.append(action_order(a,e))
                if f['center_degree']==1:
                    lam=trace(a)/dimension;assert a==e*lam;values.append(str(lam))
            cond=conductor(row['exponent'],orders);order=lcm(*orders) if orders else 1
            d=None
            if f['center_degree']==1:
                if order==1:d=1
                else:
                    candidates=[x for x in (cond,-cond) if all(int(kronecker_symbol(x,u))==int(fmpq(v)) for u,v in zip(row['unit_generators'],values))]
                    assert len(candidates)==1,(row['label'],cond,candidates)
                    d=candidates[0]
            block=dict(f,block_dimension_over_Q=dimension,character_order=order,character_conductor=cond,
                       fundamental_discriminant=d,unit_image_orders=orders,unit_center_coordinates=coords,
                       unit_values=values,acts_by_center_scalars=True)
            blocks.append(block)
    assert sum(b['block_dimension_over_Q'] for b in blocks)==row['class_count']
    # Rational classes are orbits of the full unit action, counted independently.
    seen=set();orbits=[]
    for i in range(row['class_count']):
        if i in seen:continue
        orbit={i};todo=[i]
        for j in todo:
            for pm in row['unit_power_maps']:
                x=pm[j]-1
                if x not in orbit:orbit.add(x);todo.append(x)
        seen.update(orbit);orbits.append(sorted(orbit))
    fixed=sum(b['block_dimension_over_Q'] for b in blocks if b['character_order']==1)
    assert fixed==len(orbits),(row['label'],fixed,len(orbits))
    return dict(label=row['label'],class_count=row['class_count'],unit_generators=row['unit_generators'],
                rational_class_count=len(orbits),rational_class_orbits=orbits,blocks=blocks)

def main():
    rows=json.loads((ROOT/'data/sporadic_tables.json').read_text())['groups']
    census={r['label']:r for r in json.loads((ROOT/'data/census.json').read_text())['groups']}
    out=[]
    for r in rows:
        x=group_characters(r,census[r['label']]);out.append(x)
        print(r['label'],[(b['block_dimension_over_Q'],b['fundamental_discriminant'],b['character_order'],b['character_conductor']) for b in x['blocks']],flush=True)
        (ROOT/'data/good_label_characters.json').write_text(json.dumps(dict(schema='eidetic-block-characters-v1',complete=len(out)==26,groups=out),indent=2)+'\n')
if __name__=='__main__':main()
