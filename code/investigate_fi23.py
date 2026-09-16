"""Exact independent check of the Fi23 principal-block assumption."""
import json
from pathlib import Path
from flint import fmpq, fmpq_mat
from eidetic import *
from runners.experiments.sporadic_adams_star_census_20260724.run_census import exact_pullback,exact_adjoint
root=Path(__file__).resolve().parent
r=next(x for x in json.loads((root/'data/sporadic_tables.json').read_text())['groups'] if x['label']=='Fi23')
labels,h=stable_diagonal_partition(r['class_sizes'],r['power_maps']);cells=partition_cells(labels);reps=[x[0] for x in cells]
ps=sorted(map(int,factorint(r['exponent'])));mapping=dict(zip(r['directions'],r['power_maps']))
mats=[exact_pullback(mapping[p]) for p in ps];mats += [exact_adjoint(a,r['class_sizes']) for a in mats.copy()]
_, modular_generators = generators(r, PRIMES[0])
c = cyclic_words(modular_generators, [int(o==1) for o in r['class_orders']], PRIMES[0], reps)
(root/'data/fi23_orbit_mod.json').write_text(json.dumps(c,indent=2)+'\n')
vs=[]
for word in c['words']:
 v=fmpq_mat([[int(o==1)] for o in r['class_orders']])
 for j in word:v=mats[j]*v
 vs.append([v[i,0] for i in reps])
a=fmpq_mat(vs);rr,rank=a.rref();print('Exact rank',rank,flush=True)
piv=[]
for i in range(rank):piv.append(next(j for j in range(len(reps)) if rr[i,j]))
free=[i for i in range(len(reps)) if i not in piv];assert len(free)==1
v=[fmpq(0)]*len(reps);v[free[0]]=fmpq(1)
for i,j in enumerate(piv):v[j]=-rr[i,free[0]]
# annihilator covector becomes vector by inverse quotient weights
w=[sum(r['class_sizes'][j] for j in cell) for cell in cells]
u=[v[i]/w[i] for i in range(len(cells))];scale=next(x for x in u if x);u=[x/scale for x in u]
full=fmpq_mat([[u[labels[i]]] for i in range(r['class_count'])]); eigen=[]
idx=next(i for i in range(r['class_count']) if full[i,0])
for a in mats:
 image=a*full;lam=image[idx,0]/full[idx,0]
 assert image==full*lam
 eigen.append(str(lam))
result={'exact_orbit_rank':rank,'proposed_principal_dimension':len(cells),'invariant_vector':[str(full[i,0]) for i in range(r['class_count'])], 'nonzero_entries':{r['class_names'][i]:str(full[i,0]) for i in range(r['class_count']) if full[i,0]},'generator_labels':[f'P{p}' for p in ps]+[f'P{p}*' for p in ps],'eigenvalues':eigen}
(root/'data/fi23_independent_check.json').write_text(json.dumps(result,indent=2)+'\n');print(result['nonzero_entries'],eigen)
