"""Fresh replayable Monster saturation certificate; no old result inputs."""
from pathlib import Path
from types import SimpleNamespace
from collections import Counter
import json
import numpy as np
from eidetic import *
from symmetry_helpers import _refined_cells,_swap_system,_rank_mod
ROOT=Path(__file__).resolve().parent

def certify(row):
    k=row['class_count'];names=row['class_names'];weights=row['class_sizes']
    primes=sorted(map(int,factorint(row['exponent'])))
    maps=[[x-1 for x in row['power_maps'][row['directions'].index(p)]] for p in primes]
    registry=SimpleNamespace(names=names,orders=row['class_orders'],class_sizes=weights)
    _,cells,rounds=_refined_cells(registry,primes,maps)
    swaps=_swap_system(cells,maps)
    assert not swaps.pinned_components
    components=[[list(swaps.pair_cells[j]) for j in comp] for comp in swaps.components]
    table_bits=[]
    for perm in row['table_automorphisms']:
        assert all(weights[i]==weights[perm[i]-1] for i in range(k))
        table_bits.append([int(perm[c[0][0]]-1==c[0][1]) for c in components])
    table_rank=_rank_mod(table_bits,2)
    good_bits=[[int(pm[c[0][0]]-1==c[0][1]) for c in components] for pm in row['unit_power_maps']]
    good_rank=_rank_mod(good_bits,2)
    principal=principal_certificate(row)
    profiles=root_profiles(row)
    receipts=[];scalar=[]
    for component in components:
        if len(component)==1:
            a,b=component[0];seed=np.zeros(k,dtype=np.int64);seed[a]=1;seed[b]=-1
            _,gs=generators(row,PRIMES[0]);values=[]
            for g in gs:
                v=g@seed%PRIMES[0];lam=int(v[a]);assert np.all(v==seed*lam%PRIMES[0]);values.append(lam)
            scalar.append(dict(pair=[names[a],names[b]],eigenvalues=values));continue
        isolated=[]
        for a,b in component:
            support=[i for i,p in enumerate(profiles) if p==profiles[a]]
            if support==[a,b]:isolated.append((a,b))
        assert isolated,component
        a,b=isolated[0];seed=np.zeros(k,dtype=np.int64);seed[a]=1;seed[b]=-1
        cs=[]
        for prime in PRIMES:
            _,gs=generators(row,prime)
            c=cyclic_words(gs,seed,prime,[a for a,b in component])
            assert c['rank']==len(component)
            replay_cyclic(row,c,seed);cs.append(c)
        receipts.append(dict(seed_pair=[names[a],names[b]],seed_indices=[a,b],dimension=len(component),
                            root_profile=[str(x) for x in profiles[a]],certificates=cs))
    assert len({tuple(s['eigenvalues']) for s in scalar})==len(scalar)
    # Verify the paper's sharper scalar separation by P2,P5,P7.
    inds=[principal['generator_labels'].index(f'P{p}') for p in (2,5,7)]
    assert len({tuple(s['eigenvalues'][i] for i in inds) for s in scalar})==len(scalar)
    exceptional=next(c for c in components if any(names[a].lower()=='16b' for a,b in c))
    bit=[int(c==exceptional) for c in components]
    assert _rank_mod(table_bits+[bit],2)>table_rank
    degrees=[len(cells)]+[len(c) for c in components]
    return dict(status='verified',group=row['label'],symmetry_rank=len(components),
                stable_cells=len(cells),refinement_rounds_including_fixed_round=rounds,
                table_automorphism_rank=table_rank,good_label_image_rank=good_rank,
                swap_components=[[[names[a],names[b]] for a,b in c] for c in components],
                exceptional_component=[[names[a],names[b]] for a,b in exceptional],
                block_degrees=degrees,algebra_dimension=sum(d*d for d in degrees),
                certificate_vector_count=principal['principal_dimension']+sum(x['dimension'] for x in receipts),
                principal=principal,nontrivial_cyclic_certificates=receipts,scalar_channels=scalar)

def main():
    rows=json.loads((ROOT/'data/sporadic_tables.json').read_text())['groups']
    result=certify(next(r for r in rows if r['label']=='Monster'))
    (ROOT/'data/monster_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:result[k] for k in ['status','symmetry_rank','table_automorphism_rank','good_label_image_rank','block_degrees','algebra_dimension','certificate_vector_count']})
if __name__=='__main__':main()
