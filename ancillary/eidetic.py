"""Reconstructed certificate layer for the current Eidetic paper.

Conventions: stored GAP maps are one-based. Matrices act on columns, with
P_n[i, power_n(i)] = 1 and P_n^* = D^-1 P_n^T D.
No historical JSON result is used as evidence by this module.
"""
from collections import Counter
from fractions import Fraction
from math import gcd
import numpy as np
from flint import nmod_mat
from sympy import factorint
from proof_construction_assembled.experiments.finite_group_adams_star_census_20260724.run_census import ModularBasis
from proof_construction_assembled.experiments.sporadic_adams_star_census_20260724.run_census import stable_diagonal_partition, partition_cells

PRIMES = (1000003, 1000033)

def validate_input(row):
    k = row['class_count']; weights = row['class_sizes']; maps = row['power_maps']
    assert k == len(weights) == len(row['class_orders'])
    assert all(isinstance(w,int) and w > 0 for w in weights)
    assert sum(weights) == row['order']
    assert len(maps) == len(row['directions'])
    expected = {int(p)**a for p,e in factorint(row['exponent']).items() for a in range(1,int(e)+1)}
    assert set(row['directions']) == expected
    by_direction = dict(zip(row['directions'], maps))
    for p,e in factorint(row['exponent']).items():
        base = by_direction[int(p)]
        current = list(range(1,k+1))
        for a in range(1,int(e)+1):
            current = [base[t-1] for t in current]
            assert current == by_direction[int(p)**a]
    identity = [i for i,o in enumerate(row['class_orders']) if o == 1]
    assert len(identity) == 1 and weights[identity[0]] == 1
    for n,mapping in zip(row['directions'],maps):
        assert len(mapping) == k and all(1 <= t <= k for t in mapping)
        assert all(row['class_orders'][mapping[i]-1] == o//gcd(o,n) for i,o in enumerate(row['class_orders']))
    for a in maps:
        for b in maps:
            assert all(a[b[i]-1] == b[a[i]-1] for i in range(k))
    return identity[0]

def generators(row, prime):
    """Only prime divisors generate P_ram; higher prime powers are redundant."""
    weights = row['class_sizes']; k = row['class_count']
    if any(w % prime == 0 for w in weights):
        raise ValueError('bad reduction prime divides a weight')
    primes = sorted(map(int,factorint(row['exponent'])))
    mapping = dict(zip(row['directions'],row['power_maps']))
    mats=[]; labels=[]
    for p in primes:
        a=np.zeros((k,k),dtype=np.int64)
        for i,t in enumerate(mapping[p]): a[i,t-1]=1
        mats.append(a);labels.append(f'P{p}')
    for p,a in zip(primes,mats.copy()):
        b=np.zeros((k,k),dtype=np.int64)
        for i in range(k):
            for j in range(k):
                if a[j,i]:b[i,j]=(weights[j]%prime)*pow(weights[i]%prime,-1,prime)%prime
        mats.append(b);labels.append(f'P{p}*')
    return labels,mats

def cyclic_words(mats, seed, prime, coordinates=None):
    """Store actual generator words and a square nonzero minor, not just rank."""
    coordinates=list(range(len(seed))) if coordinates is None else coordinates
    basis=ModularBasis(len(coordinates),prime)
    selected=[];vectors=[];queue=[]
    def add(v,w):
        if basis.add(v[coordinates]):
            vectors.append(v.copy());selected.append(w);queue.append((v,w))
    add(np.asarray(seed,dtype=np.int64)%prime,[])
    at=0
    while at<len(queue):
        v,w=queue[at];at+=1
        for j,a in enumerate(mats):add(a@v%prime,w+[j])
    pivots=sorted(basis.rows)
    square=[[int(v[coordinates[i]]) for v in vectors] for i in pivots]
    determinant=int(nmod_mat(square,prime).det()) if square else 1
    assert determinant != 0
    return dict(prime=prime,rank=len(vectors),words=selected,
                coordinate_rows=[coordinates[i] for i in pivots],minor_determinant=determinant,
                maximum_word_length=max(map(len,selected),default=0))

def replay_cyclic(row, certificate, seed):
    prime=certificate['prime'];_,mats=generators(row,prime)
    vs=[]
    for word in certificate['words']:
        v=np.asarray(seed,dtype=np.int64)%prime
        for j in word:v=mats[j]@v%prime
        vs.append(v)
    square=[[int(v[i]) for v in vs] for i in certificate['coordinate_rows']]
    assert len(square)==len(vs)==certificate['rank']
    det=int(nmod_mat(square,prime).det())
    assert det==certificate['minor_determinant'] and det!=0
    return True

def principal_certificate(row, allow_proper=False):
    identity=validate_input(row)
    labels,history=stable_diagonal_partition(row['class_sizes'],row['power_maps'])
    cells=partition_cells(labels); reps=[c[0] for c in cells]
    assert cells[labels[identity]]==[identity]
    seed=[int(i==identity) for i in range(row['class_count'])]
    receipts=[]
    for prime in PRIMES:
        names,mats=generators(row,prime)
        # Exact stable partition gives invariant cell-constant space; verify reductions too.
        for a in mats:
            for cell in cells:
                totals=[tuple(int(sum(a[i,j] for j in target))%prime for target in cells) for i in cell]
                assert len(set(totals))==1
        cert=cyclic_words(mats,seed,prime,reps)
        if not allow_proper:
            assert cert['rank']==len(cells), (row['label'],cert['rank'],len(cells))
        replay_cyclic(row,cert,seed);receipts.append(cert)
    return dict(label=row['label'],principal_dimension=receipts[0]['rank'],stable_partition_dimension=len(cells),cells=cells,
                refinement_history=history,generator_labels=names,identity_index=identity,
                projection_formula='P_exp^* P_exp / group_order',certificates=receipts)

def root_profiles(row):
    weights=row['class_sizes'];k=row['class_count'];profiles=[[] for _ in weights]
    for mapping in row['power_maps']:
        totals=[0]*k
        for i,t in enumerate(mapping):totals[t-1]+=weights[i]
        for i in range(k):profiles[i].append(Fraction(totals[i],weights[i]))
    return [tuple(p) for p in profiles]

def cycle_lengths(mapping):
    """Cycles of a finite function; transient vertices contribute no det factor."""
    done=set(); cycles=[]
    for start in range(len(mapping)):
        if start in done:continue
        positions={};v=start
        while v not in done and v not in positions:
            positions[v]=len(positions);v=mapping[v]
        if v in positions:cycles.append(len(positions)-positions[v])
        done.update(positions)
    return sorted(cycles)
