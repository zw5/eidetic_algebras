from pathlib import Path
from copy import deepcopy
from collections import Counter
import json
import numpy as np
import pytest
from flint import fmpq,fmpq_mat
from eidetic import validate_input, replay_cyclic, generators, PRIMES, cycle_lengths
from monster_certificate import certify
from local_factors import group_factors,cyclic_row
from characters import group_characters
from runners.experiments.finite_group_adams_star_census_20260724.run_census import analyze_group,exact_pullback_matrix,exact_adjoint_matrix
from runners.experiments.sporadic_adams_star_census_20260724.run_census import analyze_group as sporadic_analyze,exact_pullback,exact_adjoint
ROOT=Path(__file__).resolve().parents[1]
def load(name):return json.loads((ROOT/'data'/name).read_text())
TABLES={r['label']:r for r in load('sporadic_tables.json')['groups']}
CENSUS={r['label']:r for r in load('census.json')['groups']}

@pytest.mark.parametrize('label',list(TABLES))
def test_table_integrity(label):validate_input(TABLES[label])

@pytest.mark.parametrize('label',list(TABLES))
def test_principal_certificate_replay(label):
    c=next(c for c in load('principal_certificates.json') if c['label']==label);r=TABLES[label]
    seed=[int(o==1) for o in r['class_orders']]
    for receipt in c['certificates']:assert replay_cyclic(r,receipt,seed)
    assert c['principal_dimension']==CENSUS[label]['principal_block_size']

def test_corrupted_minor_is_rejected():
    r=TABLES['M11'];c=deepcopy(load('principal_certificates.json')[0]['certificates'][0]);c['minor_determinant']+=1
    with pytest.raises(AssertionError):replay_cyclic(r,c,[int(o==1) for o in r['class_orders']])

def test_corrupted_map_is_rejected():
    r=deepcopy(TABLES['M11']);r['power_maps'][0][1]=2
    with pytest.raises(AssertionError):validate_input(r)

def test_bad_reduction_is_rejected():
    with pytest.raises(ValueError):generators(TABLES['M11'],2)

def test_monster_fresh_full_certificate():
    result=certify(TABLES['Monster'])
    assert result['certificate_vector_count']==186
    assert result['symmetry_rank']==14 and result['table_automorphism_rank']==13
    assert result['algebra_dimension']==28958
    assert sorted(result['block_degrees'])==[1]*8+[2]*4+[3,5,170]
    assert {tuple(p) for p in result['exceptional_component']}=={('16b','16c'),('32a','32b')}

def test_fi23_exact_counterexample_and_corrected_decomposition():
    row=TABLES['Fi23'];c=load('fi23_independent_check.json');v=fmpq_mat([[fmpq(x)] for x in c['invariant_vector']])
    assert c['nonzero_entries']=={'12h':'1','12l':'-1/3','12m':'-1/4','12o':'1/12'}
    assert v[0,0]==0
    for pm in row['power_maps']:
        a=exact_pullback(pm)
        for g in (a,exact_adjoint(a,row['class_sizes'])):
            x=g*v;idx=next(i for i in range(v.nrows()) if v[i,0]);lam=x[idx,0]/v[idx,0]
            assert x==v*lam
    result=sporadic_analyze(row,extra_vector=c['invariant_vector'])
    assert (result['algebra_dimension'],result['center_dimension'],result['commutant_dimension'])==(8294,6,6)
    assert sorted(f['matrix_size_over_center'] for f in result['rational_wedderburn_factors'])==[1,1,1,1,3,91]
    with pytest.raises(ValueError):sporadic_analyze(row)

def test_d8_q8_and_exact_weighted_transfer():
    path=ROOT/'runners/experiments/finite_group_adams_star_census_20260724/output/gap_groups.json'
    rows={r['label']:r for r in json.loads(path.read_text())['groups']}
    for name,expected in [('D8',17),('Q8',10)]:
        row=rows[name];r=analyze_group(row);assert r['algebra_dimension']==expected
        for pm in row['power_maps']:
            a=exact_pullback_matrix(pm);star=exact_adjoint_matrix(a,row['class_sizes']);gram=star*a
            assert all(not gram[i,j] for i in range(a.nrows()) for j in range(a.nrows()) if i!=j)
            expected_roots=[sum(row['class_sizes'][j] for j,t in enumerate(pm) if t==i+1) for i in range(len(pm))]
            assert all(gram[i,i]==fmpq(n,row['class_sizes'][i]) for i,n in enumerate(expected_roots))

def test_ramified_two_routes_all_sporadics():
    regenerated=[group_factors(r) for r in TABLES.values()]
    assert sum(len(r['bad_primes']) for r in regenerated)==173
    # Prime cyclotomic local factor at its own prime is (1-T)^2 for Q[C_p].
    for p in (2,3,5,7,11):assert group_factors(cyclic_row([p]))['bad_primes'][0]['cycle_factors']=={1:2}

def test_cycle_enumeration_ignores_transient_trees():
    assert cycle_lengths([1,2,1,3,3])==[1,2]

def test_dynamics_separates_all_26():
    packets=[tuple(sorted(tuple(sorted(Counter(cycle_lengths([x-1 for x in pm])).items())) for pm in r['power_maps'])) for r in TABLES.values()]
    assert len(set(packets))==26
    assert sum(len(r['power_maps']) for r in TABLES.values())==258

@pytest.mark.parametrize('name',['J1','Ly','Monster','Fi23'])
def test_fresh_character_reconstruction(name):
    r=group_characters(TABLES[name],CENSUS[name]);assert sum(b['block_dimension_over_Q'] for b in r['blocks'])==TABLES[name]['class_count']
    if name=='Ly':assert {(b['character_order'],b['character_conductor']) for b in r['blocks'] if b['center_degree']>1}=={(3,67),(5,31)}
    if name=='Monster':assert r['rational_class_count']==172

def test_all_wedderburn_identities_and_centers():
    fields=Counter()
    for row in CENSUS.values():
        fs=row['rational_wedderburn_factors'];fields.update(f['center_field'] for f in fs)
        assert all(f['schur_index']==1 for f in fs)
        assert sum(f['center_degree']*f['matrix_size_over_center']**2 for f in fs)==row['algebra_dimension']
        assert sum(f['center_degree']*f['module_multiplicity']**2 for f in fs)==row['commutant_dimension']
        assert sum(f['center_degree']*f['module_multiplicity']*f['matrix_size_over_center'] for f in fs)==row['class_count']
    assert fields['Q(sqrt(-3))']==9 and fields['Q(zeta_5)']==1
