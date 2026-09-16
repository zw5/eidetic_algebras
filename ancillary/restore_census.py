"""Recompute the certified census, including the independently detected Fi23 correction."""
import json
from pathlib import Path
from eidetic import principal_certificate
from proof_construction_assembled.experiments.sporadic_adams_star_census_20260724.run_census import analyze_group
ROOT=Path(__file__).resolve().parent

def main():
    source=json.loads((ROOT/'data/sporadic_tables.json').read_text())
    correction=json.loads((ROOT/'data/fi23_independent_check.json').read_text())
    out=[];principals=[]
    for row in source['groups']:
        c=principal_certificate(row,allow_proper=row['label']=='Fi23');principals.append(c)
        extra=correction['invariant_vector'] if row['label']=='Fi23' else None
        result=analyze_group(row,extra_vector=extra)
        assert result['principal_block_size']==c['principal_dimension']
        fs=result['rational_wedderburn_factors']
        assert sum(f['center_degree']*f['matrix_size_over_center']**2 for f in fs)==result['algebra_dimension']
        assert sum(f['center_degree']*f['module_multiplicity']**2 for f in fs)==result['commutant_dimension']
        assert sum(f['center_degree']*f['matrix_size_over_center']*f['module_multiplicity'] for f in fs)==row['class_count']
        if extra is not None:result['restoration_correction']='Fi23 principal rank 91, not 92; exact extra invariant vector retained'
        out.append(result)
        (ROOT/'data/census.json').write_text(json.dumps(dict(schema='eidetic-restored-census-v1',software=source['software'],complete=len(out)==26,groups=out),indent=2)+'\n')
        print(row['label'],result['algebra_dimension'],result['center_dimension'],result['commutant_dimension'],flush=True)
    (ROOT/'data/principal_certificates.json').write_text(json.dumps(principals,indent=2)+'\n')
if __name__=='__main__':main()
