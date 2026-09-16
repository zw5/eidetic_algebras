"""Rebuild all restored computational outputs and run tests.

Default uses the frozen exact GAP exports. --gap PATH refreshes external
inputs before computation. Nothing accesses external resources.
"""
import argparse,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(*args):
    print('+',*args,flush=True)
    subprocess.run(list(args),cwd=ROOT,check=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--gap',type=Path);ap.add_argument('--skip-tests',action='store_true');args=ap.parse_args()
    if args.gap:
        gap=str(args.gap.resolve())
        for name,script,out in [('sporadic_adams_star_census_20260724','export_sporadic_tables.g','sporadic_power_maps.json'),('finite_group_adams_star_census_20260724','export_gap_groups.g','gap_groups.json')]:
            base=ROOT/'proof_construction_assembled/experiments'/name
            p=subprocess.run([gap,'-A','-q',str(base/script)],text=True,capture_output=True,check=True)
            if p.stderr:raise RuntimeError(p.stderr)
            data=json.loads(p.stdout);(base/'output'/out).write_text(json.dumps(data,indent=2)+'\n')
        run(sys.executable,'export_enriched.py','--gap',gap)
        p=subprocess.run([gap,'-A','-q',str(ROOT/'export_psl.g')],text=True,capture_output=True,check=True)
        if p.stderr:raise RuntimeError(p.stderr)
        (ROOT/'data/psl_tables.json').write_text(json.dumps(json.loads(p.stdout),indent=2)+'\n')
    for script in ['investigate_fi23.py','restore_census.py','characters.py','local_factors.py','monster_certificate.py','dynamics.py','families.py','psl_examples.py']:
        run(sys.executable,script)
    run(sys.executable,'-m','proof_construction_assembled.experiments.finite_group_adams_star_census_20260724.run_census','--reuse-gap-data')
    run(sys.executable,'compare_paper.py')
    if not args.skip_tests:run(sys.executable,'-m','pytest','-q')
    print('Computational reproduction complete. Read PAPER_DISCREPANCIES.md before using the current manuscript.')
if __name__=='__main__':main()
