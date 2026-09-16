"""Re-export all 26 tables with names, unit actions, local maps, and table symmetries."""
import argparse,json,subprocess
from pathlib import Path
from symmetry_helpers import _global_unit_generators
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'proof_construction_assembled/experiments/sporadic_adams_star_census_20260724'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--gap',required=True);args=ap.parse_args()
    source=json.loads((BASE/'output/sporadic_power_maps.json').read_text())
    s=(BASE/'export_sporadic_tables.g').read_text()
    units=[[_[0] for _ in _global_unit_generators(row['exponent'])] for row in source['groups']]
    s=s.replace('GroupRows := [','UnitRows := '+str(units)+';;\nGroupRows := [')
    mark='    Print("\\\"power_maps\\\":[");'
    extra='''    Print("\\\"class_names\\\":[");
    names := ClassNames(table);
    for j in [1..Length(names)] do
      if j>1 then Print(","); fi;
      Print("\\\"",names[j],"\\\"");
    od;
    Print("],\\\"unit_generators\\\":",UnitRows[rowIndex],",");
    Print("\\\"unit_power_maps\\\":[");
    for j in [1..Length(UnitRows[rowIndex])] do
      if j>1 then Print(","); fi;
      Print(List([1..NrConjugacyClasses(table)], i -> PowerMap(table,UnitRows[rowIndex][j] mod OrdersClassRepresentatives(table)[i])[i]));
    od;
    Print("],\\\"local_power_maps\\\":[");
    for j in [1..Maximum(OrdersClassRepresentatives(table))] do
      if j>1 then Print(","); fi;
      Print(List(PowerMap(table,j), x -> x));
    od;
    Print("],\\\"table_automorphisms\\\":",
      List(GeneratorsOfGroup(AutomorphismsOfTable(table)),
        a -> List([1..NrConjugacyClasses(table)], i -> i^a)),",");
'''
    assert mark in s
    s=s.replace(mark,extra+mark)
    s=s.replace('for rowIndex in', 'table := fail;;\nfor rowIndex in')
    (ROOT/'data').mkdir(exist_ok=True)
    gapfile=ROOT/'data/export_enriched.g';gapfile.write_text(s)
    completed=subprocess.run([args.gap,'-A','-q',str(gapfile)],text=True,capture_output=True,check=True)
    if completed.stderr:raise RuntimeError(completed.stderr)
    data=json.loads(completed.stdout)
    (ROOT/'data/sporadic_tables.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Exported',len(data['groups']),'tables')
if __name__=='__main__':main()
