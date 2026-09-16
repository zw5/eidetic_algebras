"""Compare current paper's printed sporadic dimensions with recomputed values."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent

def main():
    source=ROOT.parent/'paper/sections/sporadic.tex'
    if not source.exists():
        print('Manuscript source absent; frozen comparison report retained.');return
    expected=[tuple(map(int,m)) for m in re.findall(r'& (\d+) & (\d+) & (\d+) & (\d+) &',source.read_text())]
    rows=json.loads((ROOT/'data/census.json').read_text())['groups'];assert len(expected)==len(rows)==26
    result=[]
    for r,e in zip(rows,expected):
        actual=tuple(r[k] for k in ['class_count','algebra_dimension','center_dimension','commutant_dimension'])
        result.append(dict(group=r['label'],paper=list(e),computed=list(actual),matches=actual==e))
    output=dict(columns=['class_count','algebra_dimension','center_dimension','commutant_dimension'],
                all_rows_match=all(r['matches'] for r in result),rows=result)
    (ROOT/'data/paper_comparison.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Paper comparison:',sum(r['matches'] for r in result),'/ 26 matching rows')
    for r in result:
        if not r['matches']:print('DISCREPANCY:',r)
if __name__=='__main__':main()
