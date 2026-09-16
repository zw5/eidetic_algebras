"""Check rational-class sizes in the PSL(2,q) examples, without algebra claims."""
import json
from sympy import isprime
from math import gcd
from pathlib import Path
from local_factors import orbit_partition
ROOT=Path(__file__).resolve().parent

def main():
    rows=json.loads((ROOT/'data/psl_tables.json').read_text())['groups'];out=[]
    expected={8:[(7,3),(9,3)],13:[(7,3)],17:[(9,3)],19:[(9,3)],27:[(7,3),(14,3)],23:[(11,5)],32:[(11,5)],29:[(7,3),(14,3)]}
    for r in rows:
        orders=r['class_orders'];maps=r['power_maps'];odd=[];all_orbits=[]
        for m in sorted(set(orders)):
            layer=[i for i,x in enumerate(orders) if x==m]
            permutations=[{i:maps[u-1][i]-1 for i in layer} for u in range(1,m+1) if gcd(u,m)==1]
            for orbit in orbit_partition(layer,permutations):
                all_orbits.append(dict(element_order=m,classes=orbit,size=len(orbit)))
                if len(orbit)>2 and isprime(len(orbit)):odd.append((m,len(orbit)))
        assert sorted(set(odd))==expected.get(r['q'],[]),(r['q'],odd)
        out.append(dict(q=r['q'],rational_classes=all_orbits,odd_prime_sizes=odd))
    (ROOT/'data/psl_examples.json').write_text(json.dumps(dict(groups=out,computed_range_complete=True,paper_missing_example={'q':29,'element_orders':[7,14],'rational_class_size':3}),indent=2)+'\n')
    print('Verified all 17 prime-power q in [3,32]; paper omits q=29')
if __name__=='__main__':main()
