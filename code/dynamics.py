"""Verify the unlabelled prime-power cycle packets in Appendix A."""
from pathlib import Path
from collections import Counter
import json
from eidetic import cycle_lengths
ROOT=Path(__file__).resolve().parent

def main():
    rows=json.loads((ROOT/'data/sporadic_tables.json').read_text())['groups']
    chars={r['label']:r for r in json.loads((ROOT/'data/good_label_characters.json').read_text())['groups']}
    results=[];packets=[]
    for row in rows:
        directions=[];support=set();packet=[]
        for n,pm in zip(row['directions'],row['power_maps']):
            counts=Counter(cycle_lengths([i-1 for i in pm]));support.update(counts)
            packet.append(tuple(counts.get(i,0) for i in (1,2,3,5)))
            directions.append(dict(exponent=n,cycle_counts=dict(counts)))
        orders={1}|{b['character_order'] for b in chars[row['label']]['blocks']}
        assert support==orders and support<={1,2,3,5}
        packet=tuple(sorted(packet));packets.append(packet)
        results.append(dict(label=row['label'],directions=directions,packet=packet,cycle_support=sorted(support)))
    assert len(set(packets))==26
    assert sum(len(r['directions']) for r in rows)==258
    result=dict(schema='eidetic-dynamics-v1',groups=results,direction_count=258,distinct_packet_count=26,
                cubic_groups=[r['label'] for r in results if 3 in r['cycle_support']],
                quintic_groups=[r['label'] for r in results if 5 in r['cycle_support']])
    (ROOT/'data/dynamics.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Verified 258 maps; 26 distinct unlabelled packets')
if __name__=='__main__':main()
