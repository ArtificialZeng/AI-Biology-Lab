"""Exercise released Biopython before/after a temporary two-line-per-API repair."""
import json
import math
from Bio.Align import Alignment
from Bio.Align.analysis import calculate_dn_ds
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from independent_oracle import estimate, fixtures


def rounded(value):
    if isinstance(value, float):
        return round(value, 12)
    if isinstance(value, dict):
        return {str(k): rounded(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [rounded(v) for v in value]
    return value


rows = []
for case in fixtures():
    oracle = estimate(case['sequence1'], case['sequence2'])
    results = []
    for first, second in [(case['sequence1'], case['sequence2']),
                          (case['sequence2'], case['sequence1'])]:
        for api in ['Align', 'codonalign']:
            actual = (calculate_dn_ds(Alignment([first, second]), method='LWL85')
                      if api == 'Align' else
                      cal_dn_ds(CodonSeq(first), CodonSeq(second), method='LWL85'))
            results.append({'api': api, 'dN': float(actual[0]), 'dS': float(actual[1]),
                            'agrees': all(math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)
                                          for a, b in zip(actual, [oracle['dN'], oracle['dS']]))})
    rows.append({'name': case['name'], 'oracle': oracle, 'results': results})
print(json.dumps(rounded(rows), sort_keys=True, indent=2))
