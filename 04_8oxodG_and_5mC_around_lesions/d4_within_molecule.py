#!/usr/bin/env python3
"""Ordinary-guanine anchors on lesion-carrying molecules, derived exactly from the anchor classes.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv
from collections import defaultdict

OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-14/results'
BANDS = [(11, 100), (101, 250), (251, 500), (501, 1000)]
BIN = 25

_log = open(f'{OUT}/d4_run.log', 'w')

def log(m=''):
    print(m, flush=True)
    _log.write(str(m) + '\n')
    _log.flush()

rows = list(csv.DictReader(open(f'{OUT}/d3_profile_bin1.tsv'), delimiter='\t'))
log(f'[1] read {len(rows):,} per-base-pair rows from d3_profile_bin1.tsv')

cell = defaultdict(lambda: [0.0, 0.0])
for r in rows:
    d = int(r['dist_lo'])
    n = float(r['n_pairs'])
    cell[(r['arm'], r['cut'], r['class'], d)][0] += float(r['mean_prediction_score']) * n
    cell[(r['arm'], r['cut'], r['class'], d)][1] += n

arms = sorted({r['arm'] for r in rows})
cuts = sorted({r['cut'] for r in rows}, key=lambda c: (c != 'locked', c))
maxd = max(int(r['dist_lo']) for r in rows)

def get(arm, cut, cls, lo, hi):
    s = n = 0.0
    for d in range(lo, hi + 1):
        a, b = cell[(arm, cut, cls, d)]
        s += a
        n += b
    return s, n

def same(arm, cut, lo, hi):
    sa, na = get(arm, cut, 'canon_all', lo, hi)
    sf, nf = get(arm, cut, 'canon_free', lo, hi)
    return sa - sf, na - nf

with open(f'{OUT}/d4_profile_4class.tsv', 'w') as o:
    o.write('arm\tcut\tdist_lo\tdist_hi\t'
            'n_oxo\tmean_oxo\tn_canon_same\tmean_canon_same\t'
            'n_canon_free\tmean_canon_free\tn_canon_all\tmean_canon_all\t'
            'distance_term\tcomposition_term\n')
    for arm in arms:
        for cut in cuts:
            for lo in range(1, maxd + 1, BIN):
                hi = min(lo + BIN - 1, maxd)
                so, no = get(arm, cut, 'oxo', lo, hi)
                ss, ns = same(arm, cut, lo, hi)
                sf, nf = get(arm, cut, 'canon_free', lo, hi)
                sa, na = get(arm, cut, 'canon_all', lo, hi)
                if no == 0 and na == 0:
                    continue
                mo = so / no if no else float('nan')
                ms = ss / ns if ns else float('nan')
                mf = sf / nf if nf else float('nan')
                ma = sa / na if na else float('nan')
                o.write(f'{arm}\t{cut}\t{lo}\t{hi}\t{int(no)}\t{mo:.6f}\t{int(ns)}\t{ms:.6f}\t'
                        f'{int(nf)}\t{mf:.6f}\t{int(na)}\t{ma:.6f}\t'
                        f'{mo-ms:.6f}\t{ms-mf:.6f}\n')
log(f'[2] wrote d4_profile_4class.tsv  ({BIN} bp bins)')

with open(f'{OUT}/d4_bands.tsv', 'w') as o:
    o.write('arm\tcut\tband_lo\tband_hi\tn_oxo\tmean_oxo\tn_canon_same\tmean_canon_same\t'
            'n_canon_free\tmean_canon_free\tdistance_term\tcomposition_term\n')
    for arm in arms:
        for cut in cuts:
            log(f'\n  {arm}, cut {cut}')
            log(f'    {"band (bp)":>12}{"oxo":>19}{"canon SAME read":>23}{"canon FREE read":>23}'
                f'{"distance":>11}{"composition":>13}')
            for lo, hi in BANDS:
                so, no = get(arm, cut, 'oxo', lo, hi)
                ss, ns = same(arm, cut, lo, hi)
                sf, nf = get(arm, cut, 'canon_free', lo, hi)
                mo = so / no if no else float('nan')
                ms = ss / ns if ns else float('nan')
                mf = sf / nf if nf else float('nan')
                log(f'    {f"{lo}-{hi}":>12}{f"{mo:.4f} ({int(no)})":>19}'
                    f'{f"{ms:.4f} ({int(ns)})":>23}{f"{mf:.4f} ({int(nf)})":>23}'
                    f'{mo-ms:>+11.4f}{ms-mf:>+13.4f}')
                o.write(f'{arm}\t{cut}\t{lo}\t{hi}\t{int(no)}\t{mo:.6f}\t{int(ns)}\t{ms:.6f}\t'
                        f'{int(nf)}\t{mf:.6f}\t{mo-ms:.6f}\t{ms-mf:.6f}\n')
log('\n[3] wrote d4_bands.tsv')

log('\n  distance term    = oxo - canon_SAME   : within one molecule, is 5mC lower near the lesion')
log('  composition term = canon_SAME - canon_FREE : between molecules, are lesion-carrying')
log('                     molecules globally less methylated')
log('\n!! no threshold and no statistics anywhere - the paper reports a curve, and so does this.')
log('   The pairs behind each mean are not independent (many CpG per anchor, many anchors per read);')
log('   that is why no interval is quoted, not an oversight.')
_log.close()
