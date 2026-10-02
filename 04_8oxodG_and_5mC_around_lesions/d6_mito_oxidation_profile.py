#!/usr/bin/env python3
"""8-oxo-dG calls and opportunity along the mitochondrial genome (Figure 8b).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import json
from collections import defaultdict

import numpy as np
import pysam

P = '/data/8oxo_project'
SRC = f'{P}/allctx_5mC_test/prod_AdMSC'
OUT = f'{P}/AdMSC_paper/sessions/2026-08-14/results'
CHRM, MTLEN = 'CP068254.1', 16569
CUTS = [0.70, 0.50]
COLS = ['locked'] + [str(c) for c in CUTS]
BIN = 100

_log = open(f'{OUT}/d6_run.log', 'w')

def log(m=''):
    print(m, flush=True)
    _log.write(str(m) + '\n')
    _log.flush()

_TJ = json.load(open(f'{P}/refs/thresholds.json'))
TH = {k: v for k, v in _TJ['thresholds'].items() if k not in set(_TJ['high_fp_kmers'])}
assert len(TH) == 110
log(f'[1] allowlist {len(TH)} 5-mers; chrM = {CHRM}, {MTLEN:,} bp; bins {BIN} bp')

feats = []
for line in open(f'{P}/refs/T2T-CHM13v2.0.CAT_liftoff.gff3'):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != 'chrM' or f[2] != 'gene':
        continue
    name = ''
    for kv in f[8].split(';'):
        if kv.startswith('gene_name='):
            name = kv[10:]
    feats.append((int(f[3]), int(f[4]), f[6], name))
feats.sort()
covered = np.zeros(MTLEN + 1, dtype=bool)
for s, e, _, _ in feats:
    covered[s:e + 1] = True
gaps, i = [], 1
while i <= MTLEN:
    if not covered[i]:
        j = i
        while j <= MTLEN and not covered[j]:
            j += 1
        if j - i >= 50:
            gaps.append((i, j - 1, '.', 'non-genic (control region / spacers)'))
        i = j
    else:
        i += 1
feats += gaps
feats.sort()
log(f'[2] chrM features: {len(feats)} '
    f'({sum(1 for f in feats if f[3].startswith("MT-"))} MT- genes, '
    f'{len(gaps)} non-genic blocks >= 50 bp)')

res = {}
for arm in ('NT', 'PA'):
    log(f'\n═══════════════ {arm} ═══════════════')

    pos_by_read = defaultdict(list)
    n_rows = n_call = 0
    with open(f'{SRC}/{arm}/modcall/batch0.txt') as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            n_rows += 1
            t = TH.get(f[3].upper())
            if t is None:
                continue
            n_call += 1
            s = float(f[2])
            pos_by_read[f[0]].append((int(f[1]), s, s >= t))
    log(f'[3] modcall rows {n_rows:,} -> allowlist callable-G {n_call:,} on {len(pos_by_read):,} reads')

    guide_chrm = set()
    g = pysam.AlignmentFile(f'{SRC}/{arm}/guide_aln.bam', 'rb')
    for a in g.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        if a.reference_name == CHRM:
            guide_chrm.add(a.query_name)
    g.close()
    vN = sum(len(v) for r, v in pos_by_read.items() if r in guide_chrm)
    vH = sum(1 for r, v in pos_by_read.items() if r in guide_chrm for _, _, lk in v if lk)
    log(f'[4] VALIDATION, read space, chrM-primary in guide_aln.bam: '
        f'callable-G {vN:,} - hi-conf {vH}   '
        f'(LOCKED {"L1 = 356,221 / 34" if arm == "NT" else "L2 = 231,114 / 28"})')

    cov = {c: np.zeros(MTLEN + 2, dtype=np.int64) for c in ['callable'] + COLS}
    cov_str = {s: np.zeros(MTLEN + 2, dtype=np.int64) for s in ('+', '-')}
    locked_calls = []
    n_aln = n_proj = n_lost = n_rev = 0
    b = pysam.AlignmentFile(f'{SRC}/{arm}/esox_wg.bam', 'rb')
    for a in b.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        pl = pos_by_read.get(a.query_name)
        if pl is None or a.reference_name != CHRM:
            continue
        n_aln += 1
        rev, L = a.is_reverse, a.infer_read_length()
        n_rev += rev
        pairs = dict(a.get_aligned_pairs(matches_only=True))
        for p, sc, lk in pl:
            q = (L - 1 - p) if rev else p
            r = pairs.get(q)
            if r is None:
                n_lost += 1
                continue
            n_proj += 1
            r1 = r + 1
            cov['callable'][r1] += 1
            cov_str['-' if rev else '+'][r1] += 1
            for ci, fl in enumerate([lk] + [sc >= c for c in CUTS]):
                if fl:
                    cov[COLS[ci]][r1] += 1
            if lk:
                locked_calls.append((r1, '-' if rev else '+', sc,
                                     next((f[3] for f in feats if f[0] <= r1 <= f[1]), '.'),
                                     a.query_name))
    b.close()
    log(f'[5] esox reads with a chrM primary alignment: {n_aln:,} ({n_rev:,} reverse -> q = L-1-p)')
    log(f'[5] projected {n_proj:,} allowlist G onto chrM, lost {n_lost:,} '
        f'({100*n_lost/max(1,n_proj+n_lost):.2f}%) - the map is this set')
    log(f'[5] locked-point calls placed on chrM: {len(locked_calls)}   '
        f'at >0.70: {int(cov["0.7"].sum())}   at >0.50: {int(cov["0.5"].sum())}')
    depth = cov['callable'][1:MTLEN + 1]
    log(f'    callable-G depth per position: median {int(np.median(depth))}, '
        f'max {int(depth.max())}, positions with zero {int((depth == 0).sum()):,} of {MTLEN:,}')
    log(f'    strand split of the callable-G: + {int(cov_str["+"].sum()):,}  '
        f'- {int(cov_str["-"].sum()):,}')
    res[arm] = (cov, cov_str, locked_calls, vN, vH)

with open(f'{OUT}/d6_mito_bins.tsv', 'w') as o:
    o.write('bin_start\tbin_end\t' + '\t'.join(
        f'{k}_{a}' for a in ('NT', 'PA')
        for k in ['callable_G', 'calls_locked', 'calls_0.70', 'calls_0.50',
                  'callable_plus', 'callable_minus']) + '\n')
    for s in range(1, MTLEN + 1, BIN):
        e = min(s + BIN - 1, MTLEN)
        row = [s, e]
        for arm in ('NT', 'PA'):
            cov, cov_str, _, _, _ = res[arm]
            row += [int(cov['callable'][s:e + 1].sum()),
                    int(cov['locked'][s:e + 1].sum()),
                    int(cov['0.7'][s:e + 1].sum()),
                    int(cov['0.5'][s:e + 1].sum()),
                    int(cov_str['+'][s:e + 1].sum()),
                    int(cov_str['-'][s:e + 1].sum())]
        o.write('\t'.join(map(str, row)) + '\n')

with open(f'{OUT}/d6_mito_features.tsv', 'w') as o:
    o.write('feature\tstart\tend\tstrand\tlength_bp\t' + '\t'.join(
        f'{k}_{a}' for a in ('NT', 'PA')
        for k in ['callable_G', 'calls_locked', 'calls_0.70', 'calls_0.50']) + '\n')
    for s, e, st, name in feats:
        row = [name, s, e, st, e - s + 1]
        for arm in ('NT', 'PA'):
            cov = res[arm][0]
            row += [int(cov['callable'][s:e + 1].sum()), int(cov['locked'][s:e + 1].sum()),
                    int(cov['0.7'][s:e + 1].sum()), int(cov['0.5'][s:e + 1].sum())]
        o.write('\t'.join(map(str, row)) + '\n')

with open(f'{OUT}/d6_mito_calls_locked.tsv', 'w') as o:
    o.write('arm\tchrM_pos\tread_orientation\toxog_score\tfeature\tread_id\n')
    for arm in ('NT', 'PA'):
        for r1, st, sc, feat, rid in sorted(res[arm][2]):
            o.write(f'{arm}\t{r1}\t{st}\t{sc:.4f}\t{feat}\t{rid}\n')

log('\n═══════════ where the locked-point calls sit ═══════════')
for arm in ('NT', 'PA'):
    fc = defaultdict(int)
    for r1, st, sc, feat, rid in res[arm][2]:
        fc[feat] += 1
    log(f'  {arm} ({len(res[arm][2])} calls): ' +
        ', '.join(f'{k} {v}' for k, v in sorted(fc.items(), key=lambda kv: -kv[1])))

log('\n═══════════ callable-G opportunity, by feature (the denominator of any map) ═══════════')
log(f'  {"feature":<38}{"bp":>7}{"callable-G NT":>15}{"callable-G PA":>15}'
    f'{"calls NT/PA @locked":>22}')
for s, e, st, name in feats:
    nN = int(res['NT'][0]['callable'][s:e + 1].sum())
    nP = int(res['PA'][0]['callable'][s:e + 1].sum())
    cN = int(res['NT'][0]['locked'][s:e + 1].sum())
    cP = int(res['PA'][0]['locked'][s:e + 1].sum())
    if nN + nP == 0:
        continue
    log(f'  {name[:37]:<38}{e-s+1:>7}{nN:>15,}{nP:>15,}{f"{cN}/{cP}":>22}')

log('\nwrote d6_mito_bins.tsv - d6_mito_features.tsv - d6_mito_calls_locked.tsv')
log('\n!! THIS IS A MAP, NOT A RATE. RECIPE 12 section 3 permits projection for maps and figures and')
log('   forbids it for a rate. The reportable mitochondrial rates are the LOCKED rows L1, L2 and L3,')
log('   which are read-space and are cited, never recomputed. the detectability gate applies to them: neither arm')
log('   clears the detectability gate. 34 and 28 events over 16,569 bp cannot support a per-feature')
log('   or per-position contrast, and none is offered.')
_log.close()
