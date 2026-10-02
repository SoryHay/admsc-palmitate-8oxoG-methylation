#!/usr/bin/env python3
"""Mean 5mC prediction score of CpGs around 8-oxo-dG and ordinary guanine anchors on the same read.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import json
from collections import defaultdict

import numpy as np
import pysam

P = '/data/8oxo_project'
PREV = f'{P}/AdMSC_paper/sessions/2026-08-13/results'
OUT = f'{P}/AdMSC_paper/sessions/2026-08-14/results'
W = 1000
MASK = 10
CUTS = ['locked', 0.70]

_log = open(f'{OUT}/d3_run.log', 'w')

def log(m=''):
    print(m, flush=True)
    _log.write(str(m) + '\n')
    _log.flush()

_TJ = json.load(open(f'{P}/refs/thresholds.json'))
TH = {k: v for k, v in _TJ['thresholds'].items() if k not in set(_TJ['high_fp_kmers'])}
assert len(TH) == 110
log(f'[1] allowlist {len(TH)} 5-mers ("analysis was restricted to only 5-mers in the training set")')
log(f'[1] window +/-{W} bp, crosstalk mask +/-{MASK} bp applied to BOTH classes (deviation c)')

def cpg_and_score(fields):
    """Return (is_cpg, P_5mC) for one modkit `extract calls` row, or None if it is not a C call.

    Column convention verified 2026-08-14 on the data itself:
      0 read_id . 2 ref_pos . 3 chrom . 5 ref_strand . 12 prob OF THE CALLED STATE . 13 mod_code
      . 15 reference 5-mer, MIXED CASE . 17 canonical base
    A CpG on the reverse strand appears as 'g' at the centre of the reference-forward 5-mer with its
    partner 'c' at index 1; on the forward strand the centre is 'c' with 'g' at index 3.
    """
    if fields[17] != 'C':
        return None
    k = fields[15].lower()
    if len(k) != 5:
        return None
    if fields[5] == '+':
        is_cpg = k[2] == 'c' and k[3] == 'g'
    else:
        is_cpg = k[2] == 'g' and k[1] == 'c'
    q = float(fields[12])
    p5mc = q if fields[13] == 'm' else 1.0 - q
    return is_cpg, p5mc

results = {}
for arm in ('NT', 'PA'):
    log(f'\n═══════════════════ {arm} ═══════════════════')

    calls = defaultdict(list)
    n_rows = n_c = n_cpg = 0
    psum = 0.0
    for line in open(f'{PREV}/panel_rerio_calls_{arm}.tsv'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 18:
            continue
        n_rows += 1
        r = cpg_and_score(f)
        if r is None:
            continue
        n_c += 1
        is_cpg, p5 = r
        if not is_cpg:
            continue
        n_cpg += 1
        psum += p5
        calls[f[0]].append((int(f[2]), p5, f[3]))
    log(f'[2] Rerio calls {n_rows:,} rows -> C calls {n_c:,} -> CpG calls {n_cpg:,} '
        f'on {len(calls):,} reads')
    log(f'    baseline: mean methylation PREDICTION SCORE over all CpG calls = '
        f'{psum/max(1,n_cpg):.4f}   ("average genome-wide methylation level", the paper\'s recovery '
        f'reference)')

    pos_by_read = defaultdict(list)
    for b in (0, 1):
        with open(f'{P}/panel_oxo/modcall_{arm}/batch{b}.txt') as fh:
            next(fh)
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 4:
                    continue
                t = TH.get(f[3].upper())
                if t is None:
                    continue
                s = float(f[2])
                pos_by_read[f[0]].append((int(f[1]), s, s >= t))

    anchors = defaultdict(list)
    n_aln = n_proj = 0
    bam = pysam.AlignmentFile(f'{PREV}/esox_{arm}_t2t.bam', 'rb')
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        pl = pos_by_read.get(a.query_name)
        if pl is None or a.query_name not in calls:
            continue
        n_aln += 1
        rev, L = a.is_reverse, a.infer_read_length()
        pairs = dict(a.get_aligned_pairs(matches_only=True))
        for p, sc, lk in pl:
            q = (L - 1 - p) if rev else p
            r = pairs.get(q)
            if r is None:
                continue
            n_proj += 1
            anchors[a.query_name].append((r, sc, lk, a.reference_name))
    bam.close()
    log(f'[3] reads carrying BOTH an esox alignment and Rerio CpG calls: {n_aln:,}; '
        f'allowlist G projected on them: {n_proj:,}')

    per_cut = {}
    for cut in CUTS:
        acc = {c: (np.zeros(2 * W + 1), np.zeros(2 * W + 1))
               for c in ('oxo', 'canon_all', 'canon_free')}
        n_anchor = defaultdict(int)
        reach = []
        for rid, alist in anchors.items():
            cl = calls.get(rid)
            if not cl:
                continue
            chrom_c = cl[0][2]
            cpos = np.array([x[0] for x in cl], dtype=np.int64)
            cp = np.array([x[1] for x in cl], dtype=np.float64)
            is_oxo = [(lk if cut == 'locked' else sc >= cut) for _, sc, lk, _ in alist]
            read_has_lesion = any(is_oxo)
            for (rpos, sc, lk, chrom_a), oxo in zip(alist, is_oxo):
                if chrom_a != chrom_c:
                    continue
                d = cpos - rpos
                sel = (np.abs(d) > MASK) & (np.abs(d) <= W)
                if not sel.any():
                    continue
                idx = d[sel] + W
                cls = 'oxo' if oxo else 'canon_all'
                np.add.at(acc[cls][0], idx, cp[sel])
                np.add.at(acc[cls][1], idx, 1.0)
                n_anchor[cls] += 1
                if oxo:
                    reach.append(int(np.abs(d[sel]).max()))
                elif not read_has_lesion:
                    np.add.at(acc['canon_free'][0], idx, cp[sel])
                    np.add.at(acc['canon_free'][1], idx, 1.0)
                    n_anchor['canon_free'] += 1
        per_cut[cut] = (acc, dict(n_anchor), reach)
        tot = {c: int(acc[c][1].sum()) for c in acc}
        log(f'[4] cut {str(cut):>6}: anchors  oxo {n_anchor["oxo"]:,}  '
            f'canon_all {n_anchor["canon_all"]:,}  canon_free {n_anchor["canon_free"]:,}   |   '
            f'(anchor,CpG) pairs  oxo {tot["oxo"]:,}  canon_all {tot["canon_all"]:,}  '
            f'canon_free {tot["canon_free"]:,}')
        if reach:
            r = np.array(reach)
            log(f'    achievable distance from an 8-oxo-dG anchor (max |d| per anchor): '
                f'median {int(np.median(r))} bp, p90 {int(np.percentile(r,90))} bp, max {r.max()} bp '
                f'- the read span, not the window, is the limit')
    results[arm] = per_cut

for BIN in (1, 25, 50):
    fn = f'{OUT}/d3_profile_bin{BIN}.tsv'
    with open(fn, 'w') as o:
        o.write('arm\tcut\tclass\tdist_lo\tdist_hi\tn_pairs\tmean_prediction_score\n')
        for arm in ('NT', 'PA'):
            for cut in CUTS:
                acc = results[arm][cut][0]
                for cls in ('oxo', 'canon_all', 'canon_free'):
                    s, n = acc[cls]

                    half_s = s[W + 1:] + s[:W][::-1]
                    half_n = n[W + 1:] + n[:W][::-1]
                    for b0 in range(0, W, BIN):
                        b1 = min(b0 + BIN, W)
                        ss, nn = half_s[b0:b1].sum(), half_n[b0:b1].sum()
                        if nn == 0:
                            continue
                        o.write(f'{arm}\t{cut}\t{cls}\t{b0+1}\t{b1}\t{int(nn)}\t{ss/nn:.6f}\n')
    log(f'wrote {fn}')

log('\n═══════════ the paper\'s quantity, in coarse bands (mean prediction score) ═══════════')
for arm in ('NT', 'PA'):
    for cut in CUTS:
        acc = results[arm][cut][0]
        log(f'\n  {arm}, cut {cut}')
        log(f'    {"band (bp)":>14}' + ''.join(f'{c:>14}' for c in
                                               ('oxo', 'canon_all', 'canon_free')))
        for b0, b1 in ((11, 100), (101, 250), (251, 500), (501, 1000)):
            row = f'    {f"{b0}-{b1}":>14}'
            for cls in ('oxo', 'canon_all', 'canon_free'):
                s, n = acc[cls]
                hs = s[W + 1:] + s[:W][::-1]
                hn = n[W + 1:] + n[:W][::-1]
                ss, nn = hs[b0 - 1:b1].sum(), hn[b0 - 1:b1].sum()
                row += f'{(ss/nn if nn else float("nan")):>9.4f} ({int(nn)})' if nn else f'{"-":>14}'
            log(row)

log('\n!! WHAT THIS IS AND IS NOT')
log('  * it is the paper\'s own quantity: the mean methylation PREDICTION SCORE per base-pair offset,')
log('    with no threshold and no test, as published.')
log('  * it is NOT their experiment: DAAO-induced intranuclear H2O2 against native AdMSC, a 1 kb')
log('    window against their 10 kb, and a different 5mC model. Only the SHAPE is comparable.')
log('  * no absolute oxidation rate, no operating point adopted, no magnitude')
log('    (n = 1 per arm, condition confounded with flow cell).')
_log.close()
