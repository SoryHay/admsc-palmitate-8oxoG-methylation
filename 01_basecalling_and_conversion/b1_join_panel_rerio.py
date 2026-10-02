#!/usr/bin/env python3
"""Join Rerio CpG calls to the reads carrying esox calls.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os
import sys
import time

P = '/data/8oxo_project'
OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-13/results'
os.makedirs(OUT, exist_ok=True)

def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)

def panel_ids(arm):
    """Unique read ids across BOTH modcall batches. batch1 was missed on first inspection."""
    ids = set()
    for b in (0, 1):
        p = f'{P}/panel_oxo/modcall_{arm}/batch{b}.txt'
        if not os.path.isfile(p):
            log(f'  WARNING: {p} absent')
            continue
        with open(p) as fh:
            next(fh, None)
            for line in fh:
                i = line.find('\t')
                if i > 0:
                    ids.add(line[:i])
    return ids

gate = []
for arm in ('NT', 'PA'):
    ids = panel_ids(arm)
    log(f'{arm}: {len(ids):,} unique panel read ids from modcall batch0+batch1')

    src = f'{P}/AdMSC_DMR/readlevel_calls/{arm}.calls.tsv'
    dst = f'{OUT}/panel_rerio_calls_{arm}.tsv'
    seen = set()
    n_in = n_out = 0
    t0 = time.time()
    with open(src) as fh, open(dst, 'w') as o:
        for line in fh:
            n_in += 1
            i = line.find('\t')
            if i > 0 and line[:i] in ids:
                o.write(line)
                n_out += 1
                seen.add(line[:i])
            if n_in % 20_000_000 == 0:
                log(f'  {arm}: {n_in/1e6:.0f} M rows, {len(seen):,} panel reads hit, '
                    f'{time.time()-t0:.0f} s')
    log(f'{arm}: DONE {n_in:,} rows scanned, {n_out:,} calls kept, '
        f'{len(seen):,} of {len(ids):,} panel reads found '
        f'({100*len(seen)/len(ids):.1f} %), {time.time()-t0:.0f} s')
    gate.append((arm, len(ids), len(seen), n_out, n_in))

with open(f'{OUT}/b1_gate.tsv', 'w') as o:
    o.write('arm\tpanel_read_ids\tfound_in_readlevel\tpct_found\trerio_calls_kept\trows_scanned\n')
    for arm, nid, nfound, nout, nin in gate:
        o.write(f'{arm}\t{nid}\t{nfound}\t{100*nfound/nid:.2f}\t{nout}\t{nin}\n')
log('gate written to results/b1_gate.tsv')
