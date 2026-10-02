#!/usr/bin/env python3
"""Region-level 5mC profiles of the public Zenodo record.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import gzip, os, sys
import numpy as np, pandas as pd, pysam

P = '/data/8oxo_project'
CTX = P + '/AdMSC_DMR/trackB_rerio_allctx/{arm}/ctx_thr070.cpg.tsv.gz'
FEAT = P + '/AdMSC_DMR/features'
REF = P + '/refs/T2T-CHM13v2.0.fa'
OUT = P + '/AdMSC_paper/deposit_zenodo_profiles_v1/data'
os.makedirs(OUT, exist_ok=True)
MAP = {r.split('\t')[1]: r.split('\t')[0] for r in open(P + '/refs/chrom_map.tsv') if r.strip() and not r.startswith('#')}
EXCLUDE = {'CP068255.2', 'CP086569.2', 'CP068254.1'}
MIN_CPG, MIN_OBS = 5, 20
M006 = {'islands': 0.8894, 'exons': 0.8945, 'shelves': 0.8916, 'shores': 0.8967, 'gene bodies': 0.8971, 'promoters': 0.9130}

idx = {}
tot = {}
for arm in ('NT', 'PA'):
    df = pd.read_csv(CTX.format(arm=arm), sep='\t', usecols=['chrom', 'start', 'Nvalid', 'Nmod'],
                     dtype={'chrom': str, 'start': np.int64, 'Nvalid': np.int64, 'Nmod': np.int64})
    tot[arm] = (int(df.Nmod.sum()), int(df.Nvalid.sum()))
    for c, g in df.groupby('chrom', sort=False):
        g = g.sort_values('start')
        idx[(arm, c)] = (g.start.values, np.concatenate([[0], np.cumsum(g.Nvalid.values)]),
                         np.concatenate([[0], np.cumsum(g.Nmod.values)]))
    print(f'{arm}: {len(df):,} CpG rows, genome-wide 5mC {100 * tot[arm][0] / tot[arm][1]:.4f} %', flush=True)
    del df

def sums(arm, c, s, e):
    if (arm, c) not in idx: return 0, 0
    pos, cv, cm = idx[(arm, c)]
    i, j = np.searchsorted(pos, s, 'left'), np.searchsorted(pos, e, 'left')
    return int(cm[j] - cm[i]), int(cv[j] - cv[i])

fa = pysam.FastaFile(REF)
def ref_cpg(c, s, e):
    return fa.fetch(c, s, e).upper().count('CG')

def read_bed(name):
    rows = []
    for ln in open(f'{FEAT}/{name}.bed'):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 3 or ln.startswith(('#', 'track')): continue
        rows.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else f'{f[0]}:{f[1]}-{f[2]}', f[4] if len(f) > 4 else ''))
    return rows

def profile(rows, label, extra_col=None):
    """Per-region table (autosomes only) + class sums over all regions (every chromosome)."""
    out, csum = [], {'NT': [0, 0], 'PA': [0, 0]}
    for c, s, e, name, ex in rows:
        v = {}
        for arm in ('NT', 'PA'):
            v[arm] = sums(arm, c, s, e)
        if all(v[a][1] >= MIN_OBS for a in v):
            for a in v: csum[a][0] += v[a][0]; csum[a][1] += v[a][1]
        if c in EXCLUDE or c not in MAP: continue
        k = ref_cpg(c, s, e)
        ok = k >= MIN_CPG and all(v[a][1] >= MIN_OBS for a in v)
        rec = {'region_id': name, 'chrom_T2T': MAP[c], 'accession': c, 'start': s, 'end': e, 'n_CpG_reference': k,
               'pct_5mC_NT': round(100 * v['NT'][0] / v['NT'][1], 2) if ok else None,
               'pct_5mC_PA': round(100 * v['PA'][0] / v['PA'][1], 2) if ok else None}
        if extra_col: rec[extra_col] = ex
        out.append(rec)
    d = pd.DataFrame(out)
    d['ratio_PA_NT'] = (d.pct_5mC_PA / d.pct_5mC_NT).round(4)
    print(f'{label}: {len(rows):,} regions, {len(d):,} autosomal, {d.pct_5mC_NT.notna().sum():,} reported', flush=True)
    return d, csum

summary = []
def add_summary(label, csum):
    nt, pa = csum['NT'][0] / csum['NT'][1], csum['PA'][0] / csum['PA'][1]
    summary.append({'region_class': label, 'pct_5mC_NT': round(100 * nt, 4), 'pct_5mC_PA': round(100 * pa, 4),
                    'ratio_PA_NT': round(pa / nt, 5)})

CLASSES = [('islands', 'islands', 'cpg_islands'), ('shores', 'shores', 'cpg_shores'), ('shelves', 'shelves', 'cpg_shelves'),
           ('promoters_tss1kb', 'promoters', 'promoters_tss1kb'), ('exons_merged', 'exons', 'exons'),
           ('genebodies', 'gene bodies', 'gene_bodies')]
for bed, label, fname in CLASSES:
    d, cs = profile(read_bed(bed), label)
    add_summary(label, cs)
    if fname in ('cpg_islands', 'promoters_tss1kb', 'gene_bodies'):
        d.to_csv(f'{OUT}/methylation_profile_{fname}.tsv.gz', sep='\t', index=False, na_rep='NA', compression='gzip')
nt, pa = tot['NT'][0] / tot['NT'][1], tot['PA'][0] / tot['PA'][1]
summary.append({'region_class': 'whole genome', 'pct_5mC_NT': round(100 * nt, 4), 'pct_5mC_PA': round(100 * pa, 4),
                'ratio_PA_NT': round(pa / nt, 5)})
S = pd.DataFrame(summary)
print(S.to_string(index=False))
bad = [(r.region_class, r.ratio_PA_NT) for r in S.itertuples() if r.region_class in M006 and round(r.ratio_PA_NT, 4) != M006[r.region_class]]
if bad or abs(pa / nt - 0.89845) > 1e-5:
    sys.exit(f'STOP: class summary does not reproduce the published class ratios: {bad}, genome {pa / nt:.5f}')
S.to_csv(f'{OUT}/methylation_summary_region_classes.tsv', sep='\t', index=False)

panels = []
for bed, panel in [('adipokine_bygene', 'adipokine'), ('oxidative_bygene', 'oxidative_damage_redox'),
                   ('mito_metabolic_bygene', 'mitochondrial_metabolic'), ('senescence_bygene', 'ageing_senescence'),
                   ('innate_genebodies', 'innate_immunity')]:
    d, _ = profile(read_bed(bed), panel, extra_col='annotation')
    d.insert(0, 'panel', panel); panels.append(d)
pd.concat(panels).to_csv(f'{OUT}/methylation_profile_gene_panels.tsv', sep='\t', index=False, na_rep='NA')
print('wrote', sorted(os.listdir(OUT)))
