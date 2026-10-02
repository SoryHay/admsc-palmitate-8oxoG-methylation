#!/usr/bin/env python3
"""Remaining tables of the public Zenodo record.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv, os, shutil
import pandas as pd
from docx import Document

A = '/data/8oxo_project/AdMSC_paper'
S = A + '/sessions'
OUT = A + '/deposit_zenodo_profiles_v1/data'

rows = [
    ('L1', 'mitochondrial DNA', 'NT', 34, 356221, 95.45, '66.1-133.4'), ('L2', 'mitochondrial DNA', 'PA', 28, 231114, 121.15, '80.5-175.2'),
    ('L4', 'nuclear 82-gene panel', 'NT', 48, 772996, 62.1, ''), ('L5', 'nuclear 82-gene panel', 'PA', 45, 911550, 49.4, ''),
    ('L14', 'ribosomal DNA (RepeatMasker rRNA class)', 'NT', 75, 648674, 115.6, '90.9-144.9'),
    ('L15', 'ribosomal DNA (RepeatMasker rRNA class)', 'PA', 80, 830815, 96.3, '76.4-119.8'),
    ('L17', 'ribosomal DNA, five acrocentric arrays', 'NT', 43, 307257, 139.9, '101.3-188.5'),
    ('L18', 'ribosomal DNA, five acrocentric arrays', 'PA', 41, 404566, 101.3, '72.7-137.5')]
pd.DataFrame(rows, columns=['locked_row', 'compartment', 'arm', 'calls_8oxodG', 'callable_G', 'rate_per_million_G', 'CI95_per_million'])\
  .to_csv(f'{OUT}/oxidation_8oxodG_rates_by_compartment.tsv', sep='\t', index=False)
pd.DataFrame([('L3', 'mitochondrial DNA', 1.27), ('L6', 'nuclear 82-gene panel', 0.80), ('L16', 'ribosomal DNA (rRNA class)', 0.83),
              ('L19', 'ribosomal DNA, acrocentric arrays', 0.72)], columns=['locked_row', 'compartment', 'ratio_PA_NT'])\
  .to_csv(f'{OUT}/oxidation_8oxodG_ratios_by_compartment.tsv', sep='\t', index=False)

COPY = {f'{S}/2026-08-14/results/f5_oxo_context.tsv': 'oxidation_8oxodG_by_sequence_context.tsv',
        f'{S}/2026-08-14/results/f9_oxo_by_methylation_stratum.tsv': 'oxidation_8oxodG_by_molecule_methylation_tertile.tsv',
        f'{S}/2026-08-14/results/f2_per_read_beta.tsv': 'methylation_per_read_beta_distribution.tsv',
        f'{S}/2026-08-14/results/f6_repeat_classes.tsv': 'methylation_by_repeat_class.tsv',
        f'{S}/2026-08-14/results/f7_metagene.tsv': 'methylation_metagene.tsv',
        f'{S}/2026-09-03/results/fig_P1_5hmC_v5_values.tsv': 'methylation_5mC_5hmC_genome_wide.tsv'}
for src, dst in COPY.items():
    shutil.copyfile(src, f'{OUT}/{dst}')
d = pd.read_csv(f'{S}/2026-08-14/results/d11_inflammasome_axis.tsv', sep='\t').drop(columns=['NT_calls', 'PA_calls'])
d.to_csv(f'{OUT}/methylation_IL1B_production_axis.tsv', sep='\t', index=False)

src = f'{S}/2026-08-26b_manuscript/results/d3_profile_bin50_COPY_from_2026-08-14.tsv'
r = pd.read_csv(src, sep='\t')
out = []
for (arm, cut), g in r.groupby(['arm', 'cut']):
    a = g[g['class'] == 'canon_all'].set_index('dist_lo'); f = g[g['class'] == 'canon_free'].set_index('dist_lo')
    o = g[g['class'] == 'oxo'].set_index('dist_lo')
    for lo in sorted(a.index):
        na, nf = a.loc[lo, 'n_pairs'], f.loc[lo, 'n_pairs']
        same = (a.loc[lo, 'mean_prediction_score'] * na - f.loc[lo, 'mean_prediction_score'] * nf) / (na - nf)
        for cls, n, m in (('8-oxo-dG anchors', o.loc[lo, 'n_pairs'], o.loc[lo, 'mean_prediction_score']),
                          ('ordinary G, lesion-carrying molecules', na - nf, same),
                          ('ordinary G, lesion-free molecules', nf, f.loc[lo, 'mean_prediction_score'])):
            out.append((arm, cut, lo, int(a.loc[lo, 'dist_hi']), cls, int(n), round(float(m), 6)))
pd.DataFrame(out, columns=['arm', 'cut', 'dist_lo_bp', 'dist_hi_bp', 'anchor_class', 'anchor_CpG_pairs', 'mean_5mC_prediction_score'])\
  .to_csv(f'{OUT}/methylation_around_8oxodG_profiles.tsv', sep='\t', index=False)

sup = Document(f'{A}/Supplementary_Data_v11-3.docx')
for ti, name in ((1, 'table_S2_sequencing_run_summary.tsv'), (2, 'table_S3_5mC_5hmC_by_read_quality.tsv'),
                 (3, 'table_S4_5mC_5hmC_by_quality_decile.tsv'), (4, 'table_S5_5mC_5hmC_ribosomal_DNA.tsv')):
    with open(f'{OUT}/{name}', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        for row in sup.tables[ti].rows:
            w.writerow([c.text.strip() for c in row.cells])
print(sorted(os.listdir(OUT)))
