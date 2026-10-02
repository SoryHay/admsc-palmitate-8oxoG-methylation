#!/usr/bin/env python3
"""IL-1beta production-axis genes at five genomic scales against the matched background (Figure 11c).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv

R = '/data/8oxo_project/AdMSC_DMR/results'
OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-14/results'

AXIS = [
    ('--- canonical NLRP3 axis ---', None),
    ('NLRP3', 'sensor, NLRP3 inflammasome'),
    ('NEK7', 'licensing kinase for NLRP3'),
    ('PYCARD', 'adaptor (ASC)'),
    ('CASP1', 'protease: pro-IL-1b -> mature IL-1b'),
    ('IL1B', 'substrate'),
    ('GSDMD', 'pore, release of mature IL-1b'),
    ('--- other inflammasomes ---', None),
    ('NLRP1', 'sensor, NLRP1 inflammasome'),
    ('NLRC4', 'sensor, NLRC4/NAIP inflammasome'),
    ('AIM2', 'sensor, AIM2 inflammasome'),
    ('MEFV', 'sensor, pyrin inflammasome'),
    ('CARD8', 'sensor, CARD8 inflammasome'),
    ('NLRP6', 'sensor, regulatory'),
    ('NLRP7', 'sensor, regulatory'),
    ('NLRP12', 'sensor, negative regulator'),
    ('--- non-canonical and alternative processing ---', None),
    ('CASP4', 'non-canonical, human CASP4'),
    ('CASP5', 'non-canonical, human CASP5'),
    ('CASP8', 'alternative processing'),
    ('GSDME', 'alternative pore'),
    ('ELANE', 'extracellular processing'),
    ('PRTN3', 'extracellular processing'),
    ('CTSG', 'extracellular processing'),
    ('--- priming arm (induces pro-IL-1b) ---', None),
    ('NFKB1', 'priming'),
    ('RELA', 'priming'),
    ('MYD88', 'priming'),
    ('TLR2', 'priming'),
    ('TLR4', 'priming'),
    ('TICAM1', 'priming'),
    ('--- IL-1 signalling ---', None),
    ('IL18', 'co-substrate of CASP1'),
    ('IL1RN', 'receptor antagonist'),
    ('IL1R1', 'receptor'),
    ('IL1R2', 'decoy receptor'),
    ('IL1RAP', 'receptor accessory protein'),
]

SCALES = [
    ('span',     'layer1_cpg_genebodies.tsv',              '',      0.8971),
    ('exon',     'layer1_cpg_exons_bygene_genomewide.tsv', '',      0.8945),
    ('TSS+-1kb', 'layer1_cpg_promoters_tss1kb.tsv',        '_prom', 0.9130),
    ('exon1',    'layer1_cpg_first_exon.tsv',              '_ex1',  0.8945),
    ('intron1',  'layer1_cpg_first_intron.tsv',            '_in1',  0.8971),
]

_log = open(f'{OUT}/d11_run.log', 'w')

def log(m=''):
    print(m, flush=True)
    _log.write(str(m) + '\n')
    _log.flush()

tab = {}
for name, fn, suf, bg in SCALES:
    d = {}
    for r in csv.DictReader(open(f'{R}/{fn}'), delimiter='\t'):
        n = r['name'].upper()
        if suf and n.endswith(suf.upper()):
            n = n[:-len(suf)]
        for s in n.split(','):
            d.setdefault(s, r)
    tab[name] = (d, bg)

out_rows = []
for scale, _, _, bg in SCALES:
    d, _ = tab[scale]
    log('')
    log('=' * 122)
    log(f'  SCALE: {scale}     matched background {bg}')
    log('=' * 122)
    log(f'  {"gene":9}{"role":40}{"NT %":>8}{"PA %":>8}{"ratio":>8}{"rel":>7}'
        f'{"NT calls":>11}{"PA calls":>11}  note')
    for gene, role in AXIS:
        if role is None:
            log(f'  {gene}')
            continue
        r = d.get(gene)
        if r is None:
            log(f'  {gene:9}{role:40}{"-":>8}{"-":>8}{"-":>8}{"-":>7}{"-":>11}{"-":>11}'
                f'  no usable coverage at this scale')
            out_rows.append((scale, gene, role, '', '', '', '', '', '', 'no coverage'))
            continue
        ntr, par = 100 * float(r['NT_rate']), 100 * float(r['PA_rate'])
        ntn, pan = int(float(r['NT_valid'])), int(float(r['PA_valid']))
        try:
            ratio = float(r['ratio_PA_NT'])
        except ValueError:
            ratio = float('nan')
        rel = ratio / bg
        note = []
        if ntr < 10:
            note.append('NT below 10 % - uninterpretable')
        if min(ntn, pan) < 100:
            note.append(f'thin ({min(ntn, pan)} calls)')
        log(f'  {gene:9}{role:40}{ntr:>8.1f}{par:>8.1f}{ratio:>8.2f}{rel:>7.2f}'
            f'{ntn:>11,}{pan:>11,}  {"; ".join(note)}')
        out_rows.append((scale, gene, role, f'{ntr:.2f}', f'{par:.2f}', f'{ratio:.4f}',
                         f'{rel:.4f}', ntn, pan, '; '.join(note)))

with open(f'{OUT}/d11_inflammasome_axis.tsv', 'w') as o:
    o.write('scale\tgene\trole\tNT_pct\tPA_pct\tratio_PA_NT\trel_to_matched_background\t'
            'NT_calls\tPA_calls\tnote\n')
    for row in out_rows:
        o.write('\t'.join(str(x) for x in row) + '\n')

log('')
log('  rel = ratio / the background matched to that scale. Below 1 = loses more methylation than the')
log('        genome does; near 1 = moves with the genome; above 1 = retains more.')
log('  !! promoter and gene-body methylation associate with expression in OPPOSITE directions.')
log('     A gene body losing methylation is not evidence of higher expression. Do not pool the scales.')
log('  !! no q-value and no significance statement is available or derivable.')
log(f'  wrote {OUT}/d11_inflammasome_axis.tsv')
_log.close()
