#!/usr/bin/env python3
"""TSS and TES positions per gene.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, re
from collections import defaultdict

gff, genebodies, fasta, outdir = sys.argv[1:5]

genes = []
for line in open(genebodies):
    c, s, e, name, _, strand = line.rstrip('\n').split('\t')[:6]
    genes.append((c, int(s), int(e), name, strand))

with open(f"{outdir}/tss.bed", 'w') as ot, open(f"{outdir}/tes.bed", 'w') as oe:
    for c, s, e, name, strand in genes:
        tss, tes = (s, e - 1) if strand != '-' else (e - 1, s)
        ot.write(f"{c}\t{tss}\t{tss+1}\t{name}\t.\t{strand}\n")
        oe.write(f"{c}\t{tes}\t{tes+1}\t{name}\t.\t{strand}\n")
print(f"   tss.bed / tes.bed : {len(genes)} genes")

tx_gene, tx_exons = {}, defaultdict(list)
for line in open(gff):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9:
        continue
    if f[2] in ('mRNA', 'transcript'):
        tid = re.search(r'ID=([^;]+)', f[8])
        gname = re.search(r'gene=([^;]+)', f[8])
        if tid and gname:
            tx_gene[tid.group(1)] = (gname.group(1), f[0], f[6])
    elif f[2] == 'exon':
        par = re.search(r'Parent=([^;]+)', f[8])
        if par:
            tx_exons[par.group(1)].append((int(f[3]) - 1, int(f[4])))

best = {}
for tid, exons in tx_exons.items():
    if tid not in tx_gene:
        continue
    gname, chrom, strand = tx_gene[tid]
    span = max(e for _, e in exons) - min(s for s, _ in exons)
    key = (len(exons), span)
    if gname not in best or key > best[gname][0]:
        best[gname] = (key, tid, chrom, strand, sorted(exons))

n_fe = n_fi = 0
with open(f"{outdir}/first_exon.bed", 'w') as ofe, open(f"{outdir}/first_intron.bed", 'w') as ofi:
    for gname, (_, tid, chrom, strand, exons) in best.items():
        if strand == '-':
            exons = exons[::-1]
        s, e = exons[0]
        ofe.write(f"{chrom}\t{s}\t{e}\t{gname}_ex1\t.\t{strand}\n"); n_fe += 1
        if len(exons) > 1:
            s2, e2 = exons[1]
            a, b = (e, s2) if strand != '-' else (e2, s)
            if b > a:
                ofi.write(f"{chrom}\t{a}\t{b}\t{gname}_in1\t.\t{strand}\n"); n_fi += 1
print(f"   first_exon.bed    : {n_fe}")
print(f"   first_intron.bed  : {n_fi}")

def load_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name in wanted:
                seqs[name] = ''.join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name in wanted:
        seqs[name] = ''.join(buf).upper()
    return seqs

prom = []
for line in open(f"{outdir}/promoters_tss1kb.bed"):
    c, s, e, name, _, strand = line.rstrip('\n').split('\t')[:6]
    prom.append((c, int(s), int(e), name, strand))
seqs = load_fasta(fasta, {p[0] for p in prom})

def classify(seq):
    """sliding 500 bp, step 50; return HCP / ICP / LCP"""
    hi = lo = False
    for i in range(0, max(1, len(seq) - 500 + 1), 50):
        w = seq[i:i+500]
        if len(w) < 500:
            break
        nc, ng = w.count('C'), w.count('G')
        ncg = w.count('CG')
        if nc == 0 or ng == 0:
            continue
        oe = ncg * len(w) / (nc * ng)
        gc = (nc + ng) / len(w)
        if oe >= 0.75 and gc >= 0.55:
            hi = True
        if oe >= 0.48:
            lo = True
    return 'HCP' if hi else ('ICP' if lo else 'LCP')

counts = defaultdict(int)
with open(f"{outdir}/promoters_by_cpg_class.bed", 'w') as o:
    for c, s, e, name, strand in prom:
        seq = seqs.get(c, '')[s:e]
        cls = classify(seq) if seq else 'NA'
        counts[cls] += 1
        o.write(f"{c}\t{s}\t{e}\t{name}\t{cls}\t{strand}\n")
print("   promoters_by_cpg_class.bed : " +
      "  ".join(f"{k} {v}" for k, v in sorted(counts.items())))
