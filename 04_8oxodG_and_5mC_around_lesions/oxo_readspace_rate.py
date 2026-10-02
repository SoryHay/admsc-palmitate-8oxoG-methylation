#!/usr/bin/env python3
"""8-oxo-dG rate on the mitochondrial reads, read space (Figure 8a).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, json, glob, os, sys
import pysam

def poisson_ci(k, lo=0.025, hi=0.975):

    from scipy.stats import chi2
    lo_k = 0.0 if k == 0 else chi2.ppf(lo, 2 * k) / 2.0
    hi_k = chi2.ppf(hi, 2 * (k + 1)) / 2.0
    return lo_k, hi_k

def load_allowlist(path):
    d = json.load(open(path))
    thr = d["thresholds"]
    allow = set(thr) - set(d.get("exclude", [])) - set(d.get("high_fp_kmers", []))
    return allow, thr

def load_membership(bam_path, chrom):

    mem = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        if r.is_unmapped:
            mem.setdefault(r.query_name, (None, 0))
            continue
        mem[r.query_name] = (r.reference_name, r.mapping_quality)
    bam.close()
    return mem

def load_calls(modcall_paths, allow):

    calls = {}
    for fp in modcall_paths:
        with open(fp) as fh:
            next(fh)
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                kmer = p[3]
                if kmer not in allow:
                    continue
                calls.setdefault(p[0], []).append((kmer, float(p[2])))
    return calls

def tally(calls, thr, read_ids):

    callable_g = hi = 0
    for rid in read_ids:
        for kmer, score in calls.get(rid, ()):
            callable_g += 1
            if score >= thr[kmer]:
                hi += 1
    return callable_g, hi

def rate_block(callable_g, hi):
    rate = (hi / callable_g * 1e6) if callable_g else 0.0
    lo, hi_ci = poisson_ci(hi)
    return {
        "callable_G": callable_g,
        "hi_conf": hi,
        "rate_per_M": round(rate, 2),
        "ci95_per_M": [round(lo / callable_g * 1e6, 2) if callable_g else 0.0,
                       round(hi_ci / callable_g * 1e6, 2) if callable_g else 0.0],
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modcall", nargs="+", required=True, help="modcall .txt files or a dir")
    ap.add_argument("--bam", required=True, help="whole-genome T2T alignment of the same reads")
    ap.add_argument("--thresholds", default="/data/8oxo_project/refs/thresholds.json")
    ap.add_argument("--chrom", default="CP068254.1", help="chrM reference name")
    ap.add_argument("--label", default="sample")
    args = ap.parse_args()

    mc = []
    for x in args.modcall:
        if os.path.isdir(x):
            mc += sorted(glob.glob(os.path.join(x, "*.txt")))
        else:
            mc.append(x)

    allow, thr = load_allowlist(args.thresholds)
    mem = load_membership(args.bam, args.chrom)
    calls = load_calls(mc, allow)

    all_reads = set(calls)
    aligned = {r for r in all_reads if r in mem and mem[r][0] is not None}
    chrM_reads = {r for r in aligned if mem[r][0] == args.chrom}
    chrM_mapq20 = {r for r in chrM_reads if mem[r][1] >= 20}

    offtarget = {r for r in aligned if mem[r][0] != args.chrom}

    chrM_lowmapq = chrM_reads - chrM_mapq20

    out = {
        "label": args.label,
        "chrom": args.chrom,
        "reads": {
            "with_calls": len(all_reads),
            "aligned_primary": len(aligned),
            "chrM_primary": len(chrM_reads),
            "chrM_primary_mapq>=20": len(chrM_mapq20),
            "chrM_primary_mapq<20": len(chrM_lowmapq),
            "offtarget_primary": len(offtarget),
        },
        "envelope": {
            "a_all_reads_no_membership": rate_block(*tally(calls, thr, all_reads)),
            "b_chrM_primary_mapq0_RECOMMENDED": rate_block(*tally(calls, thr, chrM_reads)),
            "c_chrM_primary_mapq20_OLD": rate_block(*tally(calls, thr, chrM_mapq20)),
            "offtarget_primary_reads": rate_block(*tally(calls, thr, offtarget)),
        },
    }
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
