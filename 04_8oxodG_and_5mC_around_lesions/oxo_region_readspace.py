#!/usr/bin/env python3
"""8-oxo-dG rate over the reads of a region set, read space (Figure 8a).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, json, glob, os

def poisson_ci(k, lo=0.025, hi=0.975):
    from scipy.stats import chi2
    lo_k = 0.0 if k == 0 else chi2.ppf(lo, 2 * k) / 2.0
    hi_k = chi2.ppf(hi, 2 * (k + 1)) / 2.0
    return lo_k, hi_k

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modcall", nargs="+", required=True)
    ap.add_argument("--members", required=True, help="read-id list (primary overlaps region set)")
    ap.add_argument("--thresholds", default="/data/8oxo_project/refs/thresholds.json")
    ap.add_argument("--label", default="sample")
    args = ap.parse_args()

    d = json.load(open(args.thresholds))
    thr = d["thresholds"]
    allow = set(thr) - set(d.get("exclude", [])) - set(d.get("high_fp_kmers", []))
    members = {ln.strip() for ln in open(args.members) if ln.strip()}

    mc = []
    for x in args.modcall:
        mc += sorted(glob.glob(os.path.join(x, "*.txt"))) if os.path.isdir(x) else [x]

    callable_g = hi = 0
    seen = set()
    for fp in mc:
        with open(fp) as fh:
            next(fh)
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                rid, kmer = p[0], p[3]
                if rid not in members or kmer not in allow:
                    continue
                seen.add(rid)
                callable_g += 1
                if float(p[2]) >= thr[kmer]:
                    hi += 1

    rate = (hi / callable_g * 1e6) if callable_g else 0.0
    lo, hi_ci = poisson_ci(hi)
    print(json.dumps({
        "label": args.label, "member_reads_total": len(members),
        "reads_with_calls": len(seen), "callable_G": callable_g, "hi_conf": hi,
        "rate_per_M": round(rate, 1),
        "ci95_per_M": [round(lo / callable_g * 1e6, 1) if callable_g else 0.0,
                       round(hi_ci / callable_g * 1e6, 1) if callable_g else 0.0],
    }))

if __name__ == "__main__":
    main()
