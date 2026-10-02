#!/usr/bin/env python3
"""esox calls: 5-mer allowlist, per-k-mer thresholds, callable guanines, rate per million with Poisson interval.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, json, sys
def poisson_ci(k):
    try:
        from scipy.stats import chi2
        lo = chi2.ppf(0.025, 2*k)/2 if k>0 else 0.0
        hi = chi2.ppf(0.975, 2*(k+1))/2
        return lo, hi
    except Exception:
        import math
        return max(0,k-1.96*math.sqrt(k)), k+1.96*math.sqrt(k)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modcall", nargs="+", required=True)
    ap.add_argument("--thresholds", default="/data/8oxo_project/refs/thresholds.json")
    ap.add_argument("--label", default="sample")
    ap.add_argument("--pergene", help="TSV read_id<TAB>gene -> per-gene breakdown")
    a = ap.parse_args()

    th = json.load(open(a.thresholds))
    thr = th["thresholds"]
    allow = set(thr) - set(th.get("exclude", [])) - set(th.get("high_fp_kmers", []))
    gene_of = {}
    if a.pergene:
        for line in open(a.pergene):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2: gene_of[p[0]] = p[1]

    cand = callable_g = hiconf = 0
    pg = {}
    for fn in a.modcall:
        with open(fn) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4: continue
                rid, kmer = p[0], p[3]
                try: score = float(p[2])
                except ValueError: continue
                cand += 1
                if kmer in allow:
                    callable_g += 1
                    hit = score >= thr[kmer]
                    if hit: hiconf += 1
                    if a.pergene:
                        g = gene_of.get(rid, "NA"); d = pg.setdefault(g, [0, 0])
                        d[0] += 1; d[1] += int(hit)

    rate = hiconf/callable_g*1e6 if callable_g else 0.0
    lo, hi = poisson_ci(hiconf)
    print(f"[{a.label}] candidate-G={cand:,} callable-G={callable_g:,} hi-conf={hiconf} "
          f"-> {rate:.1f}/M callable-G  (95% CI {lo/callable_g*1e6:.1f}-{hi/callable_g*1e6:.1f})"
          if callable_g else f"[{a.label}] no callable-G")
    if a.pergene:
        print("gene\tcallable_G\thi_conf\trate_per_M")
        for g, (c, h) in sorted(pg.items(), key=lambda x: -x[1][0]):
            print(f"{g}\t{c}\t{h}\t{h/c*1e6:.1f}" if c else f"{g}\t0\t0\tNA")

if __name__ == "__main__":
    main()
