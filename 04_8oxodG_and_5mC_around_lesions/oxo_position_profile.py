#!/usr/bin/env python3
"""8-oxo-dG call rate by position in the read.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, glob, json, os, sys
from collections import defaultdict

END_BINS = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 200), (200, 500), (500, 10**12)]
N_DECILES = 10

def load_allowlist(path):
    d = json.load(open(path))
    thr = d["thresholds"]
    allow = set(thr) - set(d.get("exclude", [])) - set(d.get("high_fp_kmers", []))
    return allow, thr

def read_lengths_from_bam(bam_path):
    import pysam
    L = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        q = r.query_length or r.infer_read_length() or 0
        if q:
            L[r.query_name] = max(L.get(r.query_name, 0), q)
    bam.close()
    return L

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modcall", nargs="+", required=True)
    ap.add_argument("--bam", default=None, help="optional: true read lengths")
    ap.add_argument("--thresholds", default="/data/8oxo_project/refs/thresholds.json")
    ap.add_argument("--label", default="sample")
    args = ap.parse_args()

    files = []
    for x in args.modcall:
        files += sorted(glob.glob(os.path.join(x, "*.txt"))) if os.path.isdir(x) else [x]

    allow, thr = load_allowlist(args.thresholds)

    if args.bam:
        length = read_lengths_from_bam(args.bam)
        src = f"bam:{os.path.basename(args.bam)}"
    else:
        length = defaultdict(int)
        for fp in files:
            with open(fp) as fh:
                next(fh)
                for ln in fh:
                    p = ln.split("\t", 3)
                    if len(p) < 3:
                        continue
                    rid = p[0]; pos = int(p[1])
                    if pos > length[rid]:
                        length[rid] = pos
        src = "modcall max-pos proxy"

    end_cg = [0] * len(END_BINS); end_hi = [0] * len(END_BINS)
    dec_cg = [0] * N_DECILES;     dec_hi = [0] * N_DECILES
    skipped = 0
    for fp in files:
        with open(fp) as fh:
            next(fh)
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                rid, kmer = p[0], p[3]
                if kmer not in allow:
                    continue
                L = length.get(rid, 0)
                if L <= 0:
                    skipped += 1
                    continue
                pos = int(p[1])
                is_hi = float(p[2]) >= thr[kmer]

                d = min(pos, L - pos)
                for i, (lo, hi) in enumerate(END_BINS):
                    if lo <= d < hi:
                        end_cg[i] += 1
                        if is_hi:
                            end_hi[i] += 1
                        break

                k = min(int(pos / L * N_DECILES), N_DECILES - 1) if L else 0
                dec_cg[k] += 1
                if is_hi:
                    dec_hi[k] += 1

    def rate(h, c):
        return h / c * 1e6 if c else 0.0

    print(f"# label={args.label}  reads={len(length)}  length_source={src}  skipped_rows={skipped}")
    print(f"# thresholds={args.thresholds}")
    print("\n[A] distance from NEAREST read end")
    print("bin_bp\tcallable_G\thi_conf\trate_per_M")
    for (lo, hi), c, h in zip(END_BINS, end_cg, end_hi):
        nm = f"{lo}-{'inf' if hi > 10**11 else hi}"
        print(f"{nm}\t{c}\t{h}\t{rate(h, c):.1f}")
    tot_c, tot_h = sum(end_cg), sum(end_hi)
    print(f"ALL\t{tot_c}\t{tot_h}\t{rate(tot_h, tot_c):.1f}")

    o_c = sum(end_cg[:4]); o_h = sum(end_hi[:4])
    i_c = sum(end_cg[4:]); i_h = sum(end_hi[4:])
    print(f"\n[B] outer 100 bp vs interior")
    print("zone\tcallable_G\thi_conf\trate_per_M")
    print(f"outer<100\t{o_c}\t{o_h}\t{rate(o_h, o_c):.1f}")
    print(f"interior>=100\t{i_c}\t{i_h}\t{rate(i_h, i_c):.1f}")
    if i_c and o_c and i_h:
        print(f"ratio_outer_over_interior\t{rate(o_h, o_c) / rate(i_h, i_c):.3f}")

    print("\n[C] relative position along read (deciles, 5'->3')")
    print("decile\tcallable_G\thi_conf\trate_per_M")
    for i, (c, h) in enumerate(zip(dec_cg, dec_hi)):
        print(f"{i*10}-{(i+1)*10}%\t{c}\t{h}\t{rate(h, c):.1f}")

if __name__ == "__main__":
    main()
