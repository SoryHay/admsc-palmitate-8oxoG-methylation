#!/usr/bin/env python3
"""Aggregate per-position input for PoreMeth2.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys
import os

in_dir = sys.argv[1]
out_dir = sys.argv[2]
K = int(sys.argv[3])
MAXGAP = int(sys.argv[4]) if len(sys.argv) > 4 else 20000

os.makedirs(out_dir, exist_ok=True)
SUF = ".modkit_adapted_sorted.entropy.file.tsv"

def load(arm):
    """chrom -> {pos: (entropy, entropy_cov, beta, beta_cov)}"""
    d = {}
    with open(os.path.join(in_dir, arm + SUF)) as f:
        for line in f:
            c, p, e, ec, b, bc = line.rstrip("\n").split("\t")
            d.setdefault(c, {})[int(p)] = (float(e), int(ec), float(b), int(bc))
    return d

print(f"k = {K}, max gap = {MAXGAP} bp", flush=True)
NT = load("NT")
PA = load("PA")
print(f"  NT {sum(len(v) for v in NT.values()):,} positions   "
      f"PA {sum(len(v) for v in PA.values()):,} positions", flush=True)

out = {"NT": [], "PA": []}
n_shared = n_bins = 0

for c in sorted(set(NT) & set(PA), key=lambda x: (len(x), x)):
    shared = sorted(set(NT[c]) & set(PA[c]))
    n_shared += len(shared)

    runs, cur = [], [shared[0]] if shared else []
    for p in shared[1:]:
        if p - cur[-1] <= MAXGAP:
            cur.append(p)
        else:
            runs.append(cur); cur = [p]
    if cur:
        runs.append(cur)

    for run in runs:
        for i in range(0, len(run) - K + 1, K):
            grp = run[i:i + K]
            n_bins += 1
            for arm, src in (("NT", NT), ("PA", PA)):
                se = sec = sb = sbc = 0.0
                for p in grp:
                    e, ec, b, bc = src[c][p]
                    se += e * ec; sec += ec
                    sb += b * bc; sbc += bc
                out[arm].append(
                    f"{c}\t{grp[0]}\t{(se/sec if sec else 0):.10f}\t{int(sec)}"
                    f"\t{(sb/sbc if sbc else 0):.10f}\t{int(sbc)}")

for arm in ("NT", "PA"):
    path = os.path.join(out_dir, arm + SUF)
    with open(path, "w") as f:
        f.write("\n".join(out[arm]) + "\n")
    print(f"  wrote {len(out[arm]):,} bins -> {path}", flush=True)

print(f"\nshared CpG positions {n_shared:,}   bins emitted {n_bins:,}   "
      f"CpG per bin {K}   CpG used {n_bins*K:,} ({100*n_bins*K/max(1,n_shared):.1f}% of shared)")

import statistics
db = [float(a.split("\t")[4]) - float(b.split("\t")[4]) for a, b in zip(out["PA"], out["NT"])]
de = [float(a.split("\t")[2]) - float(b.split("\t")[2]) for a, b in zip(out["PA"], out["NT"])]

def mad(x):
    m = statistics.median(x)
    return statistics.median([abs(v - m) for v in x])

for name, v in (("delta beta", db), ("delta entropy", de)):
    z = sum(1 for t in v if t == 0) / len(v)
    print(f"{name:>14}: zero-fraction {z:.3f}   MAD {mad(v):.6f}   sd {statistics.pstdev(v):.4f}"
          f"   {'-> MAD is 0, PoreMeth2DMR will still fail' if mad(v) == 0 else '-> MAD > 0, usable'}")
