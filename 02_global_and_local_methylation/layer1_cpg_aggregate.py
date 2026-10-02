#!/usr/bin/env python3
"""Pooled fraction sum(modified)/sum(valid) per region; regions with >= MIN_OBS valid calls in both arms.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, math

try:
    from scipy.stats import fisher_exact, chi2_contingency
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

nt_f, pa_f, out_f = sys.argv[1:4]
MIN_OBS = int(sys.argv[4]) if len(sys.argv) > 4 else 20

def load(path):
    d = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        v, m = f[-2], f[-1]
        v = 0 if v in ('.', '') else int(float(v))
        m = 0 if m in ('.', '') else int(float(m))
        key = (f[0], f[1], f[2], f[3] if len(f) > 3 else '.')
        d[key] = (v, m)
    return d

def wilson(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))

def newcombe(k1, n1, k2, n2):
    """Hybrid-score CI on p2 - p1 (Newcombe 1998, method 10)."""
    l1, u1 = wilson(k1, n1)
    l2, u2 = wilson(k2, n2)
    d = (k2 / n2 if n2 else 0) - (k1 / n1 if n1 else 0)
    lo = d - math.sqrt((k2 / n2 - l2) ** 2 + (u1 - k1 / n1) ** 2) if n1 and n2 else 0
    hi = d + math.sqrt((u2 - k2 / n2) ** 2 + (k1 / n1 - l1) ** 2) if n1 and n2 else 0
    return d, lo, hi

NT, PA = load(nt_f), load(pa_f)
rows, excluded = [], 0
for key in NT:
    if key not in PA:
        continue
    vN, mN = NT[key]
    vP, mP = PA[key]
    if vN < MIN_OBS or vP < MIN_OBS:
        excluded += 1
        continue
    rN, rP = mN / vN, mP / vP
    diff, lo, hi = newcombe(mN, vN, mP, vP)
    tab = [[mN, vN - mN], [mP, vP - mP]]
    if HAVE_SCIPY:
        if min(vN, vP) < 100:
            p = fisher_exact(tab)[1]
        else:
            try:
                p = chi2_contingency(tab)[1]
            except ValueError:
                p = 1.0
    else:
        p = float('nan')
    rows.append([*key, vN, mN, rN, vP, mP, rP, diff, lo, hi,
                 (rP / rN if rN > 0 else float('nan')), p])

rows.sort(key=lambda r: (r[-1] if r[-1] == r[-1] else 1.0))
n = len(rows)
qprev = 1.0
for i in range(n - 1, -1, -1):
    p = rows[i][-1]
    q = min(qprev, p * n / (i + 1)) if p == p else float('nan')
    qprev = q if q == q else qprev
    rows[i].append(q)

rows.sort(key=lambda r: (r[0], int(r[1])))
with open(out_f, 'w') as o:
    o.write("chrom\tstart\tend\tname\tNT_valid\tNT_mod\tNT_rate\t"
            "PA_valid\tPA_mod\tPA_rate\tdiff_PA_minus_NT\tdiff_lo95\tdiff_hi95\t"
            "ratio_PA_NT\tp\tq_BH\n")
    for r in rows:
        o.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.6g}"
                          for x in r) + "\n")

tv_n = sum(r[4] for r in rows); tm_n = sum(r[5] for r in rows)
tv_p = sum(r[7] for r in rows); tm_p = sum(r[8] for r in rows)
sig = sum(1 for r in rows if r[-1] == r[-1] and r[-1] < 0.05)
up = sum(1 for r in rows if r[-1] == r[-1] and r[-1] < 0.05 and r[10] > 0)
print(f"regions tested   : {len(rows):,}   excluded (<{MIN_OBS} obs in an arm): {excluded:,}")
if rows:
    print(f"pooled NT rate   : {100*tm_n/tv_n:.4f}%  ({tm_n:,} / {tv_n:,})")
    print(f"pooled PA rate   : {100*tm_p/tv_p:.4f}%  ({tm_p:,} / {tv_p:,})")
    print(f"pooled PA - NT   : {100*(tm_p/tv_p - tm_n/tv_n):+.4f} pp   "
          f"ratio {(tm_p/tv_p)/(tm_n/tv_n):.4f}x")
    print(f"BH q < 0.05      : {sig:,} regions  ({up:,} up in PA, {sig-up:,} down)")
