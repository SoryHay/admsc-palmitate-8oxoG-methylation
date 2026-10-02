#!/usr/bin/env python3
"""5mC and 5hmC fractions per arm and their palmitate-minus-untreated differences (Figure 13a-b).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, collections
nt_bed, pa_bed, out = sys.argv[1:4]; label = sys.argv[4] if len(sys.argv) > 4 else 'P-1'

def load(path):
    tot = collections.Counter()
    pos = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            key = c[0] + ':' + c[1] + ':' + c[5]
            valid, nmod, code = int(c[9]), int(c[11]), c[3]
            rec = pos.get(key)
            if rec is None:
                rec = pos[key] = [valid, 0, 0]; tot['valid'] += valid; tot['positions'] += 1
            if code == 'm': rec[1] += nmod; tot['m'] += nmod
            elif code == 'h': rec[2] += nmod; tot['h'] += nmod
    return tot, pos

def frac(t): return (t['m'] / t['valid'], t['h'] / t['valid']) if t['valid'] else (float('nan'),) * 2

def verdict(dm, dh):
    if dh <= 0: return 'PASSIVE route supported (delta-5hmC <= 0)'
    if dh >= 0.5 * abs(dm): return 'ACTIVE route supported (delta-5hmC > 0 and >= 50 % of |delta-5mC|)'
    return 'NEITHER (delta-5hmC > 0 but < 50 % of |delta-5mC|); no mechanism claimed'

nt_t, nt_p = load(nt_bed); pa_t, pa_p = load(pa_bed)
shared = nt_p.keys() & pa_p.keys()
def restrict(p):
    t = collections.Counter()
    for k in shared:
        v, m, h = p[k]; t['valid'] += v; t['m'] += m; t['h'] += h; t['positions'] += 1
    return t
nt_r, pa_r = restrict(nt_p), restrict(pa_p)

L = []
L.append(f'# {label} — delta-5mC and delta-5hmC, combined model 5mCG_5hmCG@v0\n')
L.append(f'Inputs: `{nt_bed}`, `{pa_bed}` (modkit pileup --cpg --filter-threshold C:0.7). Command: `p1_step8_delta.py`.\n')
for name, a, b in (('All covered CpG positions, per arm (primary)', nt_t, pa_t),
                   ('Positions covered in both arms (check only)', nt_r, pa_r)):
    fm_nt, fh_nt = frac(a); fm_pa, fh_pa = frac(b); dm, dh = fm_pa - fm_nt, fh_pa - fh_nt
    L.append(f'\n## {name}\n')
    L.append('| arm | CpG positions | valid calls | 5mC calls | 5hmC calls | 5mC % | 5hmC % |\n|---|---|---|---|---|---|---|')
    for arm, t, fm, fh in (('NT', a, fm_nt, fh_nt), ('PA', b, fm_pa, fh_pa)):
        L.append(f"| {arm} | {t['positions']:,} | {t['valid']:,} | {t['m']:,} | {t['h']:,} | {100*fm:.3f} | {100*fh:.3f} |")
    L.append(f'\ndelta-5mC (PA - NT) = **{100*dm:+.3f} pp** ; delta-5hmC = **{100*dh:+.3f} pp** ; '
             f'delta-5hmC / |delta-5mC| = **{(dh/abs(dm) if dm else float("nan")):+.3f}** ; '
             f'ratio 5mC PA/NT = {fm_pa/fm_nt if fm_nt else float("nan"):.4f} ; ratio 5hmC PA/NT = {fh_pa/fh_nt if fh_nt else float("nan"):.4f}')
    L.append(f'\n**Verdict ({ "primary" if "primary" in name else "check"}):** {verdict(dm, dh)}')
L.append('\nNo p-value and no interval: the two arms are one library each and the read-level dispersion is unaddressed. '
         'Direction and the fraction ratio only. No 5mC value here is comparable with a Rerio value (different model).')
open(out, 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
