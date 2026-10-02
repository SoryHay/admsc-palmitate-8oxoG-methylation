#!/usr/bin/env python3
"""CpG islands on T2T (strict criteria).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, argparse
import numpy as np

def read_fasta(path):
    """Yield (name, uppercase bytes-array) per record, streaming to limit memory."""
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, np.frombuffer(b"".join(chunks).upper(), dtype=np.uint8)
                name = line[1:].split()[0].decode()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, np.frombuffer(b"".join(chunks).upper(), dtype=np.uint8)

def islands_for_chrom(seq, win, step, min_len, min_gc, min_oe):
    C, G, N = ord("C"), ord("G"), ord("N")
    isC = (seq == C)
    isG = (seq == G)
    isN = (seq == N)
    isCG = np.zeros(len(seq), dtype=bool)
    if len(seq) > 1:
        isCG[:-1] = isC[:-1] & isG[1:]

    pC = np.concatenate(([0], np.cumsum(isC, dtype=np.int64)))
    pG = np.concatenate(([0], np.cumsum(isG, dtype=np.int64)))
    pCG = np.concatenate(([0], np.cumsum(isCG, dtype=np.int64)))
    pN = np.concatenate(([0], np.cumsum(isN, dtype=np.int64)))

    starts = np.arange(0, max(len(seq) - win, 0), step, dtype=np.int64)
    if starts.size == 0:
        return []
    ends = starts + win
    nC = pC[ends] - pC[starts]
    nG = pG[ends] - pG[starts]
    nCG = pCG[ends] - pCG[starts]
    nN = pN[ends] - pN[starts]

    valid = (win - nN)
    ok = valid > (win * 0.5)
    with np.errstate(divide="ignore", invalid="ignore"):
        gc = np.where(ok, (nC + nG) / np.maximum(valid, 1), 0.0)
        oe = np.where(ok & (nC > 0) & (nG > 0),
                      (nCG * valid) / np.maximum(nC * nG, 1), 0.0)
    passing = ok & (gc >= min_gc) & (oe >= min_oe)
    if not passing.any():
        return []

    ps = starts[passing]
    regions, cur_s, cur_e = [], ps[0], ps[0] + win
    for s in ps[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, s + win)
        else:
            regions.append((cur_s, cur_e)); cur_s, cur_e = s, s + win
    regions.append((cur_s, cur_e))

    out = []
    for s, e in regions:
        L = e - s
        if L < min_len:
            continue
        c = int(pC[e] - pC[s]); g = int(pG[e] - pG[s])
        cg = int(pCG[e] - pCG[s]); n = int(pN[e] - pN[s])
        v = L - n
        if v <= 0 or c == 0 or g == 0:
            continue
        gcf = (c + g) / v
        oef = (cg * v) / (c * g)
        if gcf >= min_gc and oef >= min_oe:
            out.append((s, e, L, gcf, oef, cg))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta"); ap.add_argument("out")
    ap.add_argument("--strict", action="store_true",
                    help="Takai & Jones 2002: 500bp / GC>=0.55 / OE>=0.65")
    ap.add_argument("--win", type=int, default=200)
    ap.add_argument("--step", type=int, default=10)
    a = ap.parse_args()
    min_len, min_gc, min_oe = (500, 0.55, 0.65) if a.strict else (200, 0.50, 0.60)

    n_tot = 0
    with open(a.out, "w") as o:
        for name, seq in read_fasta(a.fasta):
            res = islands_for_chrom(seq, a.win, a.step, min_len, min_gc, min_oe)
            for s, e, L, gc, oe, ncg in res:
                n_tot += 1
                o.write(f"{name}\t{s}\t{e}\tCGI_{n_tot}\t{L}\t.\t{gc:.3f}\t{oe:.3f}\t{ncg}\n")
            print(f"  {name}: {len(res)} islands", file=sys.stderr, flush=True)
    print(f"TOTAL islands: {n_tot}  (criteria: len>={min_len}, GC>={min_gc}, O/E>={min_oe})",
          file=sys.stderr)

if __name__ == "__main__":
    main()
