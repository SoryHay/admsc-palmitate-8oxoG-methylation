"""Export the plotted values of the imaging figures.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd
from scipy import stats

D = "<DATA_DISK>/imaging_AdMSC_PA"
OUT = "/data/imaging_nellie/source_data"
os.makedirs(OUT, exist_ok=True)
FIELDS = ("F0", "F1", "F2")
BASE_T = [0, 1, 2]
MOTION = ("vel", "acc", "divergence", "convergence", "vergere", "directionality", "node_")

field = {f: pd.read_csv(glob.glob(f"{D}/prod_deconvolved_{f}/*features_image.csv")[0]).set_index("t").sort_index()
         for f in FIELDS}
for fignum, suffix in ((1, "_mean"), (2, "_sum")):
    cols = [c for c in field["F0"].columns
            if c.endswith(suffix) and not any(m in c for m in MOTION)]
    rows = []
    for c in cols:
        norm, basepts = {}, []
        for f in FIELDS:
            v = field[f][c].astype(float); b = v.loc[BASE_T].mean()
            norm[f] = v / b * 100
            basepts += list(norm[f].loc[BASE_T])
        M = pd.DataFrame(norm)
        floor = float(np.std(basepts, ddof=1))
        for t in M.index:
            rows.append({"parameter": c, "hours_after_palmitate": int(t),
                         **{f"{f}_percent_of_baseline": round(M.loc[t, f], 3) for f in FIELDS},
                         **{f"{f}_raw": field[f].loc[t, c] for f in FIELDS},
                         "mean_of_3_fields_percent": round(M.loc[t].mean(), 3),
                         "sd_between_fields_percent": round(M.loc[t].std(), 3),
                         "technical_floor_percent": round(floor, 3)})
    pd.DataFrame(rows).to_csv(f"{OUT}/figure{fignum}_source_data.csv", index=False)
    print(f"figure{fignum}_source_data.csv: {len(cols)} parameters x 24 timepoints")

pc = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{D}/per_cell/*_per_cell.csv"))],
               ignore_index=True)
P3 = ["mito_area_fraction", "mean_organelle_area_um2", "organelle_density_per_100um2",
      "branch_length_mean", "organelle_solidity_raw", "branch_solidity_mean",
      "branch_aspect_ratio_mean"]
rows = []
for c in [c for c in P3 if c in pc.columns]:
    sub = pc[~pc.touches_border] if c in ("mito_area_fraction", "organelle_density_per_100um2") else pc
    for t, g in sub.groupby("t"):
        v = g[c].dropna()
        rows.append({"parameter": c, "hours_after_palmitate": int(t), "n_cells": len(v),
                     "median": v.median(), "q1": v.quantile(.25), "q3": v.quantile(.75),
                     "mean": v.mean(), "sd": v.std(),
                     "robust_cv_iqr_over_median": (v.quantile(.75) - v.quantile(.25)) / v.median()})
pd.DataFrame(rows).to_csv(f"{OUT}/figure3_source_data.csv", index=False)
pc.to_csv(f"{OUT}/figure3_per_cell_values.csv", index=False)
print(f"figure3_source_data.csv: {len(rows)} rows; figure3_per_cell_values.csv: {len(pc)} cell-timepoints")

st = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{D}/per_cell/*_stitched.csv"))],
               ignore_index=True)
n_t = st.t.nunique()
span = st.groupby(["field", "cell_id"]).t.nunique()
complete = set(span[span == n_t].index)
b = st[st.t.isin(BASE_T)]
interior = set(b[~b.touches_border].groupby(["field", "cell_id"]).size().index)
big = set(b.groupby(["field", "cell_id"]).cell_area_um2.mean().pipe(lambda s: s[s > 800]).index)
keep = complete & interior & big
st["key"] = list(zip(st.field, st.cell_id))
d = st[st.key.isin(keep)].copy()

P4 = ["area_weighted_size_um2", "frac_area_largest", "median_organelle_area_um2",
      "organelle_density_per_100um2", "mito_area_fraction", "cv_organelle_area"]
long, summ = [], []
for c in [c for c in P4 if c in d.columns]:
    ends = []
    for k, g in d.groupby("key"):
        base = g[g.t.isin(BASE_T)][c].mean()
        if not np.isfinite(base) or base == 0:
            continue
        for _, r in g.iterrows():
            long.append({"parameter": c, "field": k[0], "cell_id": k[1],
                         "hours_after_palmitate": int(r.t), "raw_value": r[c],
                         "own_baseline": base, "percent_of_own_baseline": r[c] / base * 100})
        ends.append(g[g.t == 23][c].iloc[0] / base * 100)
    ends = np.array(ends)
    n_dn = int((ends < 100).sum()); n = len(ends)
    summ.append({"parameter": c, "n_cells": n,
                 "median_percent_at_23h": np.median(ends),
                 "q1": np.percentile(ends, 25), "q3": np.percentile(ends, 75),
                 "n_decreased": n_dn, "n_increased": n - n_dn,
                 "sign_test_p": stats.binomtest(max(n_dn, n - n_dn), n, 0.5).pvalue})
pd.DataFrame(long).to_csv(f"{OUT}/figure4_source_data.csv", index=False)
pd.DataFrame(summ).to_csv(f"{OUT}/figure4_statistics.csv", index=False)
print(f"figure4_source_data.csv: {len(long)} rows ({len(keep)} cells); figure4_statistics.csv written")

with open(f"{OUT}/README.md", "w") as fh:
    fh.write("""# Source data

One table per figure, containing the values actually plotted.

| file | contents |
|---|---|
| `figure1_source_data.csv` | 16 per-object shape parameters, per field and timepoint: raw value, percentage of that field's 0-2 h baseline, mean and SD across the three fields, and the technical floor (SD of the nine baseline measurements) |
| `figure2_source_data.csv` | the same for the 16 summed parameters |
| `figure3_source_data.csv` | per-cell distributions summarised per timepoint: n, median, quartiles, mean, SD, robust CV |
| `figure3_per_cell_values.csv` | every individual cell-timepoint measurement behind figure 3 |
| `figure4_source_data.csv` | paired single-cell values: for each cell and timepoint, the raw value, that cell's own baseline, and the percentage of it |
| `figure4_statistics.csv` | per parameter: n cells, median and quartiles at 23 h, number of cells decreasing or increasing, and the sign-test P value |

Upstream of these, the unreduced tables are:
  field level    <ontDATA>/imaging_AdMSC_PA/prod_deconvolved_{F0,F1,F2}/*features_image.csv
  per organelle  same directories, *features_organelles.csv (77-96 MB per field)
  per cell       <ontDATA>/imaging_AdMSC_PA/per_cell/*_per_cell.csv and *_stitched.csv
""")
print("README.md written")
