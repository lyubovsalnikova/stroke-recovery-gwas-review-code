"""Can GSE163256 serve the worked examples of the methodological review?

GSE163256 is the only dataset in hand with a clinical outcome (mRS at 1 year, 21 ICH patients) and
baseline severity (NIHSS, GCS, ICH volume, a severity index), plus serial sampling of hematoma and blood
myeloid cells over 1-7 days. Three uses are tested here:

  A. calibration numbers for the collider / mediator simulation of Section 3.2.2-3.2.3: the observed
     severity-outcome correlations (Spearman, with bootstrap 95% CI) that a simulation of "adjusting for
     NIHSS" should reproduce;
  B. a second real-data instance of the change-score / regression-to-the-mean identity (Box 1) on a
     different tissue and disease: the glycolysis/HIF programme score in hematoma monocytes/macrophages,
     first versus last sampling day per patient - the Oldham identity should again reproduce the
     observed baseline-vs-change correlation;
  C. the "severity adjustment moves the estimate" illustration: association of the same score with mRS,
     unadjusted and adjusted for NIHSS and ICH volume - reported with its (wide) confidence intervals
     so the reader sees why n = 21 can calibrate a simulation but cannot carry an inference.

Outputs: tables/gse163256_worked_example.csv and a printed summary.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BASE = os.environ["PROJECT_ROOT"]
rng = np.random.default_rng(7)

meta = pd.read_csv(os.path.join(BASE, "sorted_cell/tables/harmonized_metadata.csv"))
meta = meta[~meta.is_control.astype(bool)].copy()
pt = meta.groupby("patient_id").agg(nihss=("nihss", "first"), gcs=("gcs", "first"),
                                    ich_volume=("ich_volume", "first"), severity_idx=("severity_idx", "first"),
                                    mrs_1yr=("mrs_1yr", "first"))
out = []


def boot_spearman(x, y, n=2000):
    r = stats.spearmanr(x, y).statistic
    idx = np.arange(len(x))
    bs = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), replace=True)
        if len(set(b)) > 3:
            bs.append(stats.spearmanr(x[b], y[b]).statistic)
    return r, np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)


# ---- A: severity-outcome coupling
print("A. severity-outcome coupling (n = %d patients, mRS at 1 year)" % len(pt))
for v in ("nihss", "gcs", "ich_volume", "severity_idx"):
    r, lo, hi = boot_spearman(pt[v].values.astype(float), pt.mrs_1yr.values.astype(float))
    print("   %-13s rho = %+.2f  [%+.2f, %+.2f]" % (v, r, lo, hi))
    out.append(dict(block="A_severity_outcome", quantity="spearman_%s_vs_mRS1y" % v, value=r, ci_low=lo, ci_high=hi, n=len(pt)))

# ---- glycolysis / HIF score in hematoma monocytes
# the panel is the glycolysis/HIF gene list of the companion study (its Supplementary Table S9e), in data/
panel = sorted(pd.read_csv(os.path.join(ROOT, "data", "glycolysis_HIF_panel.csv")).gene.astype(str))
print("\nglycolysis/HIF panel (S9e): %d genes: %s" % (len(panel), ", ".join(panel)))
logcpm = pd.read_csv(os.path.join(BASE, "GSE163256/logCPM_matrix.csv"), index_col=0)
hem = meta[(meta.compartment == "hematoma") & (meta.cell_type == "monos")]
X = logcpm.loc[[g for g in panel if g in logcpm.index], hem.sample_id]
Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1), axis=0)
score = Z.mean(axis=0)
hem = hem.assign(score=score.loc[hem.sample_id].values)
per_day = hem.groupby(["patient_id", "day_of_trial"]).agg(score=("score", "mean"), hrs=("hrs_onset", "mean")).reset_index()

# ---- B: Oldham identity on first vs last day
first = per_day.sort_values("hrs").groupby("patient_id").first()
last = per_day.sort_values("hrs").groupby("patient_id").last()
ok = first.index[(last.hrs - first.hrs) > 24]
x, y = first.loc[ok, "score"].values, last.loc[ok, "score"].values
d = y - x
r_obs = stats.pearsonr(x, d).statistic
rxy = stats.pearsonr(x, y).statistic
sx, sy = x.std(ddof=1), y.std(ddof=1)
r_oldham = (rxy * sy - sx) / np.sqrt(sx ** 2 + sy ** 2 - 2 * rxy * sx * sy)
r_rtm_equal_var = -np.sqrt((1 - rxy) / 2)
print("\nB. change-score identity in hematoma monocytes (n = %d patients with >24 h between first and last day)" % len(ok))
print("   r(first, last) = %+.2f; sd first %.2f, sd last %.2f" % (rxy, sx, sy))
print("   observed r(baseline, change) = %+.3f; Oldham identity = %+.3f; equal-variance RTM = %+.3f" % (r_obs, r_oldham, r_rtm_equal_var))
for q, v in (("r_first_last", rxy), ("r_baseline_change_observed", r_obs), ("r_baseline_change_oldham", r_oldham),
             ("r_rtm_equal_variance", r_rtm_equal_var), ("median_hours_first", float(np.median(first.loc[ok, "hrs"]))),
             ("median_hours_last", float(np.median(last.loc[ok, "hrs"])))):
    out.append(dict(block="B_change_score_identity", quantity=q, value=v, ci_low=np.nan, ci_high=np.nan, n=len(ok)))

# ---- C: score vs mRS, unadjusted and severity-adjusted
mean_score = per_day.groupby("patient_id").score.mean()
df = pt.join(mean_score.rename("score"), how="inner").dropna()
r, lo, hi = boot_spearman(df.score.values, df.mrs_1yr.values.astype(float))
print("\nC. glycolysis/HIF score vs mRS (n = %d): unadjusted rho = %+.2f [%+.2f, %+.2f]" % (len(df), r, lo, hi))
out.append(dict(block="C_score_outcome", quantity="spearman_unadjusted", value=r, ci_low=lo, ci_high=hi, n=len(df)))
def ols(yv, Xv):
    """slope and SE of the first regressor with an intercept (plain least squares)."""
    A = np.column_stack([np.ones(len(yv)), Xv])
    beta, res, rank, _ = np.linalg.lstsq(A, yv, rcond=None)
    resid = yv - A @ beta
    s2 = (resid @ resid) / (len(yv) - A.shape[1])
    cov = s2 * np.linalg.inv(A.T @ A)
    return beta[1], np.sqrt(cov[1, 1])


yv = df.mrs_1yr.values.astype(float)
b0, se0 = ols(yv, df[["score"]].values.astype(float))
out.append(dict(block="C_score_outcome", quantity="ols_slope_unadjusted", value=b0, ci_low=b0 - 1.96 * se0,
                ci_high=b0 + 1.96 * se0, n=len(df)))
for adj in (["nihss"], ["ich_volume"], ["nihss", "ich_volume"]):
    b, se = ols(yv, df[["score"] + adj].values.astype(float))
    print("   OLS slope of mRS on score: unadjusted %+.2f (SE %.2f); adjusted for %s: %+.2f (SE %.2f)"
          % (b0, se0, "+".join(adj), b, se))
    out.append(dict(block="C_score_outcome", quantity="ols_slope_adj_" + "+".join(adj), value=b, ci_low=b - 1.96 * se,
                    ci_high=b + 1.96 * se, n=len(df)))
pd.DataFrame(out).to_csv(os.path.join(ROOT, "results", "tables", "gse163256_worked_example.csv"), index=False)
print("\nwrote tables/gse163256_worked_example.csv")
