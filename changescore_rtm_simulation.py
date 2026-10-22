# -*- coding: utf-8 -*-
"""Reference simulation for the review's change-score / regression-to-the-mean section.

Three modules, each with a closed form checked by Monte Carlo (n = 2,000,000 per scenario):

  A. Change-score artefact for a genetic variant. A variant G shifts TRUE baseline severity by b0
     (SD units per allele) and has NO effect on the recovery process (b1 = 0). True severity tracks
     from baseline to follow-up with correlation rho_T; NIHSS is measured with reliability lam
     (var(true)/var(observed)). Three analyses of the same data:
        change score   dNIHSS = N0 - N1 (improvement positive)  ~ G
        ANCOVA         N1 ~ G + N0
        follow-up      N1 ~ G
     Closed forms (per allele, SD units of true severity):
        change score : +b0 (1 - rho_T)        the severity allele "predicts recovery"
        ANCOVA       : +rho_T b0 (1 - lam)    residual severity leakage through the noisy baseline
        follow-up    : +rho_T b0              the honest carried-over severity association
     With a true recovery effect b1 (T1 = rho_T T0 + b1 G + ...): add -b1, +b1, +b1 respectively.

  B. The correlation identity. For any two measurements X0, X1 of the same quantity,
        corr(X0, X1 - X0) = (rho s1 - s0) / sqrt(s0^2 + s1^2 - 2 rho s0 s1)   (Oldham 1962)
     exactly; the equal-variance no-effect value is -sqrt((1 - rho)/2). Checked on the three
     GSE296792 composite indices (data/T7_recovery_index_feasibility.csv, from derivations/01).

  C. Case-only (collider) selection. Liability L = a G + W, case if L > threshold(K); the other
     risk factors W also shift severity (T0 = d W + noise). Among cases, G becomes associated with
     severity although it has no severity effect. Reported as SD per allele relative to a and d.

Outputs: ../tables/A_changescore_bias.csv, B_identity_check.csv, C_collider_induced.csv;
         ../figures/changescore_rtm_simulation.png/.pdf
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "..", "results", "tables")
FIG = os.path.join(HERE, "..", "results", "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "mathtext.default": "regular", "pdf.fonttype": 42})
rng = np.random.default_rng(20260911)
N = 2_000_000


def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + list(X))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta[1:]


# ------------------------------------------------------------------ A. change-score artefact
rows = []
maf = 0.3
for rho_T in (0.5, 0.7, 0.85):
    for lam in (0.7, 0.8, 0.9, 1.0):
        for b0, b1 in ((0.10, 0.0), (0.10, -0.05)):
            g = rng.binomial(2, maf, N)
            G = (g - g.mean()) / g.std()
            T0 = b0 * G + rng.standard_normal(N)
            T1 = rho_T * T0 + b1 * G + np.sqrt(1 - rho_T ** 2) * rng.standard_normal(N)
            se = np.sqrt((1 - lam) / lam)              # var(true) = 1 -> var(error) = (1-lam)/lam
            N0 = T0 + se * rng.standard_normal(N)
            N1 = T1 + se * rng.standard_normal(N)
            b_change = ols(N0 - N1, [G])[0]
            b_ancova = ols(N1, [G, N0])[0]
            b_follow = ols(N1, [G])[0]
            rows.append(dict(rho_T=rho_T, reliability=lam, b0=b0, b1=b1,
                             change_score_sim=b_change, change_score_closed=b0 * (1 - rho_T) - b1,
                             ancova_sim=b_ancova, ancova_closed=rho_T * b0 * (1 - lam) + b1,
                             followup_sim=b_follow, followup_closed=rho_T * b0 + b1))
A = pd.DataFrame(rows)
A.round(4).to_csv(os.path.join(TAB, "A_changescore_bias.csv"), index=False)
err = max(abs(A.change_score_sim - A.change_score_closed).max(), abs(A.ancova_sim - A.ancova_closed).max(),
          abs(A.followup_sim - A.followup_closed).max())
print("A: max |simulation - closed form| = %.4f (Monte-Carlo noise ~ 0.001)" % err)
print(A[(A.b1 == 0)].pivot_table(index="rho_T", columns="reliability", values="ancova_sim").round(3))

# ------------------------------------------------------------------ B. identity on the manuscript's indices
T7 = os.path.join(HERE, "..", "data", "T7_recovery_index_feasibility.csv")
t7 = pd.read_csv(T7)
B = []
for _, r in t7.iterrows():
    s0, s1, rho = r.sd_between_patient_24h, r.sd_between_patient_168h, r.rho_24h_168h
    ident = (rho * s1 - s0) / np.sqrt(s0 ** 2 + s1 ** 2 - 2 * rho * s0 * s1)
    B.append(dict(index=r["index"], sd_24h=s0, sd_168h=s1, rho_24h_168h=rho,
                  observed_r_baseline_vs_change=r.r_baseline_vs_change_OVERLAPPING,
                  identity_value=round(ident, 3), equal_variance_null=round(-np.sqrt((1 - rho) / 2), 3),
                  nonoverlapping_r=r.r_baseline_vs_change_NONOVERLAPPING, nonoverlapping_p=r.p_nonoverlapping))
B = pd.DataFrame(B)
B.to_csv(os.path.join(TAB, "B_identity_check.csv"), index=False)
print("\nB:\n", B.to_string(index=False))
assert np.allclose(B.observed_r_baseline_vs_change, B.identity_value, atol=0.002), "identity must reproduce the observed r"

# ------------------------------------------------------------------ C. collider selection
rows = []
d = 0.3                       # other risk factors shift severity by 0.3 SD per SD
for K in (0.05, 0.15):
    thr = norm.ppf(1 - K)
    for a in (0.02, 0.05, 0.10, 0.20, 0.40):
        g = rng.binomial(2, maf, N)
        G = (g - g.mean()) / g.std()
        W = rng.standard_normal(N)
        L = a * G + W
        case = L > thr * np.sqrt(1 + a ** 2)
        T0 = d * W + rng.standard_normal(N)      # G has NO direct effect on severity
        Gc, T0c, Wc = G[case], T0[case], W[case]
        induced = ols(T0c, [Gc])[0]
        corrGW = np.corrcoef(Gc, Wc)[0, 1]
        rows.append(dict(lifetime_risk_K=K, liability_effect_a=a, n_cases=int(case.sum()),
                         corr_G_W_in_cases=round(corrGW, 4), induced_severity_effect_sd_per_allele=round(induced, 4),
                         induced_over_a=round(induced / a, 3)))
C = pd.DataFrame(rows)
C.to_csv(os.path.join(TAB, "C_collider_induced.csv"), index=False)
print("\nC:\n", C.to_string(index=False))

# ------------------------------------------------------------------ figure
# v14 (2026-10-12): drawn at the print width (6.5 in), text 8 pt; lines = closed forms, open points = Monte Carlo
# (panels A, B); panel C has no closed form, so its points (Monte Carlo) are joined by lines.
FS = 8
plt.rcParams.update({"font.size": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
                     "legend.fontsize": FS})
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.6))
ax = axes[0]
sub = A[(A.b1 == 0) & (A.reliability == 0.8)]
lab = {"change_score": "change score", "ancova": "adjusted follow-up", "followup": "unadjusted follow-up"}
for col, c in zip(lab, ("#c0392b", "#2b6cb0", "#7f8c8d")):
    ax.plot(sub.rho_T, sub[col + "_closed"] / 0.10, "-", color=c, lw=1.5, label=lab[col])
    ax.plot(sub.rho_T, sub[col + "_sim"] / 0.10, "o", mfc="white", mec=c, ms=4, mew=1.0)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("Tracking of true severity, $\\rho_T$")
ax.set_ylabel("Estimate / $b_0$")
ax.set_ylim(-0.05, 1.3)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.legend(frameon=False, loc="upper left", handlelength=1.2, borderaxespad=0.1, labelspacing=0.3)
ax = axes[1]
sub = A[(A.b1 == 0) & (A.rho_T == 0.7)]
ax.plot(sub.reliability, sub.ancova_closed / 0.10, "-", color="#2b6cb0", lw=1.5)
ax.plot(sub.reliability, sub.ancova_sim / 0.10, "o", mfc="white", mec="#2b6cb0", ms=4, mew=1.0)
ax.set_xlabel("NIHSS reliability, $\\lambda$")
ax.set_ylabel("Leak after adjustment / $b_0$")
ax.set_ylim(bottom=0)
ax = axes[2]
for K, c in ((0.05, "#c0392b"), (0.15, "#2b6cb0")):
    sub = C[C.lifetime_risk_K == K]
    ax.plot(sub.liability_effect_a, sub.induced_severity_effect_sd_per_allele, "o-", color=c, lw=1.5, ms=3.5,
            label="lifetime risk %d%%" % (100 * K))
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax.set_xlabel("Liability effect of G, $a$ (SD)")
ax.set_ylabel("Induced severity effect (SD)")
ax.legend(frameon=False, loc="lower left", handlelength=1.2, borderaxespad=0.1)
for ax, letter in zip(axes, "ABC"):  # Frontiers: uppercase panel letters
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2, pad=1)
    ax.text(-0.3, 1.04, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
fig.tight_layout(pad=0.3, w_pad=1.2)
fig.savefig(os.path.join(FIG, "changescore_rtm_simulation.png"), dpi=300)
fig.savefig(os.path.join(FIG, "changescore_rtm_simulation.pdf"))
print("\nwrote tables A/B/C and the figure")
