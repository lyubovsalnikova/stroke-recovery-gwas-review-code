# -*- coding: utf-8 -*-
"""Figure 1 of the review (v11): three structural constraints in stroke recovery genetics (schematic).

Author request 2026-10-12 (Conversation(6)): draw Figure 1 in the house style:
Arial, little empty space, text readable without enlarging. Drawn at the width it is inserted at in the docx
(6.5 in), so the point sizes below are the printed sizes (minimum 8 pt).

(A) Design. G -> stroke <- other causes; restricting to cases induces a G-other-causes association (dashed), which
    reaches recovery through the other causes. Departure from S6: S6 drew the induced association as arrows from G
    and E to recovery; the induced association is between the causes of the collider (Section 2.1).
(B) Phenotype. Simulated cohort (illustration, not data): mRS dichotomized; DeltaNIHSS against baseline NIHSS with no
    recovery effect at all (true severity tracks with rho_T = 0.7, NIHSS reliability 0.8, as in Section 3.3 of the review);
    an ICD-10 I69 code as the biobank proxy.
(C) Data. N with genetic data against whether recovery is measured. Numbers: Table 2 (VISP drawn as analyzed range
    488-1,270 and 2,100 genotyped; ENIGMA = patients with MRI).

Output: figures/constraints.{png,pdf}
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['mathtext.default'] = 'regular'
FS, FT = 8, 9            # body text, panel titles (pt at print size)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BLUE, GREY, RED, GREEN = '#2166ac', '#4d4d4d', '#b2182b', '#1b7837'

fig = plt.figure(figsize=(6.5, 5.0))


def title(x, y, letter, text):
    fig.text(x, y, letter, fontsize=12, fontweight='bold', va='top')
    fig.text(x + 0.03, y - 0.004, text, fontsize=FT, fontweight='bold', va='top')


# ------------------------------------------------------------------------------------------------ A: design
title(0.005, 0.995, "A", "Design: every case-only study\nconditions on a collider")
axA = fig.add_axes([0.0, 0.39, 0.44, 0.55])
axA.set_xlim(0, 10)
axA.set_ylim(0, 9.5)
axA.axis('off')


def node(x, y, txt, fc='white', ec=GREY, lw=1.0, w=2.9, h=1.25):
    axA.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.08,rounding_size=0.25",
                                 fc=fc, ec=ec, lw=lw))
    axA.text(x, y, txt, ha='center', va='center', fontsize=FS)


def arrow(p, q, color=GREY, ls='-', both=False):
    axA.add_patch(FancyArrowPatch(p, q, arrowstyle='<|-|>' if both else '-|>', mutation_scale=10, color=color,
                                  lw=1.1, ls=ls, shrinkA=0, shrinkB=0))


node(1.6, 8.3, "Variant G", w=2.6)
node(8.1, 8.3, "Other causes\n(risk factors,\nage)", w=3.0, h=2.0)
node(4.1, 5.0, "Stroke\n(cases only)", fc='#f2f2f2', ec='black', lw=1.8)
node(4.1, 1.0, "Recovery\n(mRS, NIHSS)", w=3.1, h=1.5)
arrow((2.0, 7.6), (3.5, 5.7))
arrow((7.0, 7.25), (4.9, 5.7))
arrow((4.1, 4.3), (4.1, 1.8))
axA.plot([9.3, 9.3, 5.8], [7.25, 1.0, 1.0], color=GREY, lw=1.1)
arrow((6.0, 1.0), (5.72, 1.0))
arrow((1.0, 7.6), (2.6, 1.3), color=BLUE, ls=(0, (1, 2)))
axA.text(0.95, 4.0, "effect\nsought", fontsize=FS, color=BLUE, ha='center', va='center')
arrow((3.0, 8.3), (6.5, 8.3), color=RED, ls=(0, (4, 2)), both=True)
axA.text(4.75, 8.0, "induced\namong cases", fontsize=FS, color=RED, ha='center', va='top')
axA.text(5.85, 3.4, "Adjusting for\nbaseline NIHSS\n(measured after\nthe stroke) does\nnot remove it", fontsize=FS,
         ha='left', va='center', color=GREY)

# ------------------------------------------------------------------------------------------------ B: phenotype
title(0.465, 0.995, "B", "Phenotype: three measures of the same patients")
rng = np.random.default_rng(1)
n = 400
t0 = rng.normal(0, 1, n)                              # true baseline severity
t1 = 0.7 * t0 + np.sqrt(1 - 0.7 ** 2) * rng.normal(0, 1, n)  # follow-up; no recovery effect
err = np.sqrt(0.25)                                   # reliability 0.8: error var = 0.25 x true var
o0 = np.clip(np.round(9 + 5 * (t0 + err * rng.normal(0, 1, n))), 0, 42)   # observed NIHSS (illustrative)
o1 = np.clip(np.round(5 + 5 * (t1 + err * rng.normal(0, 1, n))), 0, 42)
delta = o0 - o1
r = np.corrcoef(o0, delta)[0, 1]

# B1 mRS
ax1 = fig.add_axes([0.61, 0.80, 0.37, 0.12])
mrs = np.array([12, 20, 18, 16, 15, 9, 10])
ax1.bar(range(7), mrs, color=[BLUE] * 3 + [RED] * 4, width=0.75)
ax1.set_xticks(range(7))
ax1.tick_params(labelsize=FS, length=2, pad=1)
ax1.set_yticks([])
for s in ('top', 'right', 'left'):
    ax1.spines[s].set_visible(False)
ax1.set_xlabel("mRS at 90 days", fontsize=FS, labelpad=1)
ax1.text(1, 21.5, "0\u20132", color=BLUE, ha='center', fontsize=FS, fontweight='bold')
ax1.text(4.5, 21.5, "3\u20136", color=RED, ha='center', fontsize=FS, fontweight='bold')
ax1.set_ylim(0, 25)
fig.text(0.475, 0.86, "Ordinal,\nanalyzed\nas two\ngroups", fontsize=FS, va='center')

# B2 DeltaNIHSS vs baseline
ax2 = fig.add_axes([0.61, 0.49, 0.37, 0.235])
ax2.scatter(o0 + rng.uniform(-0.3, 0.3, n), delta, s=4, color=GREY, alpha=0.5, lw=0)
b = np.polyfit(o0, delta, 1)
xs = np.array([o0.min(), o0.max()])
ax2.plot(xs, b[1] + b[0] * xs, color=RED, lw=1.3)
ax2.axhline(0, color='black', lw=0.5)
ax2.set_xlabel("Baseline NIHSS", fontsize=FS, labelpad=1)
ax2.set_ylabel("\u0394NIHSS", fontsize=FS, labelpad=1)
ax2.tick_params(labelsize=FS, length=2, pad=1)
for s in ('top', 'right'):
    ax2.spines[s].set_visible(False)
ax2.text(0.02, 0.97, "r = %.2f with no recovery effect\n(simulated; Section 3.3)" % r, transform=ax2.transAxes,
         fontsize=FS, va='top', color=RED, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))
fig.text(0.475, 0.59, "Coupled\nto its\nbaseline", fontsize=FS, va='center')

# B3 ICD-10 proxy
fig.text(0.475, 0.385, "Proxy", fontsize=FS, va='center')
fig.text(0.61, 0.385, "ICD-10 I69 code (sequelae of stroke):\nyes or no, dated by a health-care\ncontact; no scale",
         fontsize=FS, va='center', color=GREY)

# ------------------------------------------------------------------------------------------------ C: data
title(0.005, 0.335, "C", "Data: no published dataset combines a recovery measure with a large N")
axC = fig.add_axes([0.24, 0.085, 0.74, 0.20])
axC.set_xscale('log')
axC.set_xlim(400, 400000)
axC.set_ylim(-0.6, 1.75)
axC.add_patch(Rectangle((32000, 0.5), 368000, 1.25, fc='#f2f2f2', ec='#bdbdbd', lw=0.6, ls='--'))
axC.text(113000, 1.12, "no published\ndataset", fontsize=FS, ha='center', va='center', color=GREY, style='italic')
hi = [("STRONG", 763, 1.5), ("VISP", 1270, 0.75), ("ENIGMA", 2100, 1.15), ("GISCOME", 6165, 0.9),
      ("STROMICS", 10241, 1.45)]
lo = [("CKB", 13000, 0.33), ("BBJ", 19702, -0.33), ("UK Biobank", 21047, 0.0), ("FinnGen", 33871, 0.33),
      ("MVP + UKB", 91909, -0.33)]
axC.plot([488, 2100], [0.75, 0.75], color=BLUE, lw=1.2)
for name, x, y in hi:
    axC.scatter([x], [y], color=BLUE, s=22, zorder=3)
    axC.text(x * (1.9 if name == "VISP" else 1.12), y, name, fontsize=FS, va='center', color=BLUE)
for name, x, y in lo:
    axC.scatter([x], [y], color=GREY, s=22, zorder=3)
    axC.text(x * 1.12, y, name, fontsize=FS, va='center', color=GREY)
axC.axhline(0.5, color='#bdbdbd', lw=0.6)
axC.set_yticks([0.0, 1.1])
axC.set_yticklabels(["Proxy only\n(ICD-10, EHR, death)", "Recovery measured\n(scales, imaging)"], fontsize=FS)
axC.set_xticks([1000, 10000, 100000])
axC.set_xticklabels(["1,000", "10,000", "100,000"], fontsize=FS)
axC.minorticks_off()
axC.tick_params(length=2, pad=1)
axC.set_xlabel("Stroke patients with genetic data (log scale)", fontsize=FS, labelpad=1)
for s in ('top', 'right'):
    axC.spines[s].set_visible(False)

os.makedirs(os.path.join(ROOT, "results", "figures"), exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(ROOT, "results", "figures", "constraints." + ext), dpi=300)
print("wrote figures/constraints.{png,pdf}; simulated r(baseline, delta) = %.2f" % r)
