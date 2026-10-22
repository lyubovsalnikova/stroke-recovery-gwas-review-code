# -*- coding: utf-8 -*-
"""Figure 3 of the review (v9; redrawn v11): three mismatches between molecular layers and the recovery outcome.

v11 (author, 2026-10-12): Arial, minimal empty space, text readable without enlarging. Drawn at the docx insert
width (6.5 in), so point sizes are print sizes (minimum 8 pt); height 7.8 -> about 5.4 in; the measured-space overlap
footnote moved to the legend. 

(A) Phase. Roundtable phases (Bernhardt et al., 2017) as background bands; one row per data type, drawn at the
    times it is sampled. Sources: blood transcriptome = tables/phase_coverage.csv (public stroke series, Supplementary
    Table S3); plasma proteome = Nguyen et al. 2020 (60 patients, 2.0-7.6 days, 3 and 12 months); NfL = Tiedt et al.
    2018 (serial to 6 months); GFAP = Barba et al. 2026 (day 1 to days 5-7); recovery-GWAS outcomes = GENISIS 24 h
    (Ibanez et al., 2022), GISCOME 60-190 days (Soderholm et al., 2019), STRONG to 1 year (Cramer et al., 2024),
    VISP to 2 years (Aldridge et al., 2023, 2024). Imaging and scales are drawn at their typical assessment times.
(B) Compartment x construct. Each data type in the cell it measures; the outcome cell outlined.

Output: figures/multiomic_mismatch.{png,pdf}
"""
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['mathtext.default'] = 'regular'

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FS = 8

H, D = 1.0, 24.0
MO, Y = 24 * 30.44, 24 * 365.25
PHASES = [("hyperacute", 0.5, D), ("acute", D, 7 * D), ("early\nsubacute", 7 * D, 3 * MO),
          ("late\nsubacute", 3 * MO, 6 * MO), ("chronic", 6 * MO, 10 * Y)]
PHASE_COL = ['#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1', '#dd3497']
MOL, CLIN, GEN = '#2166ac', '#4d4d4d', '#1b7837'

# blood-transcriptome sampling range from the audit table (acute series vs chronic comparators)
cov = pd.read_csv(os.path.join(ROOT, "data", "phase_coverage.csv"))
acute = cov[cov.role != "chronic comparator"]
chronic = cov[cov.role == "chronic comparator"]
bt_lo, bt_hi = acute.min_h.min(), acute.max_h.max()            # 3 h - 171 h
ch_lo, ch_hi = chronic.min_h.min(), chronic.max_h.max()        # 1 - 8 years
assert bt_hi <= 7.5 * D and ch_lo >= 6 * MO, (bt_hi, ch_lo)

fig = plt.figure(figsize=(6.5, 5.2))
axA = fig.add_axes([0.215, 0.565, 0.775, 0.38])
axB = fig.add_axes([0.0, 0.0, 1.0, 0.425])

# ------------------------------------------------------------------------------------------------ panel A
for (name, lo, hi), c in zip(PHASES, PHASE_COL):
    axA.axvspan(lo, hi, color=c, alpha=0.45, lw=0)
    axA.text((lo * hi) ** 0.5, 7.5, name, ha='center', va='bottom', fontsize=FS, linespacing=0.95)

rows = ["Germline genotype", "Blood transcriptome", "Plasma proteome", "NfL, GFAP (plasma)",
        "Brain imaging", "Neurological scales", "Recovery GWAS\noutcome"]
yy = {r: 7 - i for i, r in enumerate(rows)}
G = rows[-1]


def bar(row, lo, hi, color, hatch=None, alpha=1.0, h=0.42):
    y = yy[row]
    axA.add_patch(Rectangle((lo, y - h / 2), hi - lo, h, facecolor=color if not hatch else 'white',
                            edgecolor=color, hatch=hatch, alpha=alpha, lw=0.8))


def pts(row, xs, color, marker='o', s=28, dy=0.0):
    axA.scatter(xs, [yy[row] + dy] * len(xs), color=color, marker=marker, s=s, zorder=3, edgecolor='white', lw=0.5)


def note(row, x, txt, color='black', ha='left', dy=0.0):
    axA.text(x, yy[row] + dy, txt, fontsize=FS, color=color, ha=ha, va='center')


bar("Germline genotype", 0.5, 10 * Y, GEN, alpha=0.18)
note("Germline genotype", 0.6, "fixed at conception: no sampling-time problem", GEN)

bar("Blood transcriptome", bt_lo, bt_hi, MOL)
bar("Blood transcriptome", ch_lo, ch_hi, MOL, hatch='////')
note("Blood transcriptome", 0.6, "3 h\u20137 d", MOL)
axA.annotate("", xy=(6 * MO, yy["Blood transcriptome"]), xytext=(7.6 * D, yy["Blood transcriptome"]),
             arrowprops=dict(arrowstyle='<->', color='black', lw=0.7))
axA.text((7.6 * D * 6 * MO) ** 0.5, yy["Blood transcriptome"] + 0.08, "no sample", fontsize=FS, ha='center',
         va='bottom')
note("Blood transcriptome", ch_hi, "1\u20138 y, cross-sectional", MOL, ha='right', dy=0.42)

bar("Plasma proteome", 2.03 * D, 7.58 * D, MOL)
pts("Plasma proteome", [3 * MO, 12 * MO], MOL)
note("Plasma proteome", 0.6, "one longitudinal series", MOL)

bar("NfL, GFAP (plasma)", 12, 6 * MO, MOL, alpha=0.85)
note("NfL, GFAP (plasma)", 0.6, "serial", MOL)

bar("Brain imaging", 1.5, 7 * D, CLIN, alpha=0.55)
bar("Brain imaging", 1 * MO, 1 * Y, CLIN, hatch='\\\\\\\\')
note("Brain imaging", 0.6, "acute", CLIN)
note("Brain imaging", 1.15 * Y, "follow-up MRI", CLIN)

pts("Neurological scales", [2, D, 7 * D], CLIN, marker='s', s=24)
pts("Neurological scales", [90 * D, Y], 'black', marker='D', s=26)
note("Neurological scales", 0.6, "NIHSS", CLIN)
note("Neurological scales", 1.3 * Y, "mRS 90 d, 1 y", 'black')

pts(G, [D], 'black', marker='D', s=26)
bar(G, 60 * D, 190 * D, 'black', alpha=0.85, h=0.26)
axA.plot([1 * MO, 2 * Y], [yy[G] - 0.32] * 2, color='black', lw=1.0, ls=':')
axA.scatter([2 * Y], [yy[G] - 0.32], color='black', marker='>', s=14, zorder=3)
axA.scatter([Y], [yy[G] - 0.32], color='black', marker='D', s=22, zorder=3, edgecolor='white', lw=0.5)
note(G, 0.6, "GENISIS (24 h)", 'black')
note(G, 2.2 * D, "GISCOME (60\u2013190 d)", 'black', dy=0.05)
note(G, 2.4 * Y, "STRONG,\nVISP", 'black', dy=-0.12)

axA.set_xscale('log')
axA.set_xlim(0.5, 10 * Y)
axA.set_ylim(0.5, 7.45)
axA.set_yticks([yy[r] for r in rows])
axA.set_yticklabels(rows, fontsize=FS)
for t, r in zip(axA.get_yticklabels(), rows):
    t.set_color(GEN if r.startswith("Germline") else MOL if r in rows[1:4] else CLIN if r != G else 'black')
ticks = [(1, "1 h"), (D, "24 h"), (7 * D, "7 d"), (3 * MO, "3 mo"), (6 * MO, "6 mo"), (Y, "1 y"), (5 * Y, "5 y")]
axA.set_xticks([t for t, _ in ticks])
axA.set_xticklabels([l for _, l in ticks], fontsize=FS)
axA.minorticks_off()
axA.tick_params(length=2, pad=1)
axA.set_xlabel("Time after stroke onset (log scale)", fontsize=FS, labelpad=1)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)
fig.text(0.005, 0.998, "A", fontsize=12, fontweight='bold', va='top')
fig.text(0.04, 0.993, "Phase: when each layer is sampled", fontsize=9, fontweight='bold', va='top')

# ------------------------------------------------------------------------------------------------ panel B
fig.text(0.005, 0.462, "B", fontsize=12, fontweight='bold', va='top')
fig.text(0.04, 0.457, "Compartment and construct: what each layer measures", fontsize=9, fontweight='bold', va='top')
axB.set_xlim(0, 6.58)
axB.set_ylim(-0.42, 2.66)
axB.axis('off')
cols = ["Brain tissue", "CSF", "Brain-derived\nplasma proteins", "Plasma secreted\nproteins", "Blood\nleukocytes"]
constructs = ["Molecular\nabundance", "Cellular state", "Tissue structure", "Functional\noutcome"]
x0, cw, ch = 1.2, 1.06, 0.66
cells = {
    (0, 3): ("NIHSS, mRS\nhours to years", CLIN),
    (0, 2): ("MRI: lesion,\ntracts, networks", CLIN),
    (0, 1): ("histology\n(surgery, autopsy)", CLIN),
    (0, 0): ("tissue omics\n(surgery, autopsy)", MOL),
    (1, 0): ("CSF proteins\n(rarely sampled)", MOL),
    (2, 0): ("NfL, GFAP: injury\nextent, to 6 mo", MOL),
    (3, 0): ("plasma proteome\n3 d to 12 mo", MOL),
    (4, 1): ("cell counts,\ndeconvolution", MOL),
    (4, 0): ("transcriptome\n3 h to 7 d", MOL),
}
for i, c in enumerate(cols):
    axB.text(x0 + i * cw + cw / 2, -0.03, c, ha='center', va='top', fontsize=FS)
for j, c in enumerate(constructs):
    axB.text(x0 - 0.06, j * ch + ch / 2, c, ha='right', va='center', fontsize=FS)
for i in range(len(cols)):
    for j in range(len(constructs)):
        txt = cells.get((i, j))
        axB.add_patch(Rectangle((x0 + i * cw + 0.02, j * ch + 0.02), cw - 0.04, ch - 0.04,
                                facecolor='#ececec' if txt is None else 'white', edgecolor='#bdbdbd', lw=0.6))
        if txt:
            axB.text(x0 + i * cw + cw / 2, j * ch + ch / 2, txt[0], ha='center', va='center', fontsize=FS,
                     color=txt[1])
axB.add_patch(Rectangle((x0 + 0.02, 3 * ch + 0.02), cw - 0.04, ch - 0.04, facecolor='none', edgecolor='black',
                        lw=1.8))
axB.text(x0 + cw + 0.1, 3.5 * ch, "\u2190 where recovery is measured", ha='left', va='center', fontsize=FS,
         fontweight='bold')
axB.add_patch(FancyArrowPatch((x0 + 4.25 * cw, 0.95 * ch), (x0 + 1.12 * cw, 2.85 * ch), arrowstyle='-|>',
                              mutation_scale=12, color='#b2182b', lw=1.3, alpha=0.8))
axB.text(x0 + 3.5 * cw, 2.45 * ch, "distance from the outcome:\nblood \u2192 brain, molecule \u2192 function",
         fontsize=FS, color='#b2182b', ha='center', va='center',
         bbox=dict(fc='#ececec', ec='none', pad=1))

os.makedirs(os.path.join(ROOT, "results", "figures"), exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(ROOT, "results", "figures", "multiomic_mismatch." + ext), dpi=300)
print("wrote figures/multiomic_mismatch.{png,pdf}; blood transcriptome %.0f-%.0f h, chronic %.1f-%.1f y"
      % (bt_lo, bt_hi, ch_lo / Y, ch_hi / Y))
