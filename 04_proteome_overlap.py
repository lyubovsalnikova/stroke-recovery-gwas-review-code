"""Feasibility check: how much of our transcriptomic signal is testable against the
only longitudinal plasma-proteomics stroke dataset (PXD015006, MaxQuant proteinGroups)."""
import os

import pandas as pd

# proteinGroups.txt of PRIDE PXD015006 (MaxQuant output deposited by Nguyen et al., 2020), in PXD015006_DIR
SC = os.environ["PXD015006_DIR"]
D = os.environ["PROJECT_ROOT"]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "tables")
os.makedirs(OUT, exist_ok=True)

pg = pd.read_csv(os.path.join(SC, "proteinGroups.txt"), sep="\t", low_memory=False)
print("protein groups in file: %d" % len(pg))

for col in ("Reverse", "Potential contaminant", "Only identified by site"):
    if col in pg.columns:
        n = (pg[col].astype(str) == "+").sum()
        print("  flagged %-24s %d" % (col, n))
        pg = pg[pg[col].astype(str) != "+"]
print("after MaxQuant filtering: %d protein groups" % len(pg))

# a group can list several genes; take them all, they are all "measured"
genes = set()
multi = 0
for g in pg["Gene names"].dropna():
    parts = [x.strip() for x in str(g).split(";") if x.strip()]
    if len(parts) > 1:
        multi += 1
    genes.update(parts)
print("unique gene symbols mapped: %d  (%d groups map to >1 gene)" % (len(genes), multi))

# ---- our gene space and sets
expr = pd.read_csv(os.path.join(D, "diff_expr/data/stroke_precombat_integrated.csv"),
                   index_col=0, usecols=[0])
universe = set(expr.index)

a = pd.read_csv(os.path.join(D, "diff_expr/DEG_tables/A_stroke_vs_control.csv"))
contrastA = set(a[(a.FDR < 0.05) & (a.logFC.abs() >= 0.5)].gene)

rec = pd.read_csv(os.path.join(D, "cellcomp/tables/cellclass_recovery_vs_control.csv"))
recovery_robust = set(rec[rec.category == "robust"].gene)

sp = pd.read_csv(os.path.join(D, "cellcomp/tables/spline_temporal_classification.csv"))
spline_sig = set(sp[sp.temporal_class != "not_temporal"].gene)

trp = pd.read_csv(os.path.join(D, "trajectories/temporal_programs/temporal_response_programs.csv"))
trps = {p: set(g.gene) for p, g in trp.groupby("program")}

SETS = {"our transcriptome universe": universe,
        "contrast A (stroke vs control)": contrastA,
        "recovery-robust set": recovery_robust,
        "spline time-varying (any class)": spline_sig,
        "temporal core (TRP1-3)": set(trp.gene)}
for p in sorted(trps):
    SETS["  " + p] = trps[p]

rows = []
print("\n%-36s %7s %9s %8s" % ("gene set", "n", "in MS", "%"))
print("-" * 64)
for name, s in SETS.items():
    s = s & universe if name != "our transcriptome universe" else s
    hit = s & genes
    pct = 100 * len(hit) / len(s) if s else 0
    rows.append(dict(gene_set=name.strip(), n_genes=len(s), n_measured_in_plasma=len(hit),
                     pct_measurable=round(pct, 1),
                     measured_genes=";".join(sorted(hit)) if len(hit) <= 40 else ""))
    print("%-36s %7d %9d %7.1f%%" % (name, len(s), len(hit), pct))

pd.DataFrame(rows).to_csv(os.path.join(OUT, "T11_proteome_overlap_feasibility.csv"), index=False)
print("\nwrote T11_proteome_overlap_feasibility.csv")

print("\ncontrast-A genes measurable in plasma: %s"
      % ", ".join(sorted(contrastA & genes)) or "(none)")
print("\ntemporal-core genes measurable in plasma: %s"
      % ", ".join(sorted(set(trp.gene) & genes)) or "(none)")
