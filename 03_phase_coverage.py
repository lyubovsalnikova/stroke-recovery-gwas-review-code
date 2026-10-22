# Run after 01_recovery_index_tables.py (reads its T9 output).
"""Sampling-time coverage of the public stroke blood-transcriptome series by clinical recovery phase (Stroke
Recovery and Rehabilitation Roundtable; Bernhardt et al., 2017): Figure 3A and Supplementary Table S3.

For each series (primary cohort, GSE163256 sorted cells, chronic comparators): the number of distinct sampling
times in the hyperacute, acute, early subacute, late subacute and chronic phases, and the earliest and latest
time. Inputs: results/tables/T9_longitudinal_design_audit.csv (from 01_recovery_index_tables.py), GSE163256
metadata (PROJECT_ROOT), data/comparator_sampling_lag.csv. Output: results/tables/phase_coverage.csv. The same
table is shipped as data/phase_coverage.csv, the input of figures/multiomic_mismatch_figure.py (Figure 3).
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NP = os.path.dirname(ROOT)
BASE = os.environ["PROJECT_ROOT"]
os.makedirs(os.path.join(ROOT, "results", "tables"), exist_ok=True)

H = 1.0
D = 24.0
MO = 24 * 30.44
Y = 24 * 365.25
PHASES = [("hyperacute", 0.5, 1 * D), ("acute", 1 * D, 7 * D), ("early subacute", 7 * D, 3 * MO),
          ("late subacute", 3 * MO, 6 * MO), ("chronic", 6 * MO, 3.5 * Y)]

# ---- primary cohort sampling times (design audit)
t9 = pd.read_csv(os.path.join(ROOT, "results", "tables", "T9_longitudinal_design_audit.csv"))
series = []
for _, r in t9.iterrows():
    times = [float(x) for x in str(r.distinct_times_h).split(";")]
    series.append((r.GEO_id + " (%s)" % ("ICH" if "hemorrhagic" in r.stroke_type else "ischemic"),
                   "primary cohort", times, r.n_samples))
# ---- GSE163256
m = pd.read_csv(os.path.join(BASE, "sorted_cell/tables/harmonized_metadata.csv"))
m = m[~m.is_control.astype(bool)]
series.append(("GSE163256 (ICH, sorted cells)", "cell-sorted", sorted(m.hrs_onset.dropna().unique().tolist()), len(m)))
# ---- comparators (submitted S2c)
# chronic comparators: sampling lag after the event, per sample (GEO series metadata), in data/
cm = pd.read_csv(os.path.join(ROOT, "data", "comparator_sampling_lag.csv"))
# trauma series (GSE36809, GSE70311) excluded from the review figure - author decision 2026-09-30
for g, lab in (("GSE22255", "past ischemic, PBMC"), ("GSE73378", "past aSAH")):
    d = cm[cm.GEO_id == g]
    lag = pd.to_numeric(d.lag_for_sample_collection_years, errors="coerce").dropna()
    lag = lag[lag > 0] if (lag > 0).any() else lag
    series.append(("%s (%s)" % (g, lab), "chronic comparator", sorted((lag * Y).unique().tolist()), len(d)))

# ---- phase coverage table
cov = []
for name, role, times, nsamp in series:
    tt = np.array(times)
    row = dict(series=name, role=role, n_samples=nsamp, min_h=tt.min(), max_h=tt.max())
    for ph, lo, hi in PHASES:
        row[ph] = int(((tt >= lo) & (tt < hi)).sum())
    cov.append(row)
cov = pd.DataFrame(cov)
cov.to_csv(os.path.join(ROOT, "results", "tables", "phase_coverage.csv"), index=False)
print(cov.to_string())

