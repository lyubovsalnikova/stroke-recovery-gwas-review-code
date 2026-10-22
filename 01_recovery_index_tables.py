"""Expression-index time courses in GSE296792 (Supplementary Note S1; input of the simulation's identity check)
and the sampling-time audit of the public stroke series (T7_recovery_index_feasibility.csv, T9_longitudinal_design_audit.csv).

Inputs are outputs of the companion repository stroke-transcriptomics-code, read from the directory given by the
environment variable PROJECT_ROOT (the same variable that repository uses). The output is shipped in
results/tables/ so that the review can be checked without re-running this script.
"""
import os
import pandas as pd, numpy as np
from scipy import stats
from scipy.stats import t as tdist
D=os.environ["PROJECT_ROOT"]
OUT=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"results","tables")
os.makedirs(OUT,exist_ok=True)
meta=pd.read_csv(D+"/harmonization/metadata/master_metadata_stroke.csv").set_index('sample_id')
expr=pd.read_csv(D+"/diff_expr/data/stroke_precombat_integrated.csv",index_col=0)
common=[c for c in expr.columns if c in meta.index]; expr=expr[common]; meta=meta.loc[common]
A=pd.read_csv(D+"/diff_expr/DEG_tables/A_stroke_vs_control.csv")
sig=A[(A.FDR<0.05)&(A.logFC.abs()>=0.5)]
UP=[g for g in sig[sig.logFC>0].gene if g in expr.index]
DN=[g for g in sig[sig.logFC<0].gene if g in expr.index]
rec=pd.read_csv(D+"/cellcomp/tables/cellclass_recovery_vs_control.csv")
REP=[g for g in rec[(rec.category=='robust')&(rec.logFC_adj>0)].gene if g in expr.index]
# index labels carry the actual gene counts (they were once hardcoded, and went stale
# when the GSE125512 pairing repair changed the acute signature by one gene)
SETS={f'INFL_acute_up_{len(UP)}':UP,f'LYMPH_acute_down_{len(DN)}':DN,f'REPAIR_recovery_robust_up_{len(REP)}':REP}
def score(genes, idx):
    sub=expr.loc[genes, idx]
    z=sub.sub(sub.mean(axis=1),axis=0).div(sub.std(axis=1).replace(0,np.nan),axis=0)
    return z.mean(axis=0)
W=meta[meta.GEO_id=='GSE296792'].copy()
for k,g in SETS.items(): W[k]=score(g,W.index)
rows=[]
for k,genes in SETS.items():
    p=W.pivot_table(index='subject_key',columns='time_hours',values=k).dropna()
    x1,x2,x3=p[24.0],p[72.0],p[168.0]
    d13=x3-x1; d23=x3-x2
    rho13=stats.pearsonr(x1,x3)[0]
    rtm=-np.sqrt((1-rho13)/2)
    r_ov,p_ov=stats.pearsonr(x1,d13); r_no,p_no=stats.pearsonr(x1,d23)
    rc,pc=stats.pearsonr(W[k],W['consensus_cell_PC1'].astype(float))
    rn,pn=stats.pearsonr(W[k],W['MCP_Neutrophils'].astype(float))
    rows.append(dict(index=k,n_genes=len(genes),n_complete_case_patients=len(p),
        mean_24h=round(x1.mean(),3),mean_72h=round(x2.mean(),3),mean_168h=round(x3.mean(),3),
        sd_between_patient_24h=round(x1.std(),3),sd_between_patient_168h=round(x3.std(),3),
        delta_24_to_168=round(d13.mean(),3),sd_delta=round(d13.std(),3),
        p_delta_ttest=f"{stats.ttest_1samp(d13,0).pvalue:.3g}",
        rho_24h_168h=round(rho13,3),
        r_baseline_vs_change_OVERLAPPING=round(r_ov,3), p_overlapping=f"{p_ov:.3g}",
        r_expected_from_regression_to_mean=round(rtm,3),
        r_baseline_vs_change_NONOVERLAPPING=round(r_no,3), p_nonoverlapping=f"{p_no:.3g}",
        r_index_vs_consensus_cell_PC1=round(rc,3), p_comp=f"{pc:.3g}",
        r_index_vs_MCP_neutrophils=round(rn,3), p_neut=f"{pn:.3g}"))
t7=pd.DataFrame(rows); t7.to_csv(OUT+"/T7_recovery_index_feasibility.csv",index=False)
print(t7.T.to_string())
# T9 design audit of the proposed recovery-speed analysis
rows=[]
for gse,sub in meta.groupby('GEO_id'):
    pt=sub.pivot_table(index='subject_key',columns='time_hours',values='age',aggfunc='size')
    times=sorted(sub.time_hours.dropna().unique())
    n_full=int((pt.notna().sum(axis=1)==len(times)).sum()) if len(times)>0 else 0
    rows.append(dict(GEO_id=gse,n_samples=len(sub),n_subjects=sub.subject_key.nunique(),
        stroke_type=';'.join(sorted(sub.stroke_type.dropna().unique())),
        distinct_times_h=';'.join(str(t) for t in times),
        max_time_h=(max(times) if times else np.nan),
        n_subjects_all_timepoints=n_full,
        reaches_recovery_window_168h=bool(len(times)>0 and max(times)>=168),
        has_any_clinical_outcome_column=False))
t9=pd.DataFrame(rows).sort_values('GEO_id'); t9.to_csv(OUT+"/T9_longitudinal_design_audit.csv",index=False)
print(); print(t9.to_string(index=False))
for n in [23,24,27,33]:
    tc=tdist.ppf(0.975,n-2); print(f"n={n}: |r| for p<.05 = {np.sqrt(tc**2/(tc**2+n-2)):.3f}; 80% power |r| = {np.tanh(2.802/np.sqrt(n-3)):.3f}")
