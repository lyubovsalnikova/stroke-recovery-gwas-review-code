# Minimum reporting checklist for stroke recovery GWAS

Table 6 of the review, as a file to fill in. `recovery_gwas_reporting_checklist.docx` (Word) and `recovery_gwas_reporting_checklist.csv` hold the same 24 items. For each item, state whether it is reported (Yes / No / Not applicable) and where. The checklist supplements STREGA and STROBE; it does not replace them.

| Item | Domain | Recommendation | Reported? | Section/page |
|---|---|---|---|---|
| D1 | Design | Case-only or population-based design stated; if case-only, the conditioning on stroke and the possibility of collider stratification bias acknowledged | | |
| D2 | Design | Effect of each genome-wide significant variant on stroke incidence reported from an independent incidence GWAS, or its unavailability stated as a limitation | | |
| D3 | Design | Whether baseline severity was a covariate stated and its inclusion justified; its status as a mediator or collider acknowledged | | |
| D4 | Design | At least one collider-bias correction (Slope-Hunter, Dudbridge or IPW) applied to the top hits with its assumptions and their plausibility, or the omission explained | | |
| D5 | Design | Social determinants (socioeconomic position, access to care and rehabilitation, place of residence): whether measured, and whether treated as confounder, mediator or collider, with the reason | | |
| P1 | Phenotype | Exact phenotype definition: scale, time point of assessment and the Roundtable recovery phase it falls in, and ascertainment method (prospective, record review, ICD code) | | |
| P2 | Phenotype | Dichotomization and its cut-point stated; if ordinal analysis was used, whether the proportional-odds assumption was tested and held | | |
| P3 | Phenotype | For continuous phenotypes, the correlation between baseline severity and the follow-up measure, and whether baseline was a covariate | | |
| P3a | Phenotype | For a change score, s0, s1 and rho of the two measurements, the identity value of the baseline-change correlation and its equal-variance null -sqrt((1 - rho)/2) | | |
| P3b | Phenotype | For a change score analyzed with, or adjusted for, its baseline: the reliability of the baseline instrument and the implied residual leakage rho_T b0(1 - lambda) for the top loci, or the association repeated on a non-overlapping interval | | |
| P3c | Phenotype | Whether the GWAS phenotype is a change score, a follow-up score or a baseline-adjusted follow-up score, and the sign convention of the change | | |
| P4 | Phenotype | Correlation between the primary phenotype and any secondary phenotype used in the same study | | |
| P5 | Phenotype | For cohorts treated with drugs with known pharmacogenetic interactions (e.g., clopidogrel and CYP2C19), genotype-stratified analyses, or the reason they were not performed | | |
| S1 | Dataset | Sample size per analysis (cases and controls, or effective sample size), separately by ancestry | | |
| S2 | Dataset | Ancestry composition; whether the analysis was restricted to one ancestry, or the meta-analytic method used | | |
| S3 | Dataset | Time window of outcome assessment, and whether it was fixed or variable across participants | | |
| S4 | Dataset | Availability of severity data (NIHSS or equivalent) and whether they were used | | |
| S5 | Dataset | For a proxy phenotype, its correlation with the construct of interest, or a statement that it is unknown | | |
| A1 | Analysis | Genome-wide significance threshold and any adjustment beyond 5 x 10^-8 | | |
| A2 | Analysis | Replication status of each genome-wide significant locus, with the replication sample size and ancestry; unreplicated loci stated as such | | |
| A3 | Analysis | Effect size with confidence interval, or variance explained, for each reported locus - not only the p-value | | |
| A4 | Analysis | For MR, TWAS and colocalization: the assumptions, the sensitivity analyses performed, and the conditions under which the results would not be interpretable | | |
| A5 | Analysis | For a polygenic score: the ancestry of the discovery GWAS and of the target sample, and whether the score was validated in the target ancestry | | |
| A6 | Analysis | For a multi-omic claim: the number and percentage of genes or proteins measured in both layers, and the compartment (tissue), construct, sampling time and Roundtable phase of each (Section 7) | | |

Symbols: s0, s1 = SDs of baseline and follow-up; rho = their correlation; rho_T = true baseline-follow-up tracking of severity; lambda = reliability of the baseline instrument; b0 = a variant's effect on baseline severity.
