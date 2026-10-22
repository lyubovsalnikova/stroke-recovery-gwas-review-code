"""Table 3 of the review: catalogued associations with functional outcome after ischemic stroke.

Input : data/gwas_catalog_stroke_outcome_severity.csv, the GWAS Catalog association export for the trait
        "stroke outcome severity measurement" (78 associations, three studies).
Output: results/tables/table3_catalogued_outcome_associations.csv (one row per catalog study accession).
Genome-wide significance: p < 5e-8; suggestive: all other catalogued associations.
"""
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = pd.read_csv(os.path.join(ROOT, "data", "gwas_catalog_stroke_outcome_severity.csv"))
g["pValue"] = g["pValue"].astype(float)


def label(r):
    return "%s (%s, p = %.0e)" % (r.mappedGenes, r.riskAllele, r.pValue)


rows = []
for acc, d in g.groupby("accessionId", sort=True):
    d = d.sort_values("pValue")
    gw = d[d.pValue < 5e-8]
    rows.append(dict(accession=acc,
                     study="%s (PMID %d)" % (d.author.iloc[0].split()[0], d.pubmedId.iloc[0]),
                     phenotype=d.traitName.iloc[0], n_catalogued=len(d), n_genome_wide=len(gw),
                     n_suggestive=len(d) - len(gw),
                     genome_wide_loci="; ".join(label(r) for r in gw.itertuples()) if len(gw) else "none",
                     lead_association=label(d.iloc[0])))
out = os.path.join(ROOT, "results", "tables")
os.makedirs(out, exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(out, "table3_catalogued_outcome_associations.csv"), index=False)
print(pd.DataFrame(rows)[["accession", "study", "n_catalogued", "n_genome_wide"]].to_string(index=False))
