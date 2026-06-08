import pandas as pd

# ── Load your annotated variants ───────────────────────────────────
df = pd.read_csv("data/processed/variants_annotated.csv")
print(f"Variants loaded       : {len(df):,}")

# ── Load ClinVar — filter to hg38 only ────────────────────────────
print("Loading ClinVar (this takes ~30 seconds)...")

clinvar = pd.read_csv(
    "data/raw/variant_summary.txt",
    sep="\t",
    low_memory=False
)

# Keep only GRCh38 rows
clinvar = clinvar[clinvar["Assembly"] == "GRCh38"].copy()
print(f"ClinVar hg38 rows     : {len(clinvar):,}")

# ── Keep only useful columns ───────────────────────────────────────
clinvar = clinvar[[
    "Chromosome",
    "Start",
    "ClinicalSignificance",
    "ReviewStatus",
    "PhenotypeList",
    "GeneSymbol"
]].copy()

# ── Normalise types for joining ────────────────────────────────────
# Your VCF has CHROM as int (17), ClinVar has Chromosome as string ("17")
clinvar["Chromosome"] = clinvar["Chromosome"].astype(str)
clinvar["Start"]      = clinvar["Start"].astype(int)

df["CHROM"] = df["CHROM"].astype(str)
df["POS"]   = df["POS"].astype(int)

# ── Join on chromosome + position ─────────────────────────────────
merged = df.merge(
    clinvar,
    left_on=["CHROM", "POS"],
    right_on=["Chromosome", "Start"],
    how="left"
)

# ── Fill unmatched variants ────────────────────────────────────────
merged["CLNSIG"]   = merged["ClinicalSignificance"].fillna(".")
merged["PHENOTYPE"] = merged["PhenotypeList"].fillna(".")
merged["REVIEW"]   = merged["ReviewStatus"].fillna(".")

# Drop redundant columns from ClinVar
merged = merged.drop(
    columns=["Chromosome", "Start", "ClinicalSignificance",
             "PhenotypeList", "ReviewStatus", "GeneSymbol"],
    errors="ignore"
)

# ── Summary ────────────────────────────────────────────────────────
annotated    = (merged["CLNSIG"] != ".").sum()
pathogenic   = merged["CLNSIG"].str.contains(
                   "Pathogenic", case=False, na=False).sum()
uncertain    = merged["CLNSIG"].str.contains(
                   "Uncertain", case=False, na=False).sum()

print(f"\nClinVar annotation summary:")
print(f"  Variants with ClinVar hit : {annotated:,}")
print(f"  Pathogenic / Likely path  : {pathogenic:,}")
print(f"  Uncertain significance     : {uncertain:,}")
print(f"  No ClinVar annotation      : {len(merged) - annotated:,}")

if annotated > 0:
    print(f"\nTop CLNSIG categories:")
    print(merged[merged["CLNSIG"] != "."]["CLNSIG"]
              .value_counts().head(8).to_string())

# ── Save ───────────────────────────────────────────────────────────
merged.to_csv("data/processed/variants_annotated.csv", index=False)
print(f"\n✓ Saved → data/processed/variants_annotated.csv")
