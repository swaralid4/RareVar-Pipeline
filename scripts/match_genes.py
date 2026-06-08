import pandas as pd
import json

# ── Load filtered variants ─────────────────────────────────────────
df = pd.read_csv("data/processed/variants_annotated.csv")

# ── Load disease gene database ─────────────────────────────────────
with open("data/processed/disease_genes.json") as f:
    gene_db = json.load(f)

DISEASE_GENE_SET = set(gene_db.keys())

# ── Annotate variants ──────────────────────────────────────────────
df["IS_DISEASE_GENE"] = df["GENE"].isin(DISEASE_GENE_SET)

df["DISEASE_NAME"] = df["GENE"].map(
    lambda g: gene_db.get(g, {}).get("disease", "—")
)

# ── Report matches ─────────────────────────────────────────────────
matches = df[df["IS_DISEASE_GENE"]]

print(f"Total variants     : {len(df):,}")
print(f"Disease-gene hits  : {len(matches):,}")
print(f"Non-disease        : {len(df) - len(matches):,}")

if len(matches) > 0:
    print("\nMatched variants:")
    print(matches[["CHROM", "POS", "GENE", "DISEASE_NAME", "CLNSIG"]].to_string(index=False))
else:
    print("\n⚠️  No disease-gene matches found.")
    print("   This is expected with real 1000 Genomes data — GENE column is 'UNKNOWN'.")
    print("   Prioritization will still work on AF and QUAL scores.")

# ── Save ───────────────────────────────────────────────────────────
df.to_csv("data/processed/variants_annotated.csv", index=False)
print("\n✓ Saved → data/processed/variants_annotated.csv")
