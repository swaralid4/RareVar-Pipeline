import pandas as pd
import numpy as np

# ── Load filtered variants ─────────────────────────────────────────
df = pd.read_csv("data/processed/variants_filtered.csv")
print(f"Variants to annotate: {len(df):,}")

# ── Load refGene annotation (hg19) ────────────────────────────────
cols = [
    "bin", "name", "chrom", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
]
ref = pd.read_csv("data/raw/refGene.txt", sep="\t", names=cols, low_memory=False)

# Keep only chr17, use gene symbol (name2), deduplicate overlapping transcripts
ref = ref[ref["chrom"] == "chr17"][["chrom", "txStart", "txEnd", "name2"]].copy()
ref = ref.drop_duplicates().sort_values("txStart").reset_index(drop=True)
print(f"chr17 gene records   : {len(ref):,}")

# ── Positional lookup function ─────────────────────────────────────
def assign_gene(chrom, pos, gene_df):
    """Return gene name if pos falls within any transcript, else INTERGENIC."""
    hits = gene_df[
        (gene_df["txStart"] <= pos) &
        (gene_df["txEnd"]   >= pos)
    ]
    if len(hits) > 0:
        return hits.iloc[0]["name2"]
    return "INTERGENIC"

# ── Annotate chr17 variants only (all variants here are chr17) ─────
print("Annotating positions... (this may take a few minutes)")

# Vectorised interval overlap using numpy for speed
starts = ref["txStart"].values
ends   = ref["txEnd"].values
genes  = ref["name2"].values

def fast_assign(pos_array):
    result = np.full(len(pos_array), "INTERGENIC", dtype=object)
    for i, pos in enumerate(pos_array):
        idx = np.where((starts <= pos) & (ends >= pos))[0]
        if len(idx) > 0:
            result[i] = genes[idx[0]]
        if i % 200_000 == 0:
            print(f"  {i:,} / {len(pos_array):,} processed...")
    return result

df["GENE"] = fast_assign(df["POS"].values)

# ── Summary ────────────────────────────────────────────────────────
gene_counts   = (df["GENE"] != "INTERGENIC").sum()
intergenic    = (df["GENE"] == "INTERGENIC").sum()
unique_genes  = df[df["GENE"] != "INTERGENIC"]["GENE"].nunique()

print(f"\nAnnotation complete:")
print(f"  Genic variants     : {gene_counts:,}")
print(f"  Intergenic         : {intergenic:,}")
print(f"  Unique genes hit   : {unique_genes:,}")

# ── Save ───────────────────────────────────────────────────────────
df.to_csv("data/processed/variants_annotated.csv", index=False)
print("\n✓ Saved → data/processed/variants_annotated.csv")
