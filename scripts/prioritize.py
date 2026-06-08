import pandas as pd

# ── Load annotated variants ────────────────────────────────────────
df = pd.read_csv("data/processed/variants_annotated.csv")

# ── Scoring function ───────────────────────────────────────────────
def compute_score(row):
    score = 0
    sig   = str(row["CLNSIG"]).lower()

    # Clinical significance
    if "pathogenic" in sig and "likely" not in sig:
        score += 5
    elif "likely_pathogenic" in sig:
        score += 3

    # Allele frequency
    af = row["AF"]
    if pd.notnull(af):
        if af < 0.001:
            score += 3
        if af < 0.0001:
            score += 2  # ultra-rare bonus

    # Disease gene match
    if row["IS_DISEASE_GENE"]:
        score += 4

    # Quality bonus
    if row["QUAL"] > 60:
        score += 1

    return score

# ── Apply scoring and sort ─────────────────────────────────────────
df["PRIORITY_SCORE"] = df.apply(compute_score, axis=1)
df = df.sort_values("PRIORITY_SCORE", ascending=False).reset_index(drop=True)

# ── Score distribution summary ─────────────────────────────────────
print("=== SCORE DISTRIBUTION ===")
print(df["PRIORITY_SCORE"].value_counts().sort_index(ascending=False).head(8))

# ── Top 10 candidates ──────────────────────────────────────────────
print("\n=== TOP 10 PRIORITIZED VARIANTS ===")
cols = ["CHROM", "POS", "GENE", "DISEASE_NAME", "AF", "PRIORITY_SCORE"]
print(df[cols].head(10).to_string(index=False))

# ── Save ───────────────────────────────────────────────────────────
df.to_csv("data/results/variants_prioritized.csv", index=False)
print(f"\n✓ Saved {len(df):,} scored variants → data/results/variants_prioritized.csv")
print(f"  Top score : {df['PRIORITY_SCORE'].max()}")
print(f"  Mean score: {df['PRIORITY_SCORE'].mean():.2f}")
