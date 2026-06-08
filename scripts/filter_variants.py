import pandas as pd

# ── Load raw variants ──────────────────────────────────────────────
df = pd.read_csv("data/processed/variants_raw.csv")

# ── Filter thresholds ──────────────────────────────────────────────
AF_THRESH = 0.01   # keeps rare variants (<1% population frequency)

# ── Apply filters ──────────────────────────────────────────────────
# Note: QUAL filter removed — hg38 2019 release stores QUAL=0.0
# for all variants. AF computed from AC/AN is the quality signal here.
# AF = 0 exactly means AC=0 (variant called but seen in zero individuals)
# — these are kept as ultra-rare candidates.
mask = (
    df["AF"].notnull() &          # requires frequency annotation
    (df["AF"] < AF_THRESH)         # rare variants only (AF < 1%)
)

filtered = df[mask].copy().reset_index(drop=True)

# ── Summary ────────────────────────────────────────────────────────
removed  = len(df) - len(filtered)
pct_kept = len(filtered) / len(df) * 100

print(f"Before  : {len(df):>10,} variants")
print(f"After   : {len(filtered):>10,} variants")
print(f"Removed : {removed:>10,} variants ({100 - pct_kept:.1f}% filtered out)")
print(f"\nAF breakdown in filtered set:")
print(f"  AF = 0 exactly  : {(filtered['AF'] == 0).sum():>10,}  (AC=0, ultra-rare)")
print(f"  AF 0–0.001      : {((filtered['AF'] > 0) & (filtered['AF'] < 0.001)).sum():>10,}  (very rare)")
print(f"  AF 0.001–0.01   : {((filtered['AF'] >= 0.001) & (filtered['AF'] < 0.01)).sum():>10,}  (rare)")

# ── Save ───────────────────────────────────────────────────────────
filtered.to_csv("data/processed/variants_filtered.csv", index=False)
print("\n✓ Saved → data/processed/variants_filtered.csv")
