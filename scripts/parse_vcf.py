import pandas as pd
import numpy as np
from cyvcf2 import VCF
from tqdm import tqdm

# ── Paths ──────────────────────────────────────────────────────────
VCF_PATH = "data/raw/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
OUT_PATH = "data/processed/variants_raw.csv"

# ── Helper: extract a single float from INFO field ─────────────────
def get_float(variant, field):
    val = variant.INFO.get(field, None)
    if val is None:
        return None
    if isinstance(val, (tuple, list)):
        val = val[0]
    try:
        return float(val)
    except (TypeError, ValueError):
        return None

# ── Helper: compute true AF from AC and AN ────────────────────────
def compute_af(variant):
    """
    The stored AF field rounds to 2 decimal places — losing precision
    for rare variants. AC and AN are stored as exact integers.
    AC=3, AN=5096 → true AF = 0.000589 (not 0.0)
    """
    ac = variant.INFO.get("AC", None)
    an = variant.INFO.get("AN", None)

    # unpack tuples
    if isinstance(ac, (tuple, list)): ac = ac[0]
    if isinstance(an, (tuple, list)): an = an[0]

    try:
        ac = int(ac)
        an = int(an)
        if an > 0:
            return round(ac / an, 8)
    except (TypeError, ValueError, ZeroDivisionError):
        pass

    # fallback to stored AF if AC/AN unavailable
    return get_float(variant, "AF")

# ── Parse VCF ─────────────────────────────────────────────────────
print(f"Opening: {VCF_PATH}")
vcf     = VCF(VCF_PATH)
records = []

for variant in tqdm(vcf, desc="Parsing variants"):
    ac = variant.INFO.get("AC", None)
    an = variant.INFO.get("AN", None)
    if isinstance(ac, (tuple, list)): ac = ac[0]
    if isinstance(an, (tuple, list)): an = an[0]

    records.append({
        "CHROM":  str(variant.CHROM).replace("chr", ""),
        "POS":    int(variant.POS),
        "ID":     variant.ID or ".",
        "REF":    variant.REF,
        "ALT":    ",".join(variant.ALT),
        "QUAL":   float(variant.QUAL) if variant.QUAL is not None else 0.0,
        "FILTER": "PASS" if variant.FILTER is None else str(variant.FILTER),
        "AF":     compute_af(variant),   # computed from AC/AN — full precision
        "AC":     int(ac) if ac is not None else None,
        "AN":     int(an) if an is not None else None,
        "GENE":   variant.INFO.get("GENE", "UNKNOWN"),
        "CLNSIG": variant.INFO.get("CLNSIG", "."),
    })

# ── Build DataFrame ────────────────────────────────────────────────
df = pd.DataFrame(records)

# ── Sanity checks ──────────────────────────────────────────────────
print(f"\nParsed              : {len(df):,} variants")
print(f"CHROM values        : {df['CHROM'].unique()}")
print(f"AF = 0 exactly      : {(df['AF'] == 0).sum():,}")
print(f"AF > 0              : {(df['AF'] > 0).sum():,}")
print(f"AF < 0.001 (rare)   : {((df['AF'] > 0) & (df['AF'] < 0.001)).sum():,}")
print(f"AF < 0.01           : {(df['AF'] < 0.01).sum():,}")
print(f"\nAF stats:")
print(df["AF"].describe())
print(f"\nSample rows — previously zero AF, now correctly computed:")
sample = df[(df["AC"] > 0) & (df["AF"] < 0.001)].head(5)
print(sample[["CHROM","POS","REF","ALT","AC","AN","AF"]].to_string())

# ── Save ───────────────────────────────────────────────────────────
df.to_csv(OUT_PATH, index=False)
print(f"\n✓ Saved → {OUT_PATH}")
