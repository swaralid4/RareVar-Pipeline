import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use("dark_background")

ACCENT = "#00d9ff"
FIGS   = "figures/"

# ── Load data ──────────────────────────────────────────────────────
df  = pd.read_csv("data/results/variants_prioritized.csv")
raw = pd.read_csv("data/processed/variants_raw.csv")

# ── Figure 1: Variant counts per pipeline stage ────────────────────
fig, ax = plt.subplots(figsize=(9, 4))

stages = ["Raw", "Filtered", "Disease Match", "Top Candidates"]
counts = [
    len(raw),
    len(pd.read_csv("data/processed/variants_filtered.csv")),
    int(df["IS_DISEASE_GENE"].sum()),
    int((df["PRIORITY_SCORE"] >= 8).sum()),
]
colors = ["#607b96", "#7b61ff", "#ffaa00", "#00ff9d"]

bars = ax.bar(stages, counts, color=colors, edgecolor="none", width=0.55)

for bar, count in zip(bars, counts):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 30,
        f"{count:,}",
        ha="center", va="bottom",
        fontsize=11, color="white"
    )

ax.set_title("Variant Count at Each Pipeline Stage", fontsize=14, pad=14)
ax.set_ylabel("Number of Variants")
ax.set_facecolor("#0e1620")
fig.patch.set_facecolor("#080d14")
plt.tight_layout()
plt.savefig(f"{FIGS}01_pipeline_stages.png", dpi=150)
plt.close()
print("✓ Figure 1 saved")

# ── Figure 2: fix — cast CHROM to string before reindex ───────────
fig, ax = plt.subplots(figsize=(12, 4))

chrom_order  = [str(i) for i in range(1, 23)] + ["X"]
chrom_counts = (df["CHROM"].astype(str)
                .value_counts()
                .reindex(chrom_order, fill_value=0))

ax.bar(chrom_counts.index, chrom_counts.values, color=ACCENT, alpha=0.8)
ax.set_title("Filtered Variants by Chromosome", fontsize=14, pad=14)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Variant Count")
ax.set_facecolor("#0e1620")
fig.patch.set_facecolor("#080d14")
plt.tight_layout()
plt.savefig(f"{FIGS}02_chromosome_dist.png", dpi=150)
plt.close()
print("✓ Figure 2 saved")

# ── Figure 3: fix — use distinct colors per gene ──────────────────
top = (
    df[df["IS_DISEASE_GENE"]]
    .groupby("GENE")["PRIORITY_SCORE"]
    .max()
    .sort_values(ascending=True)
    .tail(10)
)

fig, ax = plt.subplots(figsize=(9, 5))

gene_colors = ["#7b61ff", "#00d9ff", "#00ff9d", "#ffaa00", "#ff4d6d",
               "#e040fb", "#40c4ff", "#ff6e40", "#69ff47", "#18ffff"]
colors_g = gene_colors[:len(top)]

bars = ax.barh(top.index, top.values, color=colors_g, edgecolor="none")

for bar, val in zip(bars, top.values):
    ax.text(bar.get_width() - 0.3, bar.get_y() + bar.get_height() / 2,
            str(val), va="center", ha="right", fontsize=12,
            color="black", fontweight="bold")

ax.set_title("Top Candidate Genes by Priority Score", fontsize=14, pad=14)
ax.set_xlabel("Max Priority Score")
ax.set_facecolor("#0e1620")
fig.patch.set_facecolor("#080d14")
plt.tight_layout()
plt.savefig(f"{FIGS}03_top_genes.png", dpi=150)
plt.close()
print("✓ Figure 3 saved — all visualizations complete!")
