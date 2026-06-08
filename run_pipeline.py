import subprocess
import sys
import logging
import pandas as pd
from pathlib import Path

# ── Logging setup ──────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s │ %(levelname)s │ %(message)s",
    datefmt="%H:%M:%S"
)
log = logging.getLogger("RareVar")

# ── Required files and directories ────────────────────────────────
REQUIRED_INPUTS = [
    "data/raw/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
    "data/raw/refGene.txt",
    "data/raw/variant_summary.txt",
]

REQUIRED_DIRS = [
    "data/raw", "data/processed", "data/results", "figures", "scripts"
]

# ── Pipeline steps ─────────────────────────────────────────────────
# annotate_genes and match_genes are NOT listed here —
# they are always run together in run_annotate_and_match()
STEPS_BEFORE = [
    ("Build disease gene list",  "scripts/build_gene_list.py"),
    ("Parse VCF",                "scripts/parse_vcf.py"),
    ("Filter variants",          "scripts/filter_variants.py"),
]

STEPS_AFTER = [
    ("Annotate ClinVar",         "scripts/annotate_clinvar.py"),
    ("Prioritize variants",      "scripts/prioritize.py"),
    ("Generate visualizations",  "scripts/visualize.py"),
]

# ── Helper: run a single script ────────────────────────────────────
def run_step(step_name, script):
    if not Path(script).exists():
        log.error(f"Script not found: {script}")
        sys.exit(1)

    log.info(f"Starting : {step_name}")
    result = subprocess.run(
        [sys.executable, script],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        log.error(f"FAILED   : {script}")
        log.error(result.stderr)
        sys.exit(1)

    for line in result.stdout.strip().splitlines():
        log.info(f"  {line}")

    log.info(f"✓ Done   : {step_name}")

# ── Chained step: annotate_genes → match_genes ────────────────────
# These two must always run back to back.
# annotate_genes writes gene names into variants_annotated.csv.
# match_genes reads that file immediately after.
# Running them separately risks match_genes reading stale UNKNOWN names.
def run_annotate_and_match():
    log.info("Starting : Annotate genes + Match disease genes (chained)")

    for script in ["scripts/annotate_genes.py", "scripts/match_genes.py"]:
        if not Path(script).exists():
            log.error(f"Script not found: {script}")
            sys.exit(1)

        result = subprocess.run(
            [sys.executable, script],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            log.error(f"FAILED: {script}")
            log.error(result.stderr)
            sys.exit(1)

        for line in result.stdout.strip().splitlines():
            log.info(f"  {line}")

    log.info("✓ Done   : Annotate genes + Match disease genes")

# ── Pre-flight checks ──────────────────────────────────────────────
def check_environment():
    log.info("Checking environment...")

    for d in REQUIRED_DIRS:
        Path(d).mkdir(parents=True, exist_ok=True)

    missing = [f for f in REQUIRED_INPUTS if not Path(f).exists()]
    if missing:
        log.error("Missing required input files:")
        for f in missing:
            log.error(f"  ✗ {f}")
        log.error("See README.md → Quick Start for download instructions.")
        sys.exit(1)

    log.info("✓ Environment OK — all required files present")

# ── Final summary ──────────────────────────────────────────────────
def print_summary():
    try:
        raw      = pd.read_csv("data/processed/variants_raw.csv")
        filtered = pd.read_csv("data/processed/variants_filtered.csv")
        results  = pd.read_csv("data/results/variants_prioritized.csv")

        disease_hits = int(results["IS_DISEASE_GENE"].sum())
        clinsig_hits = (results["CLNSIG"] != ".").sum()
        pathogenic   = results["CLNSIG"].str.contains(
                           "Pathogenic", case=False, na=False).sum()
        top_score    = results["PRIORITY_SCORE"].max()
        top_hits     = int((results["PRIORITY_SCORE"] >= 8).sum())

        log.info("=" * 58)
        log.info("  PIPELINE SUMMARY  ·  hg38 / GRCh38")
        log.info("=" * 58)
        log.info(f"  Raw variants parsed         : {len(raw):>12,}")
        log.info(f"  After AF filtering          : {len(filtered):>12,}")
        log.info(f"  Disease gene matches        : {disease_hits:>12,}")
        log.info(f"  ClinVar annotations         : {clinsig_hits:>12,}")
        log.info(f"  Pathogenic / likely path    : {pathogenic:>12,}")
        log.info(f"  Top candidates (score ≥ 8)  : {top_hits:>12,}")
        log.info(f"  Max priority score          : {top_score:>12}")
        log.info("=" * 58)
        log.info("  Results  → data/results/variants_prioritized.csv")
        log.info("  Figures  → figures/")
        log.info("=" * 58)
        log.info("🎉 RareVar-Pipeline complete!")

    except Exception as e:
        log.warning(f"Could not generate summary: {e}")

# ── Entry point ────────────────────────────────────────────────────
if __name__ == "__main__":
    log.info("=" * 58)
    log.info("  RareVar-Pipeline  ·  hg38 build  ·  starting")
    log.info("=" * 58)

    # 1. Check all files and directories exist
    check_environment()

    # 2. Run pre-annotation steps
    for step_name, script in STEPS_BEFORE:
        run_step(step_name, script)

    # 3. Annotate gene positions and match disease genes — always chained
    run_annotate_and_match()

    # 4. Run post-annotation steps
    for step_name, script in STEPS_AFTER:
        run_step(step_name, script)

    # 5. Print final summary
    print_summary()
