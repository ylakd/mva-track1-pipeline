"""Central configuration.

Every threshold and constant used in the report lives here, so that the methods
text and the code cannot drift apart.
"""
import os
from pathlib import Path

REPO_ID = "SageBio/mva-hackathon-2026-data"          # gated dataset (accept the terms first)
VCF_NAME = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
TBI_NAME = VCF_NAME + ".tbi"
PHENO_NAME = "Challenge_Clinical_Phenotype_1.docx"

GTF_URL = os.environ.get(
    "MVA_GTF_URL",
    "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gtf.gz",
)
HPO_URL = os.environ.get(
    "MVA_HPO_URL",
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt",
)

PROBAND_ID = "PROBAND01"

# The patient's HPO terms are NOT stored in the code. Stage 1 reads them from the challenge phenotype
# document and writes them to the (git-ignored) working directory; stage 6 loads them from there.
HPO_TERMS: dict = {}


def load_hpo_terms() -> dict:
    """Return the proband's HPO terms {id: label}: in memory if set, else from the working directory."""
    if HPO_TERMS:
        return HPO_TERMS
    import json as _json
    f = workdir() / FILES["hpo_terms"]
    if not f.exists():
        raise FileNotFoundError("hpo_terms.json not found: run stage 1 first (it reads the phenotype document)")
    HPO_TERMS.update(_json.loads(f.read_text()))
    return HPO_TERMS

QC_MIN_DP = 10
QC_MIN_GQ = 30
CDS_PAD = 2
RARE_AF = 0.001

IMPACT_HIGH = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
               "splice_donor_variant", "start_lost", "stop_lost"}
IMPACT_MODERATE = {"missense_variant", "inframe_insertion", "inframe_deletion",
                   "protein_altering_variant"}
TRUNCATING = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}

# Damage weights (set after the candidate was known: NOT blind, see report)
W_STRONG = 1.0
W_OTHER = 0.5
W_SINGLE_ALLELE = 0.25
W_HOM_MULT = 2.0

PERM_N = 2000
PERM_SEED = 42
IC_TOL = 0.5
BENCH_N = 500
BENCH_SEED = 2026

MYVARIANT_BATCH = 500
MYVARIANT_FIELDS = ",".join([
    "clinvar.rcv.clinical_significance",
    "gnomad_genome.af.af", "gnomad_exome.af.af",
    "snpeff.ann.effect", "snpeff.ann.hgvs_p", "snpeff.ann.genename",
    "dbnsfp.cadd.phred", "dbnsfp.revel", "dbnsfp.genename",
])

EPCR_PRIMARY = 0.60
EPCR_SECONDARY = 0.05
MAX_ROWS = 10

LOCUS_FLANK = 50_000
BIN_SIZE = 5_000_000
BAF_MIN_DP = 20
BAF_MIN_SITES_PER_BIN = 200


def workdir() -> Path:
    """Working directory: $MVA_WORKDIR, else /kaggle/working, else ./work."""
    env = os.environ.get("MVA_WORKDIR")
    if env:
        p = Path(env)
    elif Path("/kaggle/working").exists():
        p = Path("/kaggle/working")
    else:
        p = Path("./work")
    p.mkdir(parents=True, exist_ok=True)
    return p


FILES = {
    "vcf": "data/" + VCF_NAME, "tbi": "data/" + TBI_NAME, "pheno": "data/" + PHENO_NAME,
    "gtf": "genes.gtf.gz", "hpo": "phenotype_to_genes.txt",
    "qc": "variants_pass.parquet", "qc_stats": "qc_stats.json",
    "coding": "variants_coding.parquet", "myvariant_raw": "myvariant_raw.json",
    "annotated": "annotated.parquet",
    "scan_A": "scan_A_homozygous.csv", "scan_B": "scan_B_compound_het.csv",
    "scan_clinvar": "scan_clinvar_pathogenic.csv", "scan_summary": "scan_summary.json",
    "ranking": "phenotype_ranking.csv", "loo": "phenotype_leave_one_out.csv",
    "perm": "permutation_results.json", "selected": "selected_candidates.json",
    "verification": "verification.json",
    "locus": "locus_review.csv", "locus_summary": "locus_summary.json",
    "nmd": "nmd_check.json", "gnomad": "gnomad_check.json",
    "baf_chrom": "baf_by_chromosome.csv", "baf_bins": "baf_bins.csv", "baf_summary": "baf_summary.json",
    "hpo_terms": "hpo_terms.json", "timings": "timings.json", "submission": "submission.csv", "provenance": "provenance.json",
}


def path(key: str) -> Path:
    p = workdir() / FILES[key]
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
