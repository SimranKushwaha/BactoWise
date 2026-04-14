"""
================================================================================
UNIFIED ANNOTATION PIPELINE v3
================================================================================
Multi-Tool Consensus Genome Annotation Pipeline
Combines Prokka + Bakta annotations for standard bacteria

INPUT
-----
  --input   folder containing:
              genome_bakta.gff3        (Bakta GFF3 output)
              genome_prokka.gff        (Prokka GFF output)
              genome.fasta / .fna / .fa  (assembly used for annotation)
  --output  folder for all results
  --prefix  locus tag prefix  [default: GENE]

OUTPUT (only these files are kept)
------------------------------------
  Master_Table_Annotation.xlsx   final consensus annotation table
  summary_report.txt             pipeline statistics
  pipeline.log                   full execution log
  output.faa                     protein sequences (high-confidence genes only)
  output.fna                     nucleotide CDS sequences (high-confidence only)
  output.gff                     GFF3 compatible with Geneious / SnapGene
  output.gbk                     GenBank flat file compatible with Geneious /
                                 SnapGene / Addgene

CONFIDENCE FILTER FOR SEQUENCES
---------------------------------
  Only Consensus_2/3 and Consensus_3/3 features are included in FAA / FNA.
  All features appear in GFF and GBK (with confidence noted in qualifiers).

LOCUS TAGS
-----------
  New clean locus tags generated as <PREFIX>_00001, <PREFIX>_00002 …
  assigned in genomic order (contig → strand → start coordinate).

================================================================================
"""

import os
import re
import csv
import logging
import argparse
import textwrap
from pathlib import Path
from datetime import datetime
from collections import Counter
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter


# ==============================================================================
# LOGGING
# ==============================================================================
def setup_logging(output_folder: str) -> logging.Logger:
    os.makedirs(output_folder, exist_ok=True)
    log_path = os.path.join(output_folder, "pipeline.log")
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.FileHandler(log_path, mode="w"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)

logger = logging.getLogger(__name__)


# ==============================================================================
# CONSTANTS
# ==============================================================================
VALID_FEATURE_TYPES = {
    "CDS", "rRNA", "tRNA", "tmRNA", "ncRNA",
    "CRISPR", "oriC", "oriV", "oriT", "regulatory_region",
    "pseudogene",   # NCBI PGAP writes pseudogenes as a standalone type
}

TOOL_PRIORITY = ["bakta", "prokka", "pgap"]
MAX_TOOL_SUPPORT = 5

# Codon table (standard genetic code)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
LOCAL_CONFLICT_REVIEW_MAX_CONFIDENCE = 85.0

# High confidence = any consensus level (2+ tools agreed), regardless of tool count
# Evaluated dynamically: any Consensus-Level starting with "Consensus_" qualifies
def is_high_confidence(level: str) -> bool:
    return str(level).startswith("Consensus_")

FUNCTIONAL_CATEGORIES = {
    "Replication":   ["dna", "rep", "polA", "polC", "dnaN", "dnaA", "dnaB", "dnaC"],
    "Transcription": ["rpo", "sig", "nusA", "nusB", "nusG", "greA", "greB"],
    "Translation":   ["rps", "rpl", "rpm", "tuf", "tsf", "infA", "infB", "infC", "fusA"],
    "Cell_Division": ["fts", "min", "mre", "rod", "sep", "zap"],
    "Metabolism":    ["acc", "fab", "pyr", "pur", "his", "trp", "leu", "ile", "val"],
    "Stress":        ["gro", "clp", "lon", "htpG", "ibp", "dnak", "grpE"],
    "AMR":           ["bla", "aac", "aph", "ant", "tet", "erm", "mec", "van", "sul", "dfr"],
    "Mobile":        ["tnp", "int", "ins", "tra", "mob", "rep", "ist"],
    "Transport":     ["abc", "pts", "mdr", "efflux", "permease", "transporter"],
}


# ==============================================================================
# HELPER: DETECT TOOL
# ==============================================================================
def detect_tool(filename: str) -> str:
    name = Path(filename).stem.lower()
    if "ncbi" in name:
        return "pgap"
    for tool in TOOL_PRIORITY:
        if tool in name:
            return tool
    return re.sub(r"[^a-z0-9]+", "_", name).strip("_") or "unknown"


# ==============================================================================
# HELPER: FASTA PARSER (pure Python, no BioPython)
# ==============================================================================
def parse_fasta(fasta_path: str) -> dict:
    """
    Parse a FASTA file into a dict: {seq_id: sequence_string (uppercase)}.
    Handles multi-line FASTA and both .fasta / .fna / .fa extensions.
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq).upper()
                # Take only the first word of the header as the ID
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq).upper()

    return sequences


def normalize_seq_label(seq_id: str) -> str:
    seq = str(seq_id).strip()
    seq = re.sub(r"\.(fasta|fna|fa)$", "", seq, flags=re.IGNORECASE)
    seq = re.sub(r"\.\d+$", "", seq)
    return seq


def preferred_seq_id(seq_ids: list[str]) -> str:
    def rank(seq_id: str) -> tuple[int, int, str]:
        seq = str(seq_id)
        if seq.startswith("CP"):
            return (0, len(seq), seq)
        if seq.startswith("NC_"):
            return (2, len(seq), seq)
        return (1, len(seq), seq)

    return sorted(seq_ids, key=rank)[0]


def harmonize_sequence_ids(df: pd.DataFrame, genome: dict) -> tuple[pd.DataFrame, dict]:
    """
    Align contig IDs across GFF and FASTA inputs.

    The ST474 inputs mix GenBank (`CP...`), RefSeq (`NC_...`), and FASTA-style
    headers (`*.fasta`) for the same replicons. We normalize by direct label
    match first, then by unique contig length.
    """
    df = df.copy()
    genome_lengths = {seq_id: len(seq) for seq_id, seq in genome.items()}
    genome_by_norm: dict[str, list[str]] = {}
    for seq_id in genome_lengths:
        genome_by_norm.setdefault(normalize_seq_label(seq_id), []).append(seq_id)

    seq_lengths = (
        df.groupby("Sequence_ID", dropna=False)["End"]
        .max()
        .astype(int)
        .to_dict()
    )
    length_clusters: list[dict] = []
    for seq_id, seq_len in sorted(seq_lengths.items(), key=lambda item: item[1], reverse=True):
        seq_str = str(seq_id)
        seq_len = int(seq_len)
        matched = None
        for cluster in length_clusters:
            if abs(cluster["rep_len"] - seq_len) <= 5000:
                matched = cluster
                break
        if matched is None:
            matched = {"rep_len": seq_len, "seq_ids": []}
            length_clusters.append(matched)
        matched["seq_ids"].append(seq_str)
        matched["rep_len"] = max(matched["rep_len"], seq_len)

    seq_map: dict[str, str] = {}
    canonical_lengths: list[tuple[str, int]] = []
    for cluster in length_clusters:
        canonical = preferred_seq_id(cluster["seq_ids"])
        cluster["canonical"] = canonical
        cluster["canonical_len"] = max(seq_lengths[seq_id] for seq_id in cluster["seq_ids"])
        canonical_lengths.append((canonical, int(cluster["canonical_len"])))
        for seq_id in cluster["seq_ids"]:
            seq_map[seq_id] = canonical

    genome_map: dict[str, str] = {}
    genome_target_counts: Counter = Counter()
    for genome_id, genome_len in genome_lengths.items():
        norm = normalize_seq_label(genome_id)
        candidates = []
        for cluster in length_clusters:
            for seq_id in cluster["seq_ids"]:
                if normalize_seq_label(seq_id) == norm:
                    candidates.append(seq_id)
        if len(candidates) == 1:
            target = candidates[0]
        elif canonical_lengths:
            target = min(
                canonical_lengths,
                key=lambda item: abs(item[1] - int(genome_len))
            )[0]
        else:
            target = genome_id
        genome_target_counts[target] += 1
        genome_map[genome_id] = target

    normalized_genome = {}
    for genome_id, seq in genome.items():
        target = genome_map[genome_id]
        if genome_target_counts[target] == 1:
            normalized_genome[target] = seq
        else:
            normalized_genome[genome_id] = seq

    df["Sequence_ID"] = df["Sequence_ID"].astype(str).map(lambda x: seq_map.get(x, x))
    return df, normalized_genome


# ==============================================================================
# HELPER: SEQUENCE OPERATIONS
# ==============================================================================
def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def extract_subsequence(genome: dict, seq_id: str, start: int, end: int,
                         strand: str) -> str:
    """
    Extract nucleotide sequence from genome dict.
    Coordinates are 1-based inclusive (GFF convention).
    Returns reverse complement for minus strand.
    """
    contig = genome.get(seq_id, "")
    if not contig:
        return ""
    # Convert to 0-based Python slicing
    subseq = contig[start - 1: end]
    if strand == "-":
        subseq = reverse_complement(subseq)
    return subseq


def translate(nucleotide: str) -> str:
    """
    Translate a nucleotide sequence to protein using the standard genetic code.
    - Forces first codon to Met (M) regardless of actual codon, because bacteria
      commonly use alternative start codons (GTG, TTG, CTG) which are valid but
      would otherwise be translated as Val/Leu instead of Met.
    - Stops at first in-frame stop codon.
    - Returns empty string if sequence is shorter than 3 nt.
    """
    if len(nucleotide) < 3:
        return ""
    protein = []
    for i in range(0, len(nucleotide) - 2, 3):
        codon = nucleotide[i:i + 3].upper()
        if i == 0:
            aa = "M"   # force Met at start — handles GTG/TTG/CTG start codons
        else:
            aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


# ==============================================================================
# HELPER: FIND INPUT FASTA
# ==============================================================================
def find_fasta(input_folder: str) -> str:
    """
    Locate the FASTA file in the input folder.
    Accepts .fasta, .fna, .fa extensions.
    Raises FileNotFoundError if none found.
    """
    for ext in (".fasta", ".fna", ".fa"):
        for f in os.listdir(input_folder):
            if f.lower().endswith(ext):
                return os.path.join(input_folder, f)
    raise FileNotFoundError(
        f"No FASTA file (.fasta / .fna / .fa) found in: {input_folder}"
    )


# ==============================================================================
# HELPER: ASSIGN LOCUS TAGS
# ==============================================================================
def assign_locus_tags(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """
    Assign clean locus tags (<PREFIX>_00001 …) in genomic order:
    sorted by Sequence_ID → Start → Direction.
    Tags are assigned to ALL features (CDS, rRNA, tRNA …).
    """
    df = df.copy()
    df_sorted = df.sort_values(
        ["Sequence_ID", "Start", "Direction"]
    ).reset_index(drop=True)

    tags = []
    for i in range(len(df_sorted)):
        locus_num = 1 if i == 0 else i * 5
        tags.append(f"{prefix}_{locus_num:05d}")
    df_sorted["Locus_Tag"] = tags

    # Merge tags back to original order
    df = df.merge(
        df_sorted[["Sequence_ID", "Start", "End", "Direction", "Type", "Locus_Tag"]],
        on=["Sequence_ID", "Start", "End", "Direction", "Type"],
        how="left",
    )
    return df


# ==============================================================================
# HELPER: FUNCTIONAL CATEGORY
# ==============================================================================
def assign_functional_category(gene_name: str) -> str:
    if not gene_name or str(gene_name).strip() in ("", "nan", "HP"):
        return "Unknown"
    g = str(gene_name).lower()
    for category, prefixes in FUNCTIONAL_CATEGORIES.items():
        if any(g.startswith(p.lower()) for p in prefixes):
            return category
    return "Other"


def split_merged_tokens(value: str) -> list[str]:
    parts = [x.strip() for x in str(value).split(" \\ ") if x.strip()]
    return [re.sub(r"#+", "", p).strip() for p in parts if re.sub(r"#+", "", p).strip()]


def compute_priority_confidence(row: pd.Series, total_tools: int) -> float:
    seq_tokens = split_merged_tokens(row.get("Sequence_Name", ""))
    prod_tokens = [p.strip() for p in str(row.get("Product", "")).split(" \\ ") if p.strip()]
    gene_tokens = [g.strip() for g in str(row.get("Gene", "")).split(" \\ ") if g.strip() and g.strip().lower() != "nan"]

    tool_support = min(len(set(seq_tokens)), max(total_tools, 1))
    product_agreement = 1.0 / len(set(prod_tokens)) if prod_tokens else 0.0
    gene_agreement = 1.0 / len(set(gene_tokens)) if gene_tokens else 0.0
    partial_boundary = ("#" in str(row.get("Sequence_Name", ""))) or ("##" in str(row.get("Sequence_Name", "")))

    score = 10.0 + (tool_support / max(total_tools, 1)) * 55.0
    score += 0.0 if partial_boundary else 15.0
    score += product_agreement * 15.0
    score += gene_agreement * 5.0
    if tool_support >= min(3, max(total_tools, 1)) and not partial_boundary:
        score += 10.0
    if str(row.get("Status", "")).strip() in {"Overlap", "Short-CDS", "Long-CDS"}:
        score -= 10.0
    return round(min(100.0, max(0.0, score)), 2)


def best_final_product(row: pd.Series) -> str:
    consensus = str(row.get("Product-Consensus", "")).strip()
    if consensus and consensus.lower() != "nan":
        return consensus
    products = [p.strip() for p in str(row.get("Product", "")).split(" \\ ") if p.strip()]
    for product in products:
        if product.lower() not in {"hypothetical protein", "hypothetical_protein", "nan"}:
            return product
    return products[-1] if products else ""


# ==============================================================================
# STEP 1: GFF PARSING
# ==============================================================================
def parse_gff_attributes(attr_string: str) -> dict:
    attrs = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if "=" in item:
            key, _, val = item.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs


def is_ncbi_gff(gff_path: str) -> bool:
    """
    Detect NCBI PGAP GFF by checking for the NCBI annotwriter header
    or the presence of 'gene' feature type with 'gbkey=Gene' attribute.
    """
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "NCBI annotwriter" in line or "processor NCBI" in line:
                return True
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 9 and fields[2] == "gene" and "gbkey=Gene" in fields[8]:
                return True
            break
    return False


def preprocess_ncbi_gff(gff_path: str, stem: str, tool: str) -> list:
    """
    NCBI PGAP GFF3 uses a parent-child hierarchy that causes double-counting:
      gene  → CDS      (protein-coding: gene record + CDS record at same coords)
      gene  → exon     (rRNA/tRNA: annotation lives on exon, not on gene)
      pseudogene       (standalone record with pseudo=true, no child CDS)

    This pre-processor collapses each hierarchy into ONE clean record:
      gene + CDS pair   → keep CDS only  (has product, locus_tag, GO terms)
      gene + exon pair  → keep exon, remap Type from parent's gbkey (rRNA/tRNA)
      pseudogene        → emit as CDS with Pseudogene=True
      region/sequence_feature/direct_repeat/riboswitch → discard

    Also extracts NCBI-specific rich attributes:
      go_process, go_function, inference, protein_id
    and appends them to the Product field for downstream use.
    """
    # First pass: index all rows
    rows_by_id = {}
    children   = {}   # parent_id → list of child (fields, attrs)

    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs  = parse_gff_attributes(fields[8])
            row_id = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            if row_id:
                rows_by_id[row_id] = (fields, attrs)
            if parent:
                children.setdefault(parent, []).append((fields, attrs))

    # Discard types that are pure scaffolding / not biologically meaningful
    DISCARD_TYPES = {
        "region", "sequence_feature", "direct_repeat",
        "riboswitch", "SRP_RNA", "RNase_P_RNA", "antisense_RNA",
    }

    emitted = set()
    results = []

    def make_record(fields, attrs, type_override, pseudo, product_override=""):
        try:
            seq_id    = fields[0]
            start     = int(fields[3])
            end       = int(fields[4])
            length    = end - start + 1
            direction = fields[6]
            gene      = attrs.get("gene", "")
            locus_tag = attrs.get("locus_tag", "")
            product   = product_override or attrs.get("product", "")

            # Append GO terms to product for richer annotation context
            go_parts = []
            if attrs.get("go_process"):
                # GO terms are pipe-separated: "term|id||evidence"
                go_terms = [p.split("|")[0] for p in attrs["go_process"].split(",")]
                go_parts.extend(go_terms[:2])   # keep at most 2 to avoid bloat
            extra = "; ".join(go_parts)
            if extra:
                product = f"{product} [{extra}]" if product else extra

            return (stem, tool, seq_id, type_override, start, end,
                    length, direction, gene, product, locus_tag, pseudo)
        except (ValueError, IndexError):
            return None

    # Second pass: emit one record per locus
    for row_id, (fields, attrs) in rows_by_id.items():
        ftype = fields[2]

        if ftype in DISCARD_TYPES:
            continue

        # gene: skip — handled via CDS/exon child below
        if ftype == "gene":
            continue

        # pseudogene: standalone, no CDS child — emit directly
        if ftype == "pseudogene":
            rec = make_record(fields, attrs, type_override="CDS", pseudo=True)
            if rec:
                results.append(rec)
            emitted.add(row_id)
            continue

        # exon: remap type from parent's gbkey (rRNA, tRNA etc.)
        if ftype == "exon":
            parent_id  = attrs.get("Parent", "")
            parent_row = rows_by_id.get(parent_id)
            if parent_row:
                parent_attrs = parent_row[1]
                # gbkey on RNA parent gives us the correct type: rRNA, tRNA etc.
                gbkey = parent_attrs.get("gbkey", "ncRNA")
                rec   = make_record(fields, attrs, type_override=gbkey, pseudo=False)
                if rec:
                    results.append(rec)
            emitted.add(row_id)
            continue

        # CDS: the main protein-coding record
        if ftype == "CDS" and row_id not in emitted:
            pseudo = attrs.get("pseudo", "").lower() in ("true", "1", "yes")
            rec = make_record(fields, attrs, type_override="CDS", pseudo=pseudo)
            if rec:
                results.append(rec)
            emitted.add(row_id)
            continue

    return results


def parse_standard_gff(gff_path: str, stem: str, tool: str) -> list:
    """
    Standard GFF3 parser for Bakta and Prokka.
    Each line is a self-contained record — no parent-child collapsing needed.
    """
    results = []
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            try:
                seq_id    = fields[0]
                type_     = fields[2]
                start     = int(fields[3])
                end       = int(fields[4])
                length    = end - start + 1
                direction = fields[6]
                attrs     = parse_gff_attributes(fields[8])

                gene      = attrs.get("gene", "")
                product   = attrs.get("product", "")
                locus_tag = attrs.get("locus_tag", "")
                pseudo    = attrs.get("pseudo", "").lower() in ("true", "1", "yes")

                results.append((stem, tool, seq_id, type_, start, end,
                                 length, direction, gene, product, locus_tag, pseudo))

            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping malformed GFF line in {gff_path}: {e}")

    return results


def process_gff_files_in_folder(input_folder: str, output_folder: str) -> list:
    gff_files = sorted([
        f for f in os.listdir(input_folder)
        if f.endswith(".gff") or f.endswith(".gff3")
    ])
    if not gff_files:
        raise ValueError(f"No GFF/GFF3 files found in: {input_folder}")

    os.makedirs(output_folder, exist_ok=True)
    converted_stems = []

    for filename in gff_files:
        gff_path = os.path.join(input_folder, filename)
        stem     = Path(filename).stem
        tool     = detect_tool(filename)

        # Route to correct parser based on GFF source
        if is_ncbi_gff(gff_path):
            records = preprocess_ncbi_gff(gff_path, stem, tool)
            logger.info(f"  Parsed {filename}: {len(records)} features  (tool=ncbi/pgap, NCBI pre-processor applied)")
        else:
            records = parse_standard_gff(gff_path, stem, tool)
            logger.info(f"  Parsed {filename}: {len(records)} features  (tool={tool})")

        if not records:
            logger.warning(f"No records parsed from {filename} — skipping")
            continue

        df = pd.DataFrame(records, columns=[
            "Sequence_Name", "Tool", "Sequence_ID", "Type",
            "Start", "End", "Length", "Direction",
            "Gene", "Product", "Original_Locus_Tag", "Pseudogene",
        ])

        # Write per-tool intermediate (temp, not kept in final output)
        tmp_path = os.path.join(output_folder, f"_tmp_{stem}.xlsx")
        df.to_excel(tmp_path, index=False)
        converted_stems.append(stem)

    return converted_stems


# ==============================================================================
# STEP 2: MERGE
# ==============================================================================
def merge_excel(output_folder: str, converted_stems: list) -> pd.DataFrame:
    frames = []
    for stem in converted_stems:
        path = os.path.join(output_folder, f"_tmp_{stem}.xlsx")
        if os.path.exists(path):
            frames.append(pd.read_excel(path))

    if not frames:
        raise FileNotFoundError("No tool Excel files to merge.")

    merged = pd.concat(frames, ignore_index=True)
    merged.drop_duplicates(inplace=True)
    merged.sort_values(["Sequence_Name", "Sequence_ID", "Start"], inplace=True)
    logger.info(f"Step 2: Merged {len(frames)} tool files → {len(merged)} total features")
    return merged


# ==============================================================================
# STEP 3: FILTER VALID FEATURE TYPES
# ==============================================================================
def filter_feature_types(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)
    filtered = df[df["Type"].isin(VALID_FEATURE_TYPES)].copy()
    removed_types = df[~df["Type"].isin(VALID_FEATURE_TYPES)]["Type"].value_counts()

    logger.info(f"Step 3: {len(filtered)} kept, {before - len(filtered)} removed")
    for ftype, count in removed_types.items():
        logger.info(f"  Removed: {ftype} ({count})")

    return filtered


# ==============================================================================
# STEP 4: CLEAN GENE NAMES
# ==============================================================================
def clean_gene_name(gene_name) -> str:
    if pd.isna(gene_name) or str(gene_name).strip() in ("", "nan"):
        return ""
    name = str(gene_name).strip()
    # Remove only tool-appended trailing _1, _2 … suffixes
    name = re.sub(r"_\d+$", "", name)
    return name


# ==============================================================================
# STEP 5: GROUP V1 — IDENTICAL COORDINATES + STRAND
# ==============================================================================
def grouping_v1(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)
    grouped = df.groupby(
        ["Sequence_ID", "Type", "Start", "End", "Length", "Direction"],
        sort=False,
    ).agg({
        "Sequence_Name":     lambda x: " \\ ".join(x.fillna("").astype(str)),
        "Tool":              lambda x: " \\ ".join(x.fillna("").astype(str)),
        "Gene":              lambda x: " \\ ".join(x.fillna("").astype(str)),
        "Product":           lambda x: " \\ ".join(x.fillna("").astype(str)),
        "Original_Locus_Tag":lambda x: " \\ ".join(x.fillna("").astype(str)),
        "Pseudogene":        lambda x: any(x.fillna(False)),
    }).reset_index()

    grouped.sort_values(["Sequence_ID", "Start"], inplace=True)
    logger.info(f"Step 5: V1 grouping {before} → {len(grouped)} unique loci")
    return grouped


# ==============================================================================
# STEP 6: GROUP V2 — BOUNDARY RESOLUTION
# ==============================================================================
def grouping_v2(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(
        ["Sequence_ID", "Type", "Direction", "Start", "End"]
    ).reset_index(drop=True)

    result_rows = []
    used = set()

    for i, row_i in df.iterrows():
        if i in used:
            continue

        mask = (
            (~df.index.isin(used)) &
            (df.index != i) &
            (df["Sequence_ID"]  == row_i["Sequence_ID"]) &
            (df["Type"]         == row_i["Type"]) &
            (df["Direction"]    == row_i["Direction"]) &
            (
                ((df["Start"] == row_i["Start"]) & (df["End"] != row_i["End"])) |
                ((df["End"]   == row_i["End"])   & (df["Start"] != row_i["Start"]))
            )
        )
        partners = df[mask]

        if partners.empty:
            result_rows.append(row_i.to_dict())
            used.add(i)
            continue

        group = pd.concat([row_i.to_frame().T, partners]).reset_index(drop=True)

        coord_support = Counter()
        for _, r in group.iterrows():
            n = str(r["Sequence_Name"]).count("\\") + 1
            coord_support[(int(r["Start"]), int(r["End"]))] += n

        best_start, best_end = coord_support.most_common(1)[0][0]

        merged = {
            "Sequence_ID":      row_i["Sequence_ID"],
            "Type":             row_i["Type"],
            "Start":            best_start,
            "End":              best_end,
            "Length":           best_end - best_start + 1,
            "Direction":        row_i["Direction"],
            "Sequence_Name":    " \\ ".join(group["Sequence_Name"].fillna("").astype(str)),
            "Tool":             " \\ ".join(group["Tool"].fillna("").astype(str)),
            "Gene":             " \\ ".join(group["Gene"].fillna("").astype(str)),
            "Product":          " \\ ".join(group["Product"].fillna("").astype(str)),
            "Original_Locus_Tag": " \\ ".join(group["Original_Locus_Tag"].fillna("").astype(str)),
            "Pseudogene":       any(group["Pseudogene"]),
        }
        result_rows.append(merged)
        used.add(i)
        used.update(partners.index.tolist())

    result = pd.DataFrame(result_rows)
    result.sort_values(["Sequence_ID", "Start"], inplace=True)
    logger.info(f"Step 6: V2 boundary resolution → {len(result)} loci")
    return result


def grouping_v2_non_cds_near_match(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(
        ["Sequence_ID", "Type", "Direction", "Start", "End"]
    ).reset_index(drop=True)

    result_rows = []
    used = set()

    for i, row_i in df.iterrows():
        if i in used:
            continue

        if str(row_i["Type"]) == "CDS":
            result_rows.append(row_i.to_dict())
            used.add(i)
            continue

        partner_ids = []
        start_i = int(row_i["Start"])
        end_i = int(row_i["End"])
        len_i = max(int(row_i["Length"]), 1)

        for j, row_j in df.iterrows():
            if j <= i or j in used:
                continue
            if (row_j["Sequence_ID"] != row_i["Sequence_ID"] or
                    row_j["Type"] != row_i["Type"] or
                    row_j["Direction"] != row_i["Direction"]):
                continue

            start_j = int(row_j["Start"])
            end_j = int(row_j["End"])
            len_j = max(int(row_j["Length"]), 1)
            overlap = min(end_i, end_j) - max(start_i, start_j) + 1
            if overlap <= 0:
                continue

            recip_i = overlap / len_i
            recip_j = overlap / len_j
            if recip_i >= 0.95 and recip_j >= 0.95:
                partner_ids.append(j)

        if not partner_ids:
            result_rows.append(row_i.to_dict())
            used.add(i)
            continue

        group = pd.concat([row_i.to_frame().T, df.loc[partner_ids]]).reset_index(drop=True)
        best_idx = group["Length"].astype(int).idxmax()
        best_row = group.loc[best_idx]

        merged = {
            "Sequence_ID":        best_row["Sequence_ID"],
            "Type":               best_row["Type"],
            "Start":              int(best_row["Start"]),
            "End":                int(best_row["End"]),
            "Length":             int(best_row["Length"]),
            "Direction":          best_row["Direction"],
            "Sequence_Name":      " \\ ".join(group["Sequence_Name"].fillna("").astype(str)),
            "Tool":               " \\ ".join(group["Tool"].fillna("").astype(str)),
            "Gene":               " \\ ".join(group["Gene"].fillna("").astype(str)),
            "Product":            " \\ ".join(group["Product"].fillna("").astype(str)),
            "Original_Locus_Tag": " \\ ".join(group["Original_Locus_Tag"].fillna("").astype(str)),
            "Pseudogene":         any(group["Pseudogene"]),
        }
        result_rows.append(merged)
        used.add(i)
        used.update(partner_ids)

    result = pd.DataFrame(result_rows)
    result.sort_values(["Sequence_ID", "Start"], inplace=True)
    logger.info(f"Step 6b: non-CDS near-match merge → {len(result)} loci")
    return result


# ==============================================================================
# STEP 7: GROUP V3 — OVERLAP DETECTION (O(n) sorted scan)
# ==============================================================================
def grouping_v3(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(
        ["Sequence_ID", "Type", "Direction", "Start", "End"]
    ).reset_index(drop=True)
    df["Status"] = ""

    for i in range(len(df) - 1):
        r1 = df.iloc[i]
        r2 = df.iloc[i + 1]
        if (r1["Sequence_ID"] != r2["Sequence_ID"] or
                r1["Type"] != r2["Type"] or
                r1["Direction"] != r2["Direction"]):
            continue
        # Strict nesting
        if (r1["Start"] < r2["Start"] and r2["End"] < r1["End"]) or \
           (r2["Start"] < r1["Start"] and r1["End"] < r2["End"]):
            df.at[i,     "Status"] = "Overlap"
            df.at[i + 1, "Status"] = "Overlap"

    overlap_n = (df["Status"] == "Overlap").sum()
    logger.info(f"Step 7: V3 overlap detection → {overlap_n} overlapping features flagged")
    return df


# ==============================================================================
# STEP 8: BEST GENE SELECTION
# ==============================================================================
def find_best_gene(row) -> str:
    genes = [g.strip() for g in str(row.get("Gene", "")).split(" \\ ")
             if g.strip() and g.strip().lower() != "nan"]
    tools = [t.strip().lower() for t in str(row.get("Tool", "")).split(" \\ ")]

    if not genes:
        return ""

    # 1. Consensus: gene in 2+ entries
    counts = Counter(genes)
    consensus = [g for g, c in counts.items() if c > 1]
    if consensus:
        return consensus[0]

    # 2. Tool priority fallback: Bakta > Prokka
    for preferred in TOOL_PRIORITY:
        for i, tool in enumerate(tools):
            if preferred in tool and i < len(genes) and genes[i]:
                return genes[i]

    # 3. Last resort: first non-empty
    return genes[0]


# ==============================================================================
# STEP 9: FLAG HYPOTHETICAL PROTEINS + PSEUDOGENES
# ==============================================================================
def flag_hp_and_pseudogenes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # Pseudogenes: preserve gene name with _pseudo suffix, set special status
    pseudo_mask = df["Pseudogene"] == True
    df.loc[pseudo_mask, "Status"] = "Pseudogene"
    df.loc[pseudo_mask, "Best-Gene"] = (
        df.loc[pseudo_mask, "Best-Gene"].apply(
            lambda g: (str(g) + "_pseudo") if g and str(g) not in ("", "nan") else "pseudo"
        )
    )

    # HP: CDS with no gene name and not pseudogene
    is_cds    = df["Type"] == "CDS"
    no_gene   = df["Best-Gene"].astype(str).str.strip().isin(["", "nan"]) | df["Best-Gene"].isna()
    not_pseudo = df["Status"] != "Pseudogene"
    df.loc[is_cds & no_gene & not_pseudo, "Best-Gene"] = "HP"

    hp_n     = (df["Best-Gene"] == "HP").sum()
    pseudo_n = (df["Status"] == "Pseudogene").sum()
    logger.info(f"Step 9: {hp_n} hypothetical proteins, {pseudo_n} pseudogenes flagged")
    return df


# ==============================================================================
# STEP 10: HP RESCUE VIA PRODUCT STRINGS
# ==============================================================================
def rescue_hp_by_product(df: pd.DataFrame) -> pd.DataFrame:
    """
    For HP-flagged CDS rows, attempt to recover a functional name
    from the product strings contributed by the annotation tools.

    Strategy:
      1. Collect all non-hypothetical product strings across tools
      2. If both tools agree on the same product → high confidence rescue
      3. If only one tool has a non-hypothetical product → rescue with
         ScreenMaybe status
      4. Try to extract a short gene-like token from the product
         (e.g. "RecA protein" → "recA")
    """
    df = df.copy()
    hp_mask = (df["Best-Gene"] == "HP") & (df["Type"] == "CDS")
    rescued = 0

    for idx, row in df[hp_mask].iterrows():
        all_products = [
            p.strip() for p in str(row.get("Product", "")).split(" \\ ")
            if p.strip() and p.strip().lower() not in ("nan", "", "hypothetical protein",
                                                         "hypothetical_protein")
        ]
        if not all_products:
            continue

        product_counts = Counter(p.lower() for p in all_products)
        best_product   = Counter({p: c for p, c in product_counts.items()}).most_common(1)[0][0]
        agreed         = product_counts[best_product] > 1

        # Try to extract gene symbol from product name
        # Patterns: "RecA recombinase", "FtsZ cell division protein", "recA protein"
        gene_match = re.match(
            r"^([a-zA-Z]{2,4}[A-Z0-9]?)\s+(protein|recombinase|synthase|kinase|"
            r"reductase|dehydrogenase|transferase|hydrolase|family|domain)",
            best_product,
            re.IGNORECASE,
        )

        generic_product = bool(re.search(
            r"(domain-containing|family protein|family|putative|uncharacterized|"
            r"hypothetical|predicted|like protein|like family)",
            best_product,
            re.IGNORECASE,
        ))
        generic_gene_heads = {
            "phage", "lysis", "lysin", "coat", "ovule", "coa",
            "transferase", "protein", "putative"
        }

        if gene_match and not generic_product:
            extracted = gene_match.group(1)
            # Normalise to bacterial convention: first letters lower, last upper if 4-char
            extracted = extracted[0].lower() + extracted[1:]
            if extracted.lower() not in generic_gene_heads:
                df.at[idx, "Best-Gene"] = extracted
                df.at[idx, "Status"]    = "Good" if agreed else "ScreenMaybe"
                rescued += 1
                continue

        # Keep generic product labels as product-only evidence rather than
        # converting them into gene symbols like "big2", "aaa", or "phage".
        df.at[idx, "Product-Consensus"] = best_product
        df.at[idx, "Status"] = "HP-ProductRescued" if agreed else "ScreenMaybe"
        rescued += 1

    logger.info(f"Step 10 (HP Rescue): {rescued} HP entries recovered from product strings")
    return df


# ==============================================================================
# STEP 11: SINGLE GENE + STATUS PER LOCUS
# ==============================================================================
def update_columns(row):
    best   = str(row.get("Best-Gene", "")).strip()
    is_empty = best in ("", "nan") or pd.isna(row.get("Best-Gene"))

    if is_empty:
        tools = str(row.get("Tool", "")).lower()
        genes = [g.strip() for g in str(row.get("Gene", "")).split("\\")
                 if g.strip() and g.strip().lower() != "nan"]
        for preferred in TOOL_PRIORITY:
            if preferred in tools and genes:
                row["Best-Gene"] = genes[0]
                if str(row.get("Status", "")).strip() == "":
                    row["Status"] = "ScreenMaybe"
                break

    status = str(row.get("Status", "")).strip()
    if status in ("", "nan") or pd.isna(row.get("Status")):
        row["Status"] = "Good"

    return row


# ==============================================================================
# STEP 12: GENE SYNONYMS
# ==============================================================================
def extract_gene_synonyms(row) -> str:
    all_genes  = {g.strip().strip('"') for g in str(row.get("Gene", "")).split("\\")
                  if g.strip() and g.strip().lower() not in ("nan", "")}
    best_genes = {g.strip().strip('"') for g in str(row.get("Best-Gene", "")).split("\\")
                  if g.strip() and g.strip().lower() not in ("nan", "")}
    return " \\ ".join(sorted(all_genes - best_genes))


# ==============================================================================
# STEP 13: CONSENSUS SCORING
# ==============================================================================
def consensus_scoring(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Source-of-Gene-Name"] = ""
    df["Consensus-Level"]    = ""
    df["Confidence-Score"]   = 0.0
    df["Product-Consensus"]  = df.get("Product-Consensus", "")

    # Detect total number of distinct tools actually present in this run
    # so labels say Consensus_2/2 (not Consensus_2/3) when only 2 tools were used
    all_tools = set()
    for tools_str in df["Tool"].dropna():
        for t in str(tools_str).split(" \\ "):
            t = t.strip().lower()
            if t and t != "nan":
                all_tools.add(t)
    total_tools = max(min(len(all_tools), MAX_TOOL_SUPPORT), 1)
    logger.info(f"  Detected {total_tools} tool(s) for consensus labeling: {', '.join(sorted(all_tools))}")

    for idx, row in df.iterrows():
        best_gene = str(row.get("Best-Gene", "")).strip()
        if not best_gene or best_gene == "nan":
            continue

        seq_names = [s.strip() for s in str(row.get("Sequence_Name", "")).split(" \\ ")]
        genes     = [g.strip() for g in str(row.get("Gene", "")).split(" \\ ")]
        products  = [p.strip() for p in str(row.get("Product", "")).split(" \\ ")]

        if best_gene == "HP":
            sources = [s.replace("#", "").strip() for s in seq_names]
        else:
            sources = [
                seq_names[i] for i, g in enumerate(genes)
                if g == best_gene and i < len(seq_names)
            ]

        source_str = ", ".join(filter(None, sources)).rstrip(", ")
        df.at[idx, "Source-of-Gene-Name"] = source_str

        num_tools = min(len([s for s in sources if s.strip()]), total_tools)

        # Labels are relative to actual tool count in this run
        if num_tools >= total_tools:
            # All tools agree
            df.at[idx, "Consensus-Level"]  = f"Consensus_{total_tools}/{total_tools}"
        elif num_tools >= 2:
            # Partial consensus (at least 2 tools agree)
            df.at[idx, "Consensus-Level"]  = f"Consensus_{num_tools}/{total_tools}"
        else:
            df.at[idx, "Consensus-Level"]  = "Single-Tool"

        # Product consensus
        if not str(df.at[idx, "Product-Consensus"]).strip():
            valid_products = [
                p for p in products
                if p and p.lower() not in ("nan", "", "hypothetical protein",
                                            "hypothetical_protein")
            ]
            if valid_products:
                df.at[idx, "Product-Consensus"] = Counter(valid_products).most_common(1)[0][0]
            elif best_gene == "HP":
                df.at[idx, "Product-Consensus"] = "hypothetical protein"

    df["Confidence-Score"] = df.apply(
        lambda row: compute_priority_confidence(row, total_tools), axis=1
    )

    return df


# ==============================================================================
# STEP 14: LENGTH OUTLIER + STRAND BIAS FLAGS
# ==============================================================================
def flag_length_outliers(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cds = df["Type"] == "CDS"
    # < 90 bp (< 30 aa) or > 15 000 bp (> 5000 aa) — flag for review
    df.loc[cds & (df["Length"] < 90),    "Status"] = "Short-CDS"
    df.loc[cds & (df["Length"] > 15000), "Status"] = "Long-CDS"
    short_n = (df["Status"] == "Short-CDS").sum()
    long_n  = (df["Status"] == "Long-CDS").sum()
    if short_n or long_n:
        logger.warning(f"Length outliers: {short_n} Short-CDS (<90 bp), {long_n} Long-CDS (>15 kb)")
    return df


def check_strand_bias(df: pd.DataFrame):
    cds = df[df["Type"] == "CDS"]
    for strand in ["+", "-"]:
        sub = cds[cds["Direction"] == strand]
        if len(sub) == 0:
            continue
        hp_rate = (sub["Best-Gene"] == "HP").sum() / len(sub) * 100
        if hp_rate > 35:
            logger.warning(
                f"High HP rate on {strand} strand: {hp_rate:.1f}% "
                f"— possible contig orientation issue in assembly"
            )


# ==============================================================================
# OUTPUT: GFF3
# ==============================================================================
def write_gff(df: pd.DataFrame, output_folder: str, prefix: str):
    """
    Write a clean GFF3 file compatible with Geneious, SnapGene, and standard
    GFF3 parsers. Includes all features (not just high-confidence).
    Confidence level and score are written as GFF attributes.
    """
    out_path = os.path.join(output_folder, f"{prefix}.gff3")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"# Generated by Unified Annotation Pipeline v3  {datetime.now().date()}\n")
        fh.write(f"# Locus tag prefix: {prefix}\n")

        for _, row in df.sort_values(["Sequence_ID", "Start"]).iterrows():
            seq_id    = row["Sequence_ID"]
            ftype     = row["Type"]
            start     = int(row["Start"])
            end       = int(row["End"])
            strand    = row["Direction"] if row["Direction"] in ("+", "-") else "."
            locus_tag = row.get("Locus_Tag", "")
            best_gene = str(row.get("Best-Gene", "")).strip()
            product   = str(row.get("Product-Consensus", "")).strip()
            conf      = str(row.get("Consensus-Level", "")).strip()
            score_val = str(row.get("Confidence-Score", "."))
            status    = str(row.get("Status", "")).strip()

            # Build attributes
            attrs = []
            if locus_tag:
                attrs.append(f"ID={locus_tag}")
                attrs.append(f"locus_tag={locus_tag}")
            if best_gene and best_gene not in ("HP", "nan", ""):
                attrs.append(f"gene={best_gene}")
            if product and product not in ("nan", ""):
                attrs.append(f"product={product}")
            if conf:
                attrs.append(f"consensus_level={conf}")
            if status:
                attrs.append(f"status={status}")
            if row.get("Pseudogene"):
                attrs.append("pseudo=true")

            attr_str = ";".join(attrs) if attrs else "."

            fh.write(
                f"{seq_id}\tUnifiedPipeline\t{ftype}\t{start}\t{end}\t"
                f"{score_val}\t{strand}\t.\t{attr_str}\n"
            )

    logger.info(f"  GFF3 written → {out_path}  ({len(df)} features)")
    return out_path


# ==============================================================================
# OUTPUT: GenBank (GBK)
# ==============================================================================
def write_gbk(df: pd.DataFrame, genome: dict, output_folder: str,
              prefix: str):
    """
    Write a GenBank flat file (.gbk) compatible with Geneious, SnapGene,
    Addgene, and NCBI Sequin.

    One LOCUS record per contig/chromosome in the genome FASTA.
    Features are written in genomic order within each record.
    """
    out_path = os.path.join(output_folder, f"{prefix}.gbk")
    today    = datetime.now().strftime("%d-%b-%Y").upper()

    with open(out_path, "w", encoding="utf-8") as fh:
        for seq_id, seq in genome.items():
            seq_len  = len(seq)
            features = df[df["Sequence_ID"] == seq_id].sort_values("Start")

            # --- LOCUS line ---
            fh.write(
                f"LOCUS       {seq_id:<16} {seq_len:>10} bp    DNA     linear   BCT {today}\n"
            )
            fh.write(f"DEFINITION  .\n")
            fh.write(f"ACCESSION   {seq_id}\n")
            fh.write(f"VERSION     {seq_id}\n")
            fh.write(f"KEYWORDS    .\n")
            fh.write(f"SOURCE      .\n")
            fh.write(f"  ORGANISM  .\n")
            fh.write(f"            Bacteria.\n")
            fh.write("FEATURES             Location/Qualifiers\n")

            # --- source feature ---
            fh.write(f"     source          1..{seq_len}\n")
            fh.write(f"                     /mol_type=\"genomic DNA\"\n")

            # --- annotation features ---
            for _, row in features.iterrows():
                ftype     = row["Type"]
                start     = int(row["Start"])
                end       = int(row["End"])
                strand    = row["Direction"]
                locus_tag = str(row.get("Locus_Tag", "")).strip()
                best_gene = str(row.get("Best-Gene", "")).strip()
                product   = str(row.get("Product-Consensus", "")).strip()
                conf      = str(row.get("Consensus-Level", "")).strip()
                status    = str(row.get("Status", "")).strip()
                is_pseudo = row.get("Pseudogene", False)
                synonyms  = str(row.get("Gene-Putative-Synonyms", "")).strip()

                # Location string
                if strand == "-":
                    location = f"complement({start}..{end})"
                else:
                    location = f"{start}..{end}"

                # GBK spec: feature type col is 16 chars, then location.
                # If type is exactly 16 chars (e.g. "regulatory_region") we must
                # still leave at least one space before the location.
                ftype_field = f"{ftype:<16}"
                if len(ftype) >= 16:
                    ftype_field = ftype + " "   # force at least one space
                fh.write(f"     {ftype_field}{location}\n")

                def write_qualifier(key, value):
                    if not value or str(value).strip() in ("nan", "", "None"):
                        return
                    # GenBank spec: qualifier lines max 79 chars total.
                    # Prefix "                     /key="" is 21 + len(key) + 2 chars.
                    # Continuation lines start with 21 spaces.
                    # Closing quote must appear on the LAST line.
                    full = f'                     /{key}="{value}"'
                    if len(full) <= 79:
                        fh.write(full + "\n")
                        return
                    # How many chars fit on the first line after the opening quote
                    prefix_len = 21 + len(key) + 2   # "                     /key=\""
                    first_chunk = 79 - prefix_len
                    chunks = []
                    remaining = value
                    # First line
                    chunks.append(remaining[:first_chunk])
                    remaining = remaining[first_chunk:]
                    # Continuation lines: 21 spaces prefix, 58 chars of value
                    while remaining:
                        chunks.append(remaining[:58])
                        remaining = remaining[58:]
                    # Write first line without closing quote
                    fh.write(f'                     /{key}="{chunks[0]}\n')
                    # Write middle lines
                    for chunk in chunks[1:-1]:
                        fh.write(f"                     {chunk}\n")
                    # Write last line with closing quote
                    fh.write(f'                     {chunks[-1]}"\n')

                if locus_tag:
                    write_qualifier("locus_tag", locus_tag)
                if best_gene and best_gene not in ("HP", "nan", ""):
                    write_qualifier("gene", best_gene)
                if product and product not in ("nan", ""):
                    write_qualifier("product", product)
                elif ftype == "CDS":
                    fh.write('                     /product="hypothetical protein"\n')
                if ftype == "CDS":
                    write_qualifier("protein_id", locus_tag)
                if conf:
                    write_qualifier("note", f"consensus_level={conf}")
                if status and status not in ("Good", "nan", ""):
                    write_qualifier("note", f"status={status}")
                if synonyms and synonyms not in ("nan", ""):
                    write_qualifier("note", f"gene_synonyms={synonyms}")
                if is_pseudo:
                    fh.write("                     /pseudo\n")

                # CDS translation qualifier
                if ftype == "CDS" and not is_pseudo:
                    nt_seq = extract_subsequence(genome, seq_id, start, end, strand)
                    if nt_seq:
                        protein = translate(nt_seq)
                        if protein:
                            # Write /translation="..." wrapped at 58 chars per line
                            # with closing quote on the final line (GenBank spec)
                            chunks = [protein[i:i+58] for i in range(0, len(protein), 58)]
                            if len(chunks) == 1:
                                fh.write(f'                     /translation="{chunks[0]}"\n')
                            else:
                                fh.write(f'                     /translation="{chunks[0]}\n')
                                for chunk in chunks[1:-1]:
                                    fh.write(f"                     {chunk}\n")
                                fh.write(f'                     {chunks[-1]}"\n')

            # --- ORIGIN (sequence) ---
            fh.write("ORIGIN\n")
            for i in range(0, seq_len, 60):
                chunk = seq[i:i + 60]
                grouped = " ".join(chunk[j:j+10] for j in range(0, len(chunk), 10))
                fh.write(f"      {i + 1:>9} {grouped}\n")
            fh.write("//\n")

    logger.info(f"  GBK written → {out_path}  ({len(genome)} records)")
    return out_path


# ==============================================================================
# OUTPUT: FAA (protein sequences, high-confidence CDS only)
# ==============================================================================
def write_faa(df: pd.DataFrame, genome: dict, output_folder: str, prefix: str):
    """
    Write protein FASTA for all valid CDS features.
    Excludes pseudogenes and short CDS (<90 bp) only.
    Header format: ><locus_tag> <gene> <product> [<confidence_level>]
    """
    out_path = os.path.join(output_folder, f"{prefix}.faa")
    written  = 0

    high_conf_cds = df[
        (df["Type"] == "CDS") &
        (df["Pseudogene"] != True) &
        (df["Status"] != "Short-CDS")
    ].sort_values(["Sequence_ID", "Start"])

    skipped_no_seq = 0
    skipped_no_aa  = 0

    with open(out_path, "w", encoding="utf-8") as fh:
        for _, row in high_conf_cds.iterrows():
            locus_tag = str(row.get("Locus_Tag", "")).strip()
            gene      = str(row.get("Best-Gene", "")).strip()
            product   = str(row.get("Product-Consensus", "hypothetical protein")).strip()
            conf      = str(row.get("Consensus-Level", "")).strip()
            seq_id    = row["Sequence_ID"]
            start     = int(row["Start"])
            end       = int(row["End"])
            strand    = row["Direction"]

            nt_seq = extract_subsequence(genome, seq_id, start, end, strand)
            if not nt_seq:
                skipped_no_seq += 1
                continue
            protein = translate(nt_seq)
            if not protein:
                skipped_no_aa += 1
                logger.warning(
                    f"  FAA: translation failed for {locus_tag} ({gene}) "
                    f"{seq_id}:{start}..{end}({strand}) — skipped"
                )
                continue

            # Sanitise: '>' inside product/gene names breaks FASTA parsing
            gene_safe    = gene.replace(">", "-")
            product_safe = product.replace(">", "-")
            header = f">{locus_tag} {gene_safe} {product_safe} [{conf}]"
            fh.write(header + "\n")
            # Wrap at 60 chars per line
            for i in range(0, len(protein), 60):
                fh.write(protein[i:i + 60] + "\n")
            written += 1

    logger.info(
        f"  FAA written → {out_path}  ({written} proteins"
        f"{f', {skipped_no_seq} skipped (no FASTA seq)' if skipped_no_seq else ''}"
        f"{f', {skipped_no_aa} skipped (translation failed)' if skipped_no_aa else ''})"
    )
    return out_path


# ==============================================================================
# OUTPUT: FNA (nucleotide CDS sequences, high-confidence only)
# ==============================================================================
def write_fna(df: pd.DataFrame, genome: dict, output_folder: str, prefix: str):
    """
    Write nucleotide FASTA for all valid CDS features.
    Excludes pseudogenes and short CDS (<90 bp) only.
    Header format: ><locus_tag> <gene> <start>..<end> <strand> [<confidence_level>]
    """
    out_path = os.path.join(output_folder, f"{prefix}.fna")
    written  = 0

    high_conf_cds = df[
        (df["Type"] == "CDS") &
        (df["Pseudogene"] != True) &
        (df["Status"] != "Short-CDS")
    ].sort_values(["Sequence_ID", "Start"])

    with open(out_path, "w", encoding="utf-8") as fh:
        for _, row in high_conf_cds.iterrows():
            locus_tag = str(row.get("Locus_Tag", "")).strip()
            gene      = str(row.get("Best-Gene", "")).strip()
            conf      = str(row.get("Consensus-Level", "")).strip()
            seq_id    = row["Sequence_ID"]
            start     = int(row["Start"])
            end       = int(row["End"])
            strand    = row["Direction"]

            nt_seq = extract_subsequence(genome, seq_id, start, end, strand)
            if not nt_seq:
                continue

            # Sanitise: '>' inside gene names breaks FASTA parsing
            gene_safe = gene.replace(">", "-")
            header = f">{locus_tag} {gene_safe} {seq_id}:{start}..{end}({strand}) [{conf}]"
            fh.write(header + "\n")
            for i in range(0, len(nt_seq), 60):
                fh.write(nt_seq[i:i + 60] + "\n")
            written += 1

    logger.info(f"  FNA written → {out_path}  ({written} nucleotide sequences)")
    return out_path


# ==============================================================================
# OUTPUT: MASTER TABLE (Excel)
# ==============================================================================
def write_master_table(df: pd.DataFrame, output_folder: str):
    """
    Write the final Master_Table_Annotation.xlsx with all features and columns
    in a clean, readable order.
    """
    df = df.copy()
    if "Source-Chosen-Gene" in df.columns and "Source-of-Gene-Name" not in df.columns:
        df = df.rename(columns={"Source-Chosen-Gene": "Source-of-Gene-Name"})
    review_reason_pairs = df.apply(derive_export_review_reason, axis=1, result_type="expand")
    df["Review-Reason-ID"] = review_reason_pairs[0]
    df["Review-Reason"] = review_reason_pairs[1]

    col_order = [
        "Locus_Tag", "Sequence_ID", "Type", "Start", "End", "Length",
        "Direction", "Final-Product", "Best-Gene", "Gene-Putative-Synonyms",
        "Product-Consensus",
        "Consensus-Level", "Confidence-Score", "Status", "Review-Reason-ID",
        "Review-Reason", "Source-of-Gene-Name",
        "Functional-Category", "Pseudogene",
        "Tool", "Gene", "Product", "Original_Locus_Tag",
    ]
    export_column_names = {
        "Sequence_ID": "Contig_ID",
        "Best-Gene": "Consensus-Gene-Name",
        "Gene-Putative-Synonyms": "Likely-Gene-Synonyms",
        "Confidence-Score": "Bactowise-Confidence-Score",
        "Pseudogene": "PGAP_Predicted_Pseudogene",
        "Gene": "Gene_Names_from_Input_Gffs",
        "Product": "Product_from_Input_Gffs",
        "Original_Locus_Tag": "Locus_Tag_from_Input_Gffs",
    }
    # Only keep columns that actually exist
    cols = [c for c in col_order if c in df.columns]
    out_path = os.path.join(output_folder, "Master_Table_Annotation.xlsx")
    df_out = df[cols].sort_values(["Sequence_ID", "Start"]).reset_index(drop=True)
    column_meanings = {
        "Locus_Tag": "Final locus tag assigned by the pipeline in genomic order.",
        "Sequence_ID": "Reference sequence or contig identifier used for the final locus.",
        "Type": "Final feature type, such as CDS, tRNA, rRNA, or ncRNA.",
        "Start": "Final 1-based start coordinate of the locus.",
        "End": "Final 1-based end coordinate of the locus.",
        "Length": "Feature length in base pairs.",
        "Direction": "Final strand direction, either + or -.",
        "Final-Product": "Single user-facing product description shown in the final table.",
        "Best-Gene": "Final resolved consensus gene symbol. Left blank when no trustworthy gene name is assigned.",
        "Gene-Putative-Synonyms": "Alternative gene names detected across supporting annotations.",
        "Product-Consensus": "Resolved product description before final export formatting.",
        "Consensus-Level": "Support level across tools for the final locus, such as Consensus_3/3 or Single-Tool.",
        "Confidence-Score": "Bactowise review-priority score from 0 to 100; higher means cleaner multi-tool support.",
        "Status": "Final review status, such as Good, ReviewMaybe, Overlap, Short-CDS, or Pseudogene.",
        "Review-Reason-ID": "Compact code explaining why a row needs review or carries a warning state.",
        "Review-Reason": "Plain-language explanation of the review reason for flagged rows; blank for rows with no specific review reason.",
        "Source-of-Gene-Name": "Source annotation track(s) contributing the selected gene name. This can be blank for non-gene features or product-only rows.",
        "Functional-Category": "Broad functional grouping inferred from the final gene or product.",
        "Pseudogene": "PGAP/NCBI pseudogene prediction flag. Shown as YES when true and blank otherwise.",
        "Tool": "Raw contributing tool label(s) retained for provenance.",
        "Gene": "Raw gene name field(s) from the input GFF annotations before final resolution.",
        "Product": "Raw product field(s) from the input GFF annotations before final resolution.",
        "Original_Locus_Tag": "Raw locus tag(s) from the input GFF annotations.",
    }
    readme_rows = [
        {
            "Section": "Purpose",
            "Item": "Workbook",
            "Meaning": "This workbook contains the final harmonised Bactowise annotation table plus a guide to help users interpret each exported column.",
        },
        {
            "Section": "Purpose",
            "Item": "Master_Table",
            "Meaning": "Use the second sheet for filtering, sorting, and manual review of final loci after multi-tool grouping and reconciliation.",
        },
    ]
    for i, col in enumerate(cols):
        readme_rows.append(
            {
                "Section": "Columns" if i == 0 else "",
                "Item": export_column_names.get(col, col),
                "Meaning": column_meanings.get(col, "Column retained from the final pipeline output."),
            }
        )
    readme_df = pd.DataFrame(readme_rows)
    screenmaybe_hp_mask = (
        (df_out.get("Status", "").astype(str) == "ScreenMaybe") &
        (df_out.get("Best-Gene", "").astype(str) == "HP") &
        (~df_out.get("Product-Consensus", "").astype(str).str.lower().isin(["", "nan", "hypothetical protein", "hypothetical_protein"]))
    )
    df_out.loc[screenmaybe_hp_mask, "Best-Gene"] = ""
    if "Best-Gene" in df_out.columns:
        df_out.loc[df_out["Best-Gene"].astype(str).str.strip() == "HP", "Best-Gene"] = ""
    if "Pseudogene" in df_out.columns:
        pseudo_display = df_out["Pseudogene"].fillna(False).astype(str).str.lower().isin(
            ["true", "1", "yes"]
        )
        df_out["Pseudogene"] = pseudo_display.map({True: "YES", False: ""})
    df_out["Status"] = df_out["Status"].apply(display_status)
    df_out = df_out.rename(columns=export_column_names)
    csv_path = os.path.join(output_folder, "Master_Table_Annotation.csv")
    df_out.to_csv(csv_path, index=False)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        readme_df.to_excel(writer, sheet_name="README", index=False)
        df_out.to_excel(writer, sheet_name="Master_Table", index=False)

    wb = load_workbook(out_path)
    readme_ws = wb["README"]
    ws = wb["Master_Table"]
    header_font = Font(bold=True, size=12)
    readme_font = Font(size=12)
    italic_font = Font(italic=True, color="C00000")
    centered = Alignment(horizontal="center", vertical="center")
    wrapped = Alignment(horizontal="left", vertical="top", wrap_text=True)

    for cell in readme_ws[1]:
        cell.font = header_font
        cell.alignment = centered
    for row in readme_ws.iter_rows(min_row=2):
        for cell in row:
            cell.font = readme_font
            cell.alignment = wrapped
    readme_ws.column_dimensions["A"].width = 14
    readme_ws.column_dimensions["B"].width = 24
    readme_ws.column_dimensions["C"].width = 52

    gene_style_col_idxs = []
    for cell in ws[1]:
        cell.font = header_font
        cell.alignment = centered
        if cell.value in {"Consensus-Gene-Name", "Likely-Gene-Synonyms"}:
            gene_style_col_idxs.append(cell.column)
        col_letter = get_column_letter(cell.column)
        header_width = max(len(str(cell.value)) + 2, 12)
        if ws.column_dimensions[col_letter].width is None or ws.column_dimensions[col_letter].width < header_width:
            ws.column_dimensions[col_letter].width = header_width

    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = centered
        for col_idx in gene_style_col_idxs:
            row[col_idx - 1].font = italic_font

    ws.auto_filter.ref = ws.dimensions
    wb.save(out_path)
    logger.info(f"  Master table → {out_path}  ({len(df)} features)")
    logger.info(f"  Master table CSV → {csv_path}  ({len(df)} features)")
    return out_path


def write_manual_queue(df: pd.DataFrame, output_folder: str) -> tuple[str, str]:
    full_path = os.path.join(output_folder, "Manual-Curation-Queue_full.csv")
    high_path = os.path.join(output_folder, "Manual-Curation-Queue.csv")

    manual_df = df[df["Status"].isin(["Overlap", "ScreenMaybe"])].copy()
    manual_df.to_csv(full_path, index=False)

    if len(manual_df):
        manual_df["priority_overlap"] = (manual_df["Status"] == "Overlap").astype(int)
        manual_df["priority_low_conf"] = (manual_df["Confidence-Score"] < 65).astype(int)
        manual_df["priority_unknown_source"] = (
            manual_df["Source-of-Gene-Name"].fillna("").astype(str).str.strip() == ""
        ).astype(int)
        manual_df = manual_df.sort_values(
            by=["priority_overlap", "priority_low_conf", "priority_unknown_source", "Confidence-Score"],
            ascending=[False, False, False, True],
        )
        high_df = manual_df.head(800).drop(
            columns=["priority_overlap", "priority_low_conf", "priority_unknown_source"]
        )
    else:
        high_df = manual_df

    high_df.to_csv(high_path, index=False)
    logger.info(f"  Manual queue → {high_path}")
    logger.info(f"  Manual queue full → {full_path}")
    return high_path, full_path


def blank_screenmaybe_hp_gene_display(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    mask = (
        (df.get("Status", "").astype(str) == "ScreenMaybe") &
        (df.get("Best-Gene", "").astype(str) == "HP") &
        (~df.get("Product-Consensus", "").astype(str).str.lower().isin(
            ["", "nan", "hypothetical protein", "hypothetical_protein"]
        ))
    )
    df.loc[mask, "Best-Gene"] = ""
    return df


def display_value(value) -> str:
    return "" if pd.isna(value) else str(value)


def display_status(value) -> str:
    value = display_value(value)
    return "ReviewMaybe" if value == "ScreenMaybe" else value


def derive_export_review_reason(row: pd.Series) -> tuple[str, str]:
    status = str(row.get("Status", "")).strip()
    best_gene = display_value(row.get("Best-Gene", "")).strip()
    product = str(row.get("Final-Product", row.get("Product-Consensus", ""))).strip().lower()
    if status == "ScreenMaybe" and best_gene == "HP" and product not in {"", "hypothetical protein", "hypothetical_protein"}:
        best_gene = ""

    if status == "Overlap":
        return "OVL-01", "Final loci still overlap after grouping and need manual resolution."
    if status == "Short-CDS":
        return "LEN-01", "Short CDS length outlier flagged for manual review."
    if status == "Long-CDS":
        return "LEN-02", "Long CDS length outlier flagged for manual review."
    if status == "Pseudogene" or bool(row.get("Pseudogene", False)):
        return "PSD-01", "Pseudogene flag present in the final output."
    if status == "ScreenMaybe":
        if not best_gene and ("putative" in product or "uncharacterized" in product or "unknown" in product):
            return "SCR-01", "Product-only review call with a generic or putative product label."
        if not best_gene:
            return "SCR-02", "Product-only review call without a trusted gene symbol."
        return "SCR-03", "Review call with weak or rescued naming that still needs checking."
    return "", ""


def derive_manual_check_reason(row: pd.Series) -> tuple[str, str]:
    check_group = str(row.get("check_group", "")).strip()
    status = str(row.get("Status", "")).strip()
    consensus = str(row.get("Consensus-Level", "")).strip()
    best_gene = display_value(row.get("Best-Gene", "")).strip()
    product = str(row.get("Final-Product", row.get("Product-Consensus", ""))).strip().lower()

    if check_group == "Overlap":
        return "OVL-01", "Final loci still overlap after grouping and need manual resolution."
    if check_group == "Length Outlier":
        if status == "Short-CDS":
            return "LEN-01", "Short CDS length outlier flagged for manual review."
        return "LEN-02", "Long CDS length outlier flagged for manual review."
    if check_group == "Local Conflict":
        if consensus == "Single-Tool":
            return "LOC-01", "Same-strand local conflict with only single-tool support."
        return "LOC-02", "Same-strand local conflict remains after grouping."
    if check_group == "ReviewMaybe":
        if not best_gene and ("putative" in product or "uncharacterized" in product or "unknown" in product):
            return "SCR-01", "Product-only ScreenMaybe call with a generic or putative product label."
        if not best_gene:
            return "SCR-02", "Product-only ScreenMaybe call without a trusted gene symbol."
        return "SCR-03", "ScreenMaybe row with weak or rescued naming that still needs checking."
    return "REV-00", "Manual review requested by dashboard prioritisation rules."


def write_qc_dashboard(df: pd.DataFrame, output_folder: str) -> str:
    df = blank_screenmaybe_hp_gene_display(df)
    out_path = os.path.join(output_folder, "QC-Dashboard.html")
    priority_screenmaybe_mask = (
        (df["Status"] == "ScreenMaybe") &
        (df["Confidence-Score"] <= 60.0)
    )
    length_outlier_mask = df["Status"].isin(["Short-CDS", "Long-CDS"])
    primary_manual_check_count = int(
        (df["Status"] == "Overlap").sum() +
        priority_screenmaybe_mask.sum() +
        length_outlier_mask.sum()
    )
    type_counts = df["Type"].value_counts().sort_index()
    status_counts = df["Status"].fillna("").replace("", "Empty").value_counts().sort_index()
    source_counts = (
        df["Source-of-Gene-Name"]
        .fillna("")
        .apply(lambda x: ",".join(sorted({p.strip() for p in str(x).split(",") if p.strip()})))
        .replace("", "Unknown")
        .value_counts()
        .head(15)
    )

    def series_table(title: str, series: pd.Series) -> str:
        rows = "".join(
            f"<tr><td>{label}</td><td>{int(value)}</td></tr>"
            for label, value in series.items()
        )
        return f"<section><h2>{title}</h2><table><tbody>{rows}</tbody></table></section>"

    summary = {
        "Total features": int(len(df)),
        "High confidence (>=85)": int((df["Confidence-Score"] >= 85).sum()),
        "Medium confidence (60-84)": int(((df["Confidence-Score"] >= 60) & (df["Confidence-Score"] < 85)).sum()),
        "Low confidence (<60)": int((df["Confidence-Score"] < 60).sum()),
        "Mean confidence": round(float(df["Confidence-Score"].mean()), 2) if len(df) else 0.0,
        "Worth reviewing annotations": primary_manual_check_count,
    }
    summary_html = "".join(f"<li><strong>{k}:</strong> {v}</li>" for k, v in summary.items())
    local_conflict_tags = set()
    for _, group in df.groupby(["Sequence_ID", "Direction", "Type"], dropna=False):
        group = group.sort_values(["Start", "End", "Locus_Tag"])
        active = []
        for _, row in group.iterrows():
            active = [prev for prev in active if int(prev["End"]) >= int(row["Start"])]
            for prev in active:
                if prev["Locus_Tag"] != row["Locus_Tag"]:
                    local_conflict_tags.add(str(prev["Locus_Tag"]))
                    local_conflict_tags.add(str(row["Locus_Tag"]))
            active.append(row)

    local_conflict_review_tags = set(
        df[
            df["Locus_Tag"].astype(str).isin(local_conflict_tags)
            & (df["Confidence-Score"] < LOCAL_CONFLICT_REVIEW_MAX_CONFIDENCE)
        ]["Locus_Tag"].astype(str)
    )

    grouped_counts = {
        "Overlap": int((df["Status"] == "Overlap").sum()),
        "ReviewMaybe (<=60)": int(priority_screenmaybe_mask.sum()),
        "Length outliers": int(length_outlier_mask.sum()),
        "Local Conflicts": int(len(local_conflict_review_tags)),
    }
    grouped_rows = "".join(
        f"<tr><td>{label}</td><td>{value}</td><td>{desc}</td></tr>"
        for label, value, desc in [
            ("Overlap", grouped_counts["Overlap"], "Rows explicitly flagged because another final locus still overlaps after grouping."),
            ("ReviewMaybe (<=60)", grouped_counts["ReviewMaybe (<=60)"], "Lower-confidence ReviewMaybe rows prioritised for review at confidence <= 60."),
            ("Length outliers", grouped_counts["Length outliers"], "Rows flagged as Short-CDS or Long-CDS structural outliers."),
            ("Local Conflicts", grouped_counts["Local Conflicts"], f"Rows in same-strand, same-type overlap neighborhoods with confidence below {LOCAL_CONFLICT_REVIEW_MAX_CONFIDENCE:.0f}, even if they are not explicitly marked as Overlap."),
        ]
    )

    manual_df = df[
        (df["Status"] == "Overlap") |
        priority_screenmaybe_mask |
        length_outlier_mask
    ].copy()
    if len(manual_df):
        manual_df["check_group"] = "Length Outlier"
        manual_df.loc[manual_df["Status"] == "Overlap", "check_group"] = "Overlap"
        manual_df.loc[priority_screenmaybe_mask.loc[manual_df.index], "check_group"] = "ReviewMaybe"
        manual_df.loc[length_outlier_mask.loc[manual_df.index], "check_group"] = "Length Outlier"
        reason_pairs = manual_df.apply(derive_manual_check_reason, axis=1, result_type="expand")
        manual_df["reason_id"] = reason_pairs[0]
        manual_df["reason_text"] = reason_pairs[1]
        manual_df["group_rank"] = manual_df["check_group"].map(
            {"Overlap": 0, "ReviewMaybe": 1, "Length Outlier": 2}
        ).fillna(9)
        manual_df = manual_df.sort_values(
            ["group_rank", "Confidence-Score", "Sequence_ID", "Start", "End"],
            ascending=[True, True, True, True, True],
        )
        manual_rows = "".join(
            (
                "<tr>"
                f"<td>{row.get('check_group', '')}</td>"
                f"<td>{row.get('Locus_Tag', '')}</td>"
                f"<td>{row.get('Sequence_ID', '')}</td>"
                f"<td>{row.get('Type', '')}</td>"
                f"<td>{row.get('Start', '')}</td>"
                f"<td>{row.get('End', '')}</td>"
                f"<td>{row.get('Direction', '')}</td>"
                f"<td>{display_value(row.get('Tool', ''))}</td>"
                f"<td>{display_value(row.get('Best-Gene', ''))}</td>"
                f"<td>{display_status(row.get('Status', ''))}</td>"
                f"<td>{display_value(row.get('reason_id', ''))}</td>"
                f"<td>{display_value(row.get('reason_text', ''))}</td>"
                f"<td>{display_value(row.get('Consensus-Level', ''))}</td>"
                f"<td>{display_value(row.get('Confidence-Score', ''))}</td>"
                f"<td>{display_value(row.get('Final-Product', row.get('Product-Consensus', '')))}</td>"
                "</tr>"
            )
            for _, row in manual_df.iterrows()
        )
    else:
        manual_rows = '<tr><td colspan="15">No manual-check rows.</td></tr>'

    reason_key_rows = "".join(
        f"<tr><td>{rid}</td><td>{desc}</td></tr>"
        for rid, desc in [
            ("OVL-01", "Final loci still overlap after grouping."),
            ("SCR-01", "Product-only ReviewMaybe call with a generic or putative product label."),
            ("SCR-02", "Product-only ReviewMaybe call without a trusted gene symbol."),
            ("SCR-03", "ReviewMaybe row with weak or rescued naming."),
            ("LEN-01", "Short CDS length outlier."),
            ("LEN-02", "Long CDS length outlier."),
        ]
    )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Bactowise QC Dashboard</title>
  <style>
    body {{ font-family: Helvetica, Arial, sans-serif; margin: 2rem; color: #182026; background: #f7f5ef; }}
    h1, h2 {{ color: #103c3a; }}
    p {{ max-width: 72rem; line-height: 1.55; }}
    section {{ margin: 1.2rem 0 1.5rem; }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ border-collapse: collapse; width: min(42rem, 100%); background: #fff; }}
    table.wide {{ width: 100%; min-width: 72rem; }}
    td, th {{ border: 1px solid #d7d2c8; padding: 0.45rem 0.6rem; text-align: left; vertical-align: top; }}
    ul {{ background: #fff; border: 1px solid #d7d2c8; padding: 1rem 1.4rem; width: min(42rem, 100%); }}
  </style>
</head>
<body>
  <h1>Bactowise QC Dashboard</h1>
  <p>This workflow is designed as a consensus harmonisation layer rather than a replacement annotation engine. It reads multiple annotation tracks, aligns them into a common schema, groups loci by shared coordinates, and keeps source provenance so every final call can still be traced back to the supporting tools.</p>
  <p><strong>Confidence score</strong> is a prioritisation metric for review, not a calibrated probability. It increases with broader tool support, cleaner coordinate agreement, and better gene or product agreement, and it is penalised for overlap or obvious structural warnings. A lower score therefore usually means “inspect this locus first”, not “this annotation is necessarily wrong”.</p>
  <p>The most useful points to watch are loci marked as overlap cases, lower-confidence <code>ReviewMaybe</code> rows, and length outliers. Local conflict neighborhoods are still counted in the summary, but they are not all expanded into the table below.</p>
  <section><h2>Summary</h2><ul>{summary_html}</ul></section>
  {series_table("Feature Types", type_counts)}
  {series_table("Statuses", status_counts)}
  {series_table("Source Support", source_counts)}
  <section>
    <h2>Review Summary</h2>
    <table><tbody>{grouped_rows}</tbody></table>
  </section>
  <section>
    <h2>Reason Key</h2>
    <table><tbody>{reason_key_rows}</tbody></table>
  </section>
  <section>
    <h2>Worth Reviewing Annotations</h2>
    <p>Only the top low-confidence <code>ReviewMaybe</code> rows are shown here, together with overlap and length-outlier cases. Full details remain in <code>Master_Table_Annotation.xlsx</code>.</p>
    <div class="table-wrap">
      <table class="wide">
        <thead>
          <tr>
            <th>Check-Group</th>
            <th>Locus_Tag</th>
            <th>Sequence_ID</th>
            <th>Type</th>
            <th>Start</th>
            <th>End</th>
            <th>Direction</th>
            <th>Tool</th>
            <th>Best-Gene</th>
            <th>Status</th>
            <th>Reason-ID</th>
            <th>Reason</th>
            <th>Consensus-Level</th>
            <th>Confidence-Score</th>
            <th>Final-Product</th>
          </tr>
        </thead>
        <tbody>{manual_rows}</tbody>
      </table>
    </div>
  </section>
</body>
</html>
"""
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html_doc)
    logger.info(f"  QC dashboard → {out_path}")
    return out_path


# ==============================================================================
# OUTPUT: SUMMARY REPORT (plain text)
# ==============================================================================
def generate_summary(df: pd.DataFrame, output_folder: str):
    total = len(df)
    if total == 0:
        logger.warning("Summary skipped — final DataFrame is empty")
        return

    cds_df    = df[df["Type"] == "CDS"]
    total_cds = len(cds_df)

    def pct_all(n):
        return f"{n:>6}  ({100 * n / total:5.1f}% of all features)"

    def pct_cds(n):
        if total_cds == 0:
            return f"{n:>6}"
        return f"{n:>6}  ({100 * n / total_cds:5.1f}% of CDS)"

    # Consensus levels actually present in this run
    consensus_levels = sorted(
        [l for l in df["Consensus-Level"].dropna().unique() if str(l).startswith("Consensus_")],
        key=lambda x: int(x.split("_")[1].split("/")[0]),
        reverse=True,
    )

    lines = [
        "=" * 60,
        "  UNIFIED ANNOTATION PIPELINE — SUMMARY REPORT",
        f"  Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 60,
        "",
        "TOTALS",
        f"  Total features  : {total}",
        f"  Total CDS       : {total_cds}",
        "",
        "CONSENSUS  (all features)",
    ]

    for level in consensus_levels:
        n = (df["Consensus-Level"] == level).sum()
        lines.append(f"  {level:<22}: {pct_all(n)}")
    n_single = (df["Consensus-Level"] == "Single-Tool").sum()
    lines.append(f"  {'Single-Tool':<22}: {pct_all(n_single)}")

    lines += [
        "",
        "GENE QUALITY  (CDS only)",
        f"  Hypothetical proteins (HP)  : {pct_cds((cds_df['Best-Gene'] == 'HP').sum())}",
        f"  HP rescued via product      : {pct_cds((cds_df['Status'] == 'HP-ProductRescued').sum())}",
        f"  Pseudogenes                 : {pct_cds((cds_df['Status'] == 'Pseudogene').sum())}",
        "",
        "STATUS  (all features)",
        f"  Good                        : {pct_all((df['Status'] == 'Good').sum())}",
        f"  ScreenMaybe                 : {pct_all((df['Status'] == 'ScreenMaybe').sum())}",
        f"  Overlap flagged             : {pct_all((df['Status'] == 'Overlap').sum())}",
        f"  Short CDS (<90 bp)          : {pct_all((df['Status'] == 'Short-CDS').sum())}",
        f"  Long CDS (>15 kb)           : {pct_all((df['Status'] == 'Long-CDS').sum())}",
        "",
        "SEQUENCE OUTPUT",
        f"  FAA/FNA eligible (all CDS, excl. pseudogene/short)",
        f"    : {pct_cds(((cds_df['Pseudogene'] != True) & (cds_df['Status'] != 'Short-CDS')).sum())}",
        "",
        "FEATURE TYPES  (all features)",
    ]

    for ftype, count in df["Type"].value_counts().items():
        lines.append(f"  {ftype:<30}: {pct_all(count)}")

    lines += ["", "=" * 60]

    out_path = os.path.join(output_folder, "summary_report.txt")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    logger.info(f"  Summary report → {out_path}")

    # Echo to log
    logger.info("\n" + "\n".join(lines))


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================
def validate_inputs(input_folder: str):
    if not os.path.isdir(input_folder):
        raise ValueError(f"Input folder not found: {input_folder}")

    gff_files = [f for f in os.listdir(input_folder)
                 if f.endswith(".gff") or f.endswith(".gff3")]
    if not gff_files:
        raise ValueError(f"No GFF/GFF3 files found in: {input_folder}")

    detected = {detect_tool(f) for f in gff_files}
    logger.info(f"  GFF files   : {', '.join(sorted(gff_files))}")
    logger.info(f"  Detected tools: {', '.join(sorted(detected))}")

    if len(gff_files) < 2:
        logger.warning("Only 1 GFF file — all features will be Single-Tool")
    if "unknown" in detected:
        logger.warning(
            "Could not detect tool from filename(s). "
            "Rename files to include 'bakta' or 'prokka' for best results."
        )

    fasta = find_fasta(input_folder)
    logger.info(f"  FASTA file  : {os.path.basename(fasta)}")
    return fasta


# ==============================================================================
# CLEANUP TEMP FILES
# ==============================================================================
def cleanup_temp_files(output_folder: str):
    for f in os.listdir(output_folder):
        if f.startswith("_tmp_") and f.endswith(".xlsx"):
            os.remove(os.path.join(output_folder, f))


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================
def run_annotation_pipeline(input_folder: str, output_folder: str,
                             prefix: str = "GENE") -> str | None:
    setup_logging(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    logger.info("\n" + "=" * 70)
    logger.info("UNIFIED ANNOTATION PIPELINE v3")
    logger.info(f"  Input    : {input_folder}")
    logger.info(f"  Output   : {output_folder}")
    logger.info(f"  Prefix   : {prefix}")
    logger.info("=" * 70)

    try:
        # ── Validation ────────────────────────────────────────────────────────
        logger.info("\n[1] Validating inputs...")
        fasta_path = validate_inputs(input_folder)

        # ── Load genome ───────────────────────────────────────────────────────
        logger.info("\n[2] Loading genome FASTA...")
        genome = parse_fasta(fasta_path)
        logger.info(f"  {len(genome)} contig(s) loaded, "
                    f"total {sum(len(s) for s in genome.values()):,} bp")

        # ── GFF parsing ───────────────────────────────────────────────────────
        logger.info("\n[3] Parsing GFF files...")
        stems = process_gff_files_in_folder(input_folder, output_folder)

        # ── Merge ─────────────────────────────────────────────────────────────
        logger.info("\n[4] Merging tool annotations...")
        df = merge_excel(output_folder, stems)

        logger.info("\n[4b] Harmonizing sequence IDs across tools and FASTA...")
        df, genome = harmonize_sequence_ids(df, genome)

        # ── Filter ────────────────────────────────────────────────────────────
        logger.info("\n[5] Filtering feature types...")
        df = filter_feature_types(df)

        # ── Clean gene names ──────────────────────────────────────────────────
        logger.info("\n[6] Cleaning gene names...")
        df["Gene"] = df["Gene"].apply(clean_gene_name)

        # ── Grouping ──────────────────────────────────────────────────────────
        logger.info("\n[7] Grouping V1 (identical coordinates)...")
        df = grouping_v1(df)

        logger.info("\n[8] Grouping V2 (boundary resolution)...")
        df = grouping_v2(df)

        logger.info("\n[9] Grouping V2b (non-CDS near-identical span merge)...")
        df = grouping_v2_non_cds_near_match(df)

        logger.info("\n[10] Grouping V3 (overlap detection)...")
        df = grouping_v3(df)

        # ── Best gene ─────────────────────────────────────────────────────────
        logger.info("\n[11] Selecting best gene per locus...")
        df["Best-Gene"] = df.apply(find_best_gene, axis=1)

        # ── HP + pseudogene flags ────────────────────────────────────────────
        logger.info("\n[11] Flagging hypothetical proteins and pseudogenes...")
        df = flag_hp_and_pseudogenes(df)

        # ── HP rescue ────────────────────────────────────────────────────────
        logger.info("\n[12] HP rescue via product strings...")
        df["Product-Consensus"] = ""
        df = rescue_hp_by_product(df)

        # ── Single gene + status ─────────────────────────────────────────────
        logger.info("\n[13] Ensuring single gene and status per locus...")
        df = df.apply(update_columns, axis=1)

        # ── Gene synonyms ────────────────────────────────────────────────────
        logger.info("\n[14] Extracting gene synonyms...")
        df["Gene-Putative-Synonyms"] = df.apply(extract_gene_synonyms, axis=1)

        # ── Consensus scoring ────────────────────────────────────────────────
        logger.info("\n[15] Consensus scoring...")
        df = consensus_scoring(df)

        # ── Quality flags ────────────────────────────────────────────────────
        logger.info("\n[16] Length outlier and strand bias checks...")
        df = flag_length_outliers(df)
        check_strand_bias(df)

        # ── Locus tags ───────────────────────────────────────────────────────
        logger.info("\n[17] Assigning locus tags...")
        df = assign_locus_tags(df, prefix)

        # ── Functional categories ────────────────────────────────────────────
        logger.info("\n[18] Assigning functional categories...")
        df["Functional-Category"] = df["Best-Gene"].apply(assign_functional_category)
        df["Final-Product"] = df.apply(best_final_product, axis=1)
        df["Final-Start_End_Direction"] = (
            df["Start"].astype(int).astype(str)
            + "_"
            + df["End"].astype(int).astype(str)
            + "_"
            + df["Direction"].astype(str)
        )

        # ── Write outputs ────────────────────────────────────────────────────
        logger.info("\n[19] Writing output files...")
        master = write_master_table(df, output_folder)
        write_gff(df, output_folder, prefix)
        write_gbk(df, genome, output_folder, prefix)
        write_faa(df, genome, output_folder, prefix)
        write_fna(df, genome, output_folder, prefix)
        write_manual_queue(df, output_folder)
        write_qc_dashboard(df, output_folder)

        # ── Summary ──────────────────────────────────────────────────────────
        logger.info("\n[20] Generating summary report...")
        generate_summary(df, output_folder)

        # ── Cleanup ───────────────────────────────────────────────────────────
        cleanup_temp_files(output_folder)

        logger.info("\n" + "=" * 70)
        logger.info("PIPELINE COMPLETE")
        logger.info(f"  Master_Table_Annotation.xlsx")
        logger.info(f"  Manual-Curation-Queue.csv")
        logger.info(f"  Manual-Curation-Queue_full.csv")
        logger.info(f"  QC-Dashboard.html")
        logger.info(f"  summary_report.txt")
        logger.info(f"  {prefix}.gff3  (Geneious / SnapGene / Artemis / IGV compatible)")
        logger.info(f"  {prefix}.gbk   (Geneious / SnapGene / Addgene compatible)")
        logger.info(f"  {prefix}.faa   (protein sequences, high-confidence CDS)")
        logger.info(f"  {prefix}.fna   (nucleotide sequences, high-confidence CDS)")
        logger.info(f"  pipeline.log")
        logger.info("=" * 70 + "\n")

        return master

    except Exception as e:
        logger.error(f"\n✗ PIPELINE FAILED: {e}", exc_info=True)
        return None


# ==============================================================================
# CLI
# ==============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Unified Bacterial Annotation Pipeline v3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples:
          python unified_annotation_pipeline_v3.py \\
              --input  /data/gff_files \\
              --output /data/results \\
              --prefix ECOLI

        Input folder must contain:
          *.gff / *.gff3   (Bakta and/or Prokka output — filenames must include tool name)
          *.fasta / *.fna / *.fa   (the assembly used for annotation)

        Output files:
          Master_Table_Annotation.xlsx
          summary_report.txt
          <prefix>.gff3   GFF3  — Geneious, SnapGene, Artemis, IGV
          <prefix>.gbk    GenBank — Geneious, SnapGene, Addgene
          <prefix>.faa    Protein FASTA — high-confidence CDS only
          <prefix>.fna    Nucleotide FASTA — high-confidence CDS only
          pipeline.log
        """),
    )

    parser.add_argument("--input",  required=True,
                        help="Folder with GFF + FASTA files")
    parser.add_argument("--output", required=True,
                        help="Output folder")
    parser.add_argument("--prefix", default="GENE",
                        help="Locus tag prefix  [default: GENE]")

    args   = parser.parse_args()
    result = run_annotation_pipeline(args.input, args.output, args.prefix)
    exit(0 if result else 1)
