#!/usr/bin/env python3
"""
Step 3 only: annotate ITS1/ITS2 representative sequences against marker-specific
reference resources and generate candidate lineage tables that match the frozen
Step 3 schema used by the downstream pipeline.

Fixes included:
1) candidate lineage tables use `reads` instead of `abundance`
2) candidate lineage tables explicitly retain:
   sample_id, feature_id, sequence, marker, taxon_rank, taxon_name, reads, confidence
3) ASV sample columns are validated and treated as frozen sample_id columns
4) unresolved / ambiguous / no-hit features are preserved for Step 4
5) resolution-limited, low-confidence-review, and non-specific-label flags are added
6) reference provenance is frozen into annotation_confidence.tsv / unassigned_features.tsv
7) duplicate feature_id rows in ASV tables are auto-collapsed by summing reads
8) transposed ASV tables (sample x feature) are auto-detected and transposed
9) non-feature metadata columns in transposed ASV tables (e.g. replicate) are auto-dropped
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


REQUIRED_INPUTS = {
    "ITS1_rep_seqs.fasta": "ITS1",
    "ITS2_rep_seqs.fasta": "ITS2",
    "ITS1_asv_table.tsv": "ITS1",
    "ITS2_asv_table.tsv": "ITS2",
    "dada2_stats_ITS1.tsv": "ITS1",
    "dada2_stats_ITS2.tsv": "ITS2",
}

REQUIRED_OUTPUTS = [
    "ITS1_candidate_lineages.tsv",
    "ITS2_candidate_lineages.tsv",
    "annotation_confidence.tsv",
    "annotation_rank_summary.tsv",
    "unassigned_features.tsv",
]

FROZEN_REFERENCE = {
    "ITS1": {
        "resource_name": "PLANiTS",
        "download_date": "2026-03-28",
        "frozen_package": "PLANiTS_29-03-2020.zip",
        "record_id": "NA",
    },
    "ITS2": {
        "resource_name": "ITS2 Global curated reference",
        "download_date": "2026-03-28",
        "frozen_package": "its2.global.2023-01-17.curated.tax.mc.add.fa",
        "record_id": "7968519",
    },
}

METHOD_STACK = {
    "name": "NCBI BLAST+ nucleotide alignment + ambiguity-aware rank fallback",
    "engine": "blastn / makeblastdb",
    "confidence_definition": "sortable BLAST bitscore, with pident/evalue/qcovs trace retained",
}

RESERVED_NON_SAMPLE_COLUMNS = {
    "taxonomy",
    "consensus",
    "confidence",
    "sequence",
    "taxon_rank",
    "taxon_name",
    "status",
    "status_reason",
    "marker",
}


class Step3Error(RuntimeError):
    pass


@dataclass
class ReferenceLayout:
    marker: str
    source_path: Path
    source_type: str  # zip / fasta
    fasta_path: Path
    taxonomy_table_path: Optional[Path]
    taxonomy_mode: str  # table / header_only
    frozen_match: bool
    frozen_expected_name: str
    sha256: str


@dataclass
class AnnotationConfig:
    ambiguity_policy: str
    rank_resolution_policy: str
    allow_unassigned: bool
    sample_id_regex: Optional[str]
    review_pident_below: Optional[float]
    review_qcov_below: Optional[float]
    review_evalue_above: Optional[float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step 3 reference annotation only")
    parser.add_argument(
        "--result-dir",
        default="/mnt/d/其他-小论文/数据库/zonode/陆宏龙/2/优化方案2/结果",
        help="Step 2 result directory (WSL path default)",
    )
    parser.add_argument(
        "--windows-result-dir",
        default=r"D:\其他-小论文\数据库\zonode\陆宏龙\2\优化方案2\结果",
        help="Windows mapping path printed for traceability only",
    )
    parser.add_argument(
        "--reference-root",
        default="/mnt/d/其他-小论文/数据库/zonode/陆宏龙/2/优化方案2/结果",
        help="Root folder containing frozen reference resources",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: same as result-dir)",
    )
    parser.add_argument(
        "--ambiguity-policy",
        default="bitscore_tie_resolve_if_same_taxon_else_unassigned",
        choices=["bitscore_tie_resolve_if_same_taxon_else_unassigned"],
        help="Tie handling policy for equally best BLAST hits",
    )
    parser.add_argument(
        "--rank-resolution-policy",
        default="species_then_genus_then_higher",
        choices=["species_then_genus_then_higher"],
        help="Rank fallback policy",
    )
    parser.add_argument(
        "--allow-unassigned",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep unresolved features as unassigned (recommended by Step 3 constraints)",
    )
    parser.add_argument(
        "--sample-id-regex",
        default=None,
        help="Optional regex to validate ASV sample columns as frozen sample_id values",
    )
    parser.add_argument(
        "--review-pident-below",
        type=float,
        default=None,
        help="Optional review-only threshold; does not delete rows",
    )
    parser.add_argument(
        "--review-qcov-below",
        type=float,
        default=None,
        help="Optional review-only threshold; does not delete rows",
    )
    parser.add_argument(
        "--review-evalue-above",
        type=float,
        default=None,
        help="Optional review-only threshold; does not delete rows",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise Step3Error(f"Required dependency not found in PATH: {name}")


def _read_fasta_records(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    current_id: Optional[str] = None
    seq_parts: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(seq_parts).strip()
                    if not seq:
                        raise Step3Error(f"FASTA record without sequence in {path}: {current_id}")
                    records[current_id] = seq

                header = line[1:].strip()
                if not header:
                    raise Step3Error(f"Empty FASTA header in {path}")
                current_id = header.split()[0]
                if current_id in records:
                    raise Step3Error(f"Duplicate FASTA feature_id in {path}: {current_id}")
                seq_parts = []
            else:
                if current_id is None:
                    raise Step3Error(f"Invalid FASTA: sequence before header in {path}")
                if not re.fullmatch(r"[A-Za-z*.\-]+", line):
                    raise Step3Error(f"Invalid FASTA characters in {path}: {line[:30]}...")
                seq_parts.append(line)

    if current_id is not None:
        seq = "".join(seq_parts).strip()
        if not seq:
            raise Step3Error(f"Last FASTA record without sequence in {path}: {current_id}")
        records[current_id] = seq

    if not records:
        raise Step3Error(f"No FASTA records found in {path}")

    return records


def _read_fasta_headers(path: Path, max_records: int = 1000) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                header = raw[1:].strip()
                if header:
                    out.append((header.split()[0], header))
                    if len(out) >= max_records:
                        break
    return out


def _detect_feature_column(df: pd.DataFrame) -> str:
    candidates = [
        c for c in df.columns
        if str(c).strip().lower() in {"feature_id", "feature", "asv", "asv_id", "seq_id", "id"}
    ]
    if candidates:
        return str(candidates[0])

    first_col = str(df.columns[0])
    if first_col.startswith("Unnamed"):
        df.rename(columns={df.columns[0]: "feature_id"}, inplace=True)
        return "feature_id"

    return first_col


def _validate_sample_columns(sample_cols: List[str], sample_id_regex: Optional[str]) -> None:
    if not sample_cols:
        raise Step3Error("No sample_id columns found in ASV table")

    lowered = [c.strip().lower() for c in sample_cols]
    if len(lowered) != len(set(lowered)):
        raise Step3Error("Duplicate sample_id columns detected in ASV table")

    for col in sample_cols:
        c = str(col).strip()
        if not c:
            raise Step3Error("Blank sample_id column detected in ASV table")
        if c.lower() in RESERVED_NON_SAMPLE_COLUMNS:
            raise Step3Error(f"Unexpected non-sample metadata column in ASV table: {c}")
        if sample_id_regex and re.fullmatch(sample_id_regex, c) is None:
            raise Step3Error(f"ASV column failed sample_id regex validation: {c}")


def _reorient_asv_table_if_needed(
    marker: str,
    asv_df: pd.DataFrame,
    seq_map: Dict[str, str],
    sample_id_regex: Optional[str],
) -> pd.DataFrame:
    """
    Detect whether ASV table is transposed:
      expected orientation: feature_id in rows, sample_id in columns
      transposed orientation: sample_id in rows, feature_id in columns

    If transposed, automatically convert it into the expected orientation.
    Also collapse duplicated sample rows before transpose by summing reads.
    Also auto-drop non-feature metadata columns (e.g. replicate) before transpose.
    """
    if asv_df.empty or asv_df.shape[1] < 2:
        raise Step3Error(f"ASV table is empty or malformed for {marker}")

    first_col_values = asv_df.iloc[:, 0].astype(str).str.strip().tolist()
    remaining_cols = [str(c).strip() for c in asv_df.columns[1:]]

    row_feature_hits = sum(v in seq_map for v in first_col_values)
    col_feature_hits = sum(c in seq_map for c in remaining_cols)

    if col_feature_hits > row_feature_hits:
        print(
            f"[Step3] WARN: {marker} ASV table appears transposed "
            f"(samples in rows, features in columns); auto-transposing.",
            file=sys.stderr,
        )

        tmp = asv_df.copy()
        tmp.columns = [str(c).strip() for c in tmp.columns]
        row_id_col = str(tmp.columns[0])
        tmp[row_id_col] = tmp[row_id_col].astype(str).str.strip()

        sample_ids = tmp[row_id_col].tolist()
        if any(not s for s in sample_ids):
            raise Step3Error(f"Blank sample_id rows detected in transposed ASV table for {marker}")

        if sample_id_regex:
            bad_ids = [s for s in sample_ids if re.fullmatch(sample_id_regex, s) is None]
            if bad_ids:
                raise Step3Error(
                    f"Some sample_id rows failed regex validation in transposed ASV table for {marker}: {bad_ids[:5]}"
                )

        candidate_cols = [c for c in tmp.columns if c != row_id_col]

        # Keep only real feature columns that are present in rep_seqs FASTA
        feature_cols = [c for c in candidate_cols if c in seq_map]

        # Drop metadata / non-feature columns automatically
        dropped_metadata_cols = [c for c in candidate_cols if c not in seq_map]
        if dropped_metadata_cols:
            print(
                f"[Step3] WARN: {marker} transposed ASV table contains non-feature metadata columns; "
                f"dropping before transpose: {dropped_metadata_cols[:10]}",
                file=sys.stderr,
            )

        if not feature_cols:
            raise Step3Error(f"No usable feature columns detected in transposed ASV table for {marker}")

        for col in feature_cols:
            tmp[col] = pd.to_numeric(tmp[col], errors="coerce").fillna(0)
            if (tmp[col] < 0).any():
                raise Step3Error(
                    f"Negative reads detected in transposed ASV table for {marker}, feature column {col}"
                )

        # Keep only sample_id column + true feature columns
        tmp = tmp[[row_id_col] + feature_cols].copy()

        if tmp[row_id_col].duplicated().any():
            dup_rows = int(tmp[row_id_col].duplicated(keep=False).sum())
            dup_ids = int(tmp.loc[tmp[row_id_col].duplicated(keep=False), row_id_col].nunique())
            print(
                f"[Step3] WARN: {marker} transposed ASV table contains duplicate sample rows "
                f"({dup_rows} rows across {dup_ids} duplicated sample_id values); "
                f"collapsing by summing reads before transpose.",
                file=sys.stderr,
            )
            tmp = (
                tmp.groupby(row_id_col, as_index=False)[feature_cols]
                .sum()
                .reset_index(drop=True)
            )

        tmp = tmp.set_index(row_id_col)
        tmp.index.name = "sample_id"
        tmp = tmp.T.reset_index().rename(columns={"index": "feature_id"})
        tmp.columns = [str(c).strip() for c in tmp.columns]

        return tmp

    return asv_df


def load_step2_outputs(result_dir: Path) -> Dict[str, Path]:
    if not result_dir.exists() or not result_dir.is_dir():
        raise Step3Error(f"Result directory does not exist or is not readable: {result_dir}")

    missing: List[str] = []
    file_map: Dict[str, Path] = {}

    for fname in REQUIRED_INPUTS:
        p = result_dir / fname
        if p.exists():
            file_map[fname] = p
        else:
            missing.append(fname)

    if missing:
        raise Step3Error(f"Missing Step 2 input files: {', '.join(missing)}")

    return file_map


def validate_step2_input_schema(
    file_map: Dict[str, Path],
    sample_id_regex: Optional[str],
) -> Dict[str, object]:
    out: Dict[str, object] = {}

    for marker in ("ITS1", "ITS2"):
        rep_fasta = file_map[f"{marker}_rep_seqs.fasta"]
        seq_map = _read_fasta_records(rep_fasta)
        out[f"{marker}_seq_map"] = seq_map

        asv_path = file_map[f"{marker}_asv_table.tsv"]
        asv_df = pd.read_csv(asv_path, sep="\t", dtype=str)
        if asv_df.empty or asv_df.shape[1] < 2:
            raise Step3Error(f"ASV table is empty or malformed: {asv_path}")

        asv_df = _reorient_asv_table_if_needed(
            marker=marker,
            asv_df=asv_df,
            seq_map=seq_map,
            sample_id_regex=sample_id_regex,
        )

        feature_col = _detect_feature_column(asv_df)
        asv_df[feature_col] = asv_df[feature_col].astype(str).str.strip()

        if asv_df[feature_col].eq("").any():
            raise Step3Error(f"Blank feature_id found in ASV table: {asv_path}")

        sample_cols = [str(c) for c in asv_df.columns if str(c) != feature_col]
        _validate_sample_columns(sample_cols, sample_id_regex)

        if asv_df[feature_col].duplicated().any():
            dup_rows = int(asv_df[feature_col].duplicated(keep=False).sum())
            dup_features = int(
                asv_df.loc[asv_df[feature_col].duplicated(keep=False), feature_col].nunique()
            )
            print(
                f"[Step3] WARN: {marker} ASV table contains duplicate feature_id rows "
                f"({dup_rows} rows across {dup_features} duplicated feature_id values); "
                f"collapsing by summing reads across sample_id columns.",
                file=sys.stderr,
            )

            work = asv_df[[feature_col] + sample_cols].copy()

            for col in sample_cols:
                work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0)
                if (work[col] < 0).any():
                    raise Step3Error(
                        f"Negative reads detected before duplicate-collapse in {asv_path}, column {col}"
                    )

            work = (
                work.groupby(feature_col, as_index=False)[sample_cols]
                .sum()
                .reset_index(drop=True)
            )

            for col in sample_cols:
                work[col] = work[col].round().astype(int).astype(str)

            asv_df = work.copy()

        missing_in_fasta = set(asv_df[feature_col]) - set(seq_map.keys())
        if missing_in_fasta:
            example = sorted(list(missing_in_fasta))[:5]
            raise Step3Error(
                f"ASV table contains feature_id not present in rep_seqs FASTA for {marker}: {example}"
            )

        stats_path = file_map[f"dada2_stats_{marker}.tsv"]
        stats_df = pd.read_csv(stats_path, sep="\t", dtype=str)
        if stats_df.empty:
            raise Step3Error(f"DADA2 stats is empty: {stats_path}")

        out[f"{marker}_asv_df"] = asv_df
        out[f"{marker}_feature_col"] = feature_col
        out[f"{marker}_sample_cols"] = sample_cols
        out[f"{marker}_stats_df"] = stats_df

    return out


def inspect_reference_layout(marker: str, source_path: Path, work_dir: Path) -> ReferenceLayout:
    if not source_path.exists() or not source_path.is_file():
        raise Step3Error(f"Reference source missing for {marker}: {source_path}")

    extracted_dir = work_dir / f"reference_{marker}"
    extracted_dir.mkdir(parents=True, exist_ok=True)

    if source_path.suffix.lower() == ".zip":
        source_type = "zip"
        with zipfile.ZipFile(source_path, "r") as zf:
            zf.extractall(extracted_dir)
    else:
        source_type = "fasta"
        shutil.copy2(source_path, extracted_dir / source_path.name)

    all_files = [p for p in extracted_dir.rglob("*") if p.is_file()]
    fasta_candidates = [p for p in all_files if p.suffix.lower() in {".fa", ".fasta", ".fna", ".fas"}]
    if not fasta_candidates:
        raise Step3Error(f"No FASTA found in {marker} reference resource: {source_path}")

    fasta_path = max(fasta_candidates, key=lambda p: p.stat().st_size)

    taxonomy_table_path: Optional[Path] = None
    for p in all_files:
        name = p.name.lower()
        if p.suffix.lower() in {".tsv", ".txt", ".csv"} and any(k in name for k in ("tax", "taxonomy", "lineage", "classif")):
            taxonomy_table_path = p
            break

    headers = _read_fasta_headers(fasta_path, max_records=500)
    header_has_taxonomy = any(
        any(tok in header.lower() for tok in ("k__", "p__", "c__", "o__", "f__", "g__", "s__", "taxonomy", "lineage", ";"))
        for _, header in headers
    )

    if taxonomy_table_path is not None:
        taxonomy_mode = "table"
    elif header_has_taxonomy:
        taxonomy_mode = "header_only"
    else:
        raise Step3Error(
            f"Cannot confirm taxonomy layout for {marker}; neither parseable taxonomy table nor taxonomy-like FASTA headers were found"
        )

    frozen_expected = FROZEN_REFERENCE[marker]["frozen_package"]
    return ReferenceLayout(
        marker=marker,
        source_path=source_path,
        source_type=source_type,
        fasta_path=fasta_path,
        taxonomy_table_path=taxonomy_table_path,
        taxonomy_mode=taxonomy_mode,
        frozen_match=(source_path.name == frozen_expected),
        frozen_expected_name=frozen_expected,
        sha256=_sha256(source_path),
    )


def locate_and_validate_reference_resources(reference_root: Path, work_dir: Path) -> Dict[str, ReferenceLayout]:
    if not reference_root.exists() or not reference_root.is_dir():
        raise Step3Error(f"Reference root missing: {reference_root}")

    planits_name = FROZEN_REFERENCE["ITS1"]["frozen_package"]
    its2_name = FROZEN_REFERENCE["ITS2"]["frozen_package"]

    planits_candidates = list(reference_root.rglob(planits_name))
    its2_candidates = list(reference_root.rglob(its2_name))

    if not planits_candidates:
        raise Step3Error(f"PLANiTS frozen package not found: {planits_name}")
    if not its2_candidates:
        raise Step3Error(f"ITS2 frozen FASTA not found: {its2_name}")

    layout_its1 = inspect_reference_layout("ITS1", planits_candidates[0], work_dir)
    layout_its2 = inspect_reference_layout("ITS2", its2_candidates[0], work_dir)

    if not layout_its1.frozen_match:
        raise Step3Error(
            f"ITS1 frozen version mismatch; expected {layout_its1.frozen_expected_name}, found {layout_its1.source_path.name}"
        )
    if not layout_its2.frozen_match:
        raise Step3Error(
            f"ITS2 frozen version mismatch; expected {layout_its2.frozen_expected_name}, found {layout_its2.source_path.name}"
        )
    if layout_its1.source_path.resolve() == layout_its2.source_path.resolve():
        raise Step3Error("ITS1 and ITS2 resolved to the same reference resource; forbidden by Step 3")

    return {"ITS1": layout_its1, "ITS2": layout_its2}


def freeze_reference_provenance(layouts: Dict[str, ReferenceLayout]) -> pd.DataFrame:
    rows = []
    for marker, layout in layouts.items():
        frozen = FROZEN_REFERENCE[marker]
        rows.append(
            {
                "marker": marker,
                "reference_name": frozen["resource_name"],
                "download_date_frozen": frozen["download_date"],
                "frozen_package": frozen["frozen_package"],
                "record_id": frozen["record_id"],
                "source_path": str(layout.source_path),
                "source_type": layout.source_type,
                "fasta_path": str(layout.fasta_path),
                "taxonomy_mode": layout.taxonomy_mode,
                "taxonomy_table_path": "" if layout.taxonomy_table_path is None else str(layout.taxonomy_table_path),
                "reference_sha256": layout.sha256,
            }
        )
    return pd.DataFrame(rows)


def build_reference_index(marker: str, layout: ReferenceLayout, index_dir: Path) -> Path:
    _require_tool("makeblastdb")

    marker_dir = index_dir / marker
    marker_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = marker_dir / "reference_db"

    cmd = [
        "makeblastdb",
        "-in", str(layout.fasta_path),
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", str(db_prefix),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise Step3Error(f"makeblastdb failed for {marker}: {proc.stderr.strip()}")

    expected = [db_prefix.with_suffix(s) for s in (".nhr", ".nin", ".nsq")]
    if not all(p.exists() for p in expected):
        raise Step3Error(f"BLAST index not fully created for {marker}")

    return db_prefix


def _load_taxonomy_table(layout: ReferenceLayout) -> Dict[str, str]:
    if layout.taxonomy_table_path is None:
        return {}

    path = layout.taxonomy_table_path
    try:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    except Exception as e:
        raise Step3Error(f"Failed to parse taxonomy table for {layout.marker}: {path} ({e})") from e

    if df.empty:
        raise Step3Error(f"Empty taxonomy table for {layout.marker}: {path}")

    cols_lower = {str(c).strip().lower(): c for c in df.columns}

    id_col = None
    for k in ("id", "seqid", "sequence_id", "accession", "feature_id", "name"):
        if k in cols_lower:
            id_col = cols_lower[k]
            break

    tax_col = None
    for k in ("taxonomy", "lineage", "taxon", "classification"):
        if k in cols_lower:
            tax_col = cols_lower[k]
            break

    if id_col is None or tax_col is None:
        raise Step3Error(
            f"Could not locate taxonomy id/taxonomy columns in table for {layout.marker}: {path}"
        )

    tax_map: Dict[str, str] = {}
    tmp = df[[id_col, tax_col]].dropna().copy()
    if tmp.empty:
        raise Step3Error(f"Taxonomy table has no usable rows for {layout.marker}: {path}")

    for _, row in tmp.iterrows():
        seqid = str(row[id_col]).strip()
        lineage = str(row[tax_col]).strip()
        if seqid and lineage:
            tax_map[seqid] = lineage

    if not tax_map:
        raise Step3Error(f"No taxonomy mappings were recovered for {layout.marker}: {path}")

    return tax_map


def _parse_taxonomy_from_text(text: str) -> Dict[str, str]:
    out = {r: "" for r in ("kingdom", "phylum", "class", "order", "family", "genus", "species")}
    text = (text or "").strip()
    if not text:
        return out

    rank_map = {
        "k": "kingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    for token in re.split(r"[;|]", text):
        token = token.strip()
        m = re.match(r"^([kpcofgs])__\s*(.+)$", token, flags=re.IGNORECASE)
        if m:
            out[rank_map[m.group(1).lower()]] = m.group(2).strip()

    if not out["species"]:
        m = re.search(r"([A-Z][A-Za-z._-]+\s+[a-z][A-Za-z._-]+)", text)
        if m:
            species = m.group(1).strip()
            out["species"] = species
            out["genus"] = species.split()[0]

    return out


def resolve_taxon_rank(tax: Dict[str, str], rank_resolution_policy: str) -> Dict[str, str]:
    if rank_resolution_policy != "species_then_genus_then_higher":
        raise Step3Error(f"Unsupported rank resolution policy: {rank_resolution_policy}")

    species = (tax.get("species") or "").strip()
    genus = (tax.get("genus") or "").strip()

    bad_words = r"\b(sp\.?|uncultured|unidentified|environmental|metagenome|unknown)\b"

    if species and re.search(bad_words, species.lower()):
        species = ""
    if genus and re.search(bad_words, genus.lower()):
        genus = ""

    if species:
        return {
            "taxon_rank": "species",
            "taxon_name": species,
            "status_reason": "species_resolved",
            "resolution_flag": "species_resolved",
        }

    if genus:
        return {
            "taxon_rank": "genus",
            "taxon_name": genus,
            "status_reason": "downgraded_to_genus",
            "resolution_flag": "resolution_limited_genus",
        }

    for rank in ("family", "order", "class", "phylum", "kingdom"):
        value = (tax.get(rank) or "").strip()
        if value and not re.search(bad_words, value.lower()):
            return {
                "taxon_rank": rank,
                "taxon_name": value,
                "status_reason": f"downgraded_to_{rank}",
                "resolution_flag": "resolution_limited_higher_rank",
            }

    return {
        "taxon_rank": "unassigned",
        "taxon_name": "unassigned",
        "status_reason": "no_resolvable_taxonomy",
        "resolution_flag": "unassigned",
    }


def _annotation_scope_flag(raw_tax_text: str) -> str:
    text = (raw_tax_text or "").strip().lower()
    if not text:
        return "not_evaluated"
    if re.search(r"\b(uncultured|unidentified|environmental|metagenome|unknown)\b", text):
        return "non_specific_label_review_required"
    return "in_scope_or_not_evaluated"


def _review_flag(
    pident: Optional[float],
    qcovs: Optional[float],
    evalue: Optional[float],
    config: AnnotationConfig,
) -> Tuple[str, str]:
    rules_used = any(
        x is not None for x in (
            config.review_pident_below,
            config.review_qcov_below,
            config.review_evalue_above,
        )
    )
    if not rules_used:
        return "not_evaluated", ""

    reasons: List[str] = []
    if config.review_pident_below is not None and pident is not None and pident < config.review_pident_below:
        reasons.append("pident_below_review_threshold")
    if config.review_qcov_below is not None and qcovs is not None and qcovs < config.review_qcov_below:
        reasons.append("qcov_below_review_threshold")
    if config.review_evalue_above is not None and evalue is not None and evalue > config.review_evalue_above:
        reasons.append("evalue_above_review_threshold")

    if reasons:
        return "review_required", ";".join(reasons)
    return "pass", ""


def _pick_raw_taxonomy_text(row: pd.Series, tax_map: Dict[str, str]) -> str:
    sseqid = str(row.get("sseqid", "") or "").strip()
    if sseqid and sseqid in tax_map:
        return tax_map[sseqid]

    for key in ("sscinames", "stitle"):
        value = str(row.get(key, "") or "").strip()
        if value:
            return value

    return ""


def annotate_rep_seqs_by_marker(
    marker: str,
    rep_fasta: Path,
    db_prefix: Path,
    layout: ReferenceLayout,
    config: AnnotationConfig,
    temp_dir: Path,
) -> pd.DataFrame:
    _require_tool("blastn")

    all_features = list(_read_fasta_records(rep_fasta).keys())
    blast_out = temp_dir / f"{marker}.blast6.tsv"
    outfmt = "6 qseqid sseqid pident length evalue bitscore qcovs sscinames stitle"

    cmd = [
        "blastn",
        "-query", str(rep_fasta),
        "-db", str(db_prefix),
        "-outfmt", outfmt,
        "-max_target_seqs", "25",
        "-max_hsps", "1",
        "-out", str(blast_out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise Step3Error(f"blastn failed for {marker}: {proc.stderr.strip()}")

    tax_map = _load_taxonomy_table(layout) if layout.taxonomy_mode == "table" else {}
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs", "sscinames", "stitle"]

    if blast_out.exists() and blast_out.stat().st_size > 0:
        hits = pd.read_csv(blast_out, sep="\t", names=cols, dtype=str)
        for c in ("pident", "evalue", "bitscore", "qcovs"):
            hits[c] = pd.to_numeric(hits[c], errors="coerce")
    else:
        hits = pd.DataFrame(columns=cols)

    rows: List[Dict[str, object]] = []

    for feature_id in all_features:
        feature_hits = hits[hits["qseqid"] == feature_id].copy() if not hits.empty else pd.DataFrame(columns=cols)

        if feature_hits.empty:
            rows.append(
                {
                    "feature_id": feature_id,
                    "marker": marker,
                    "status": "unassigned",
                    "status_reason": "no_blast_hit",
                    "taxon_rank": "unassigned",
                    "taxon_name": "unassigned",
                    "resolution_flag": "unassigned",
                    "ambiguity_flag": "no_hit",
                    "annotation_scope_flag": "not_evaluated",
                    "low_confidence_flag": "not_evaluated",
                    "low_confidence_reason": "",
                    "confidence": None,
                    "confidence_metric": "bitscore",
                    "pident": None,
                    "evalue": None,
                    "qcovs": None,
                    "sseqid": "",
                    "raw_taxonomy_text": "",
                    "stitle": "",
                }
            )
            continue

        feature_hits = feature_hits.sort_values(
            ["bitscore", "pident", "qcovs"],
            ascending=[False, False, False],
            na_position="last",
        ).reset_index(drop=True)

        top_bitscore = feature_hits.iloc[0]["bitscore"]
        tied = feature_hits[feature_hits["bitscore"] == top_bitscore].copy().reset_index(drop=True)

        resolved_candidates: List[Dict[str, object]] = []
        for _, hit in tied.iterrows():
            raw_tax = _pick_raw_taxonomy_text(hit, tax_map)
            parsed = _parse_taxonomy_from_text(raw_tax)
            resolved = resolve_taxon_rank(parsed, config.rank_resolution_policy)
            resolved_candidates.append(
                {
                    "sseqid": str(hit.get("sseqid", "") or ""),
                    "raw_taxonomy_text": raw_tax,
                    "resolved": resolved,
                    "pident": None if pd.isna(hit["pident"]) else float(hit["pident"]),
                    "evalue": None if pd.isna(hit["evalue"]) else float(hit["evalue"]),
                    "qcovs": None if pd.isna(hit["qcovs"]) else float(hit["qcovs"]),
                    "bitscore": None if pd.isna(hit["bitscore"]) else float(hit["bitscore"]),
                    "stitle": str(hit.get("stitle", "") or ""),
                }
            )

        selected = resolved_candidates[0]
        ambiguity_flag = "single_best_hit"

        if len(resolved_candidates) > 1:
            unique_resolved = {
                (
                    rc["resolved"]["taxon_rank"],
                    rc["resolved"]["taxon_name"],
                    rc["resolved"]["status_reason"],
                )
                for rc in resolved_candidates
            }
            only_unassigned = all(x[0] == "unassigned" for x in unique_resolved)

            if len(unique_resolved) == 1 and not only_unassigned:
                selected = resolved_candidates[0]
                ambiguity_flag = "top_bitscore_tie_same_resolved_taxon"
            else:
                rows.append(
                    {
                        "feature_id": feature_id,
                        "marker": marker,
                        "status": "unassigned",
                        "status_reason": "ambiguous_top_bitscore_tie",
                        "taxon_rank": "unassigned",
                        "taxon_name": "unassigned",
                        "resolution_flag": "unassigned",
                        "ambiguity_flag": "top_bitscore_tie_conflict",
                        "annotation_scope_flag": "review_required",
                        "low_confidence_flag": "review_required",
                        "low_confidence_reason": "ambiguous_top_bitscore_tie",
                        "confidence": resolved_candidates[0]["bitscore"],
                        "confidence_metric": "bitscore",
                        "pident": resolved_candidates[0]["pident"],
                        "evalue": resolved_candidates[0]["evalue"],
                        "qcovs": resolved_candidates[0]["qcovs"],
                        "sseqid": ";".join(str(x["sseqid"]) for x in resolved_candidates),
                        "raw_taxonomy_text": "; ".join(
                            str(x["raw_taxonomy_text"]) for x in resolved_candidates if x["raw_taxonomy_text"]
                        ),
                        "stitle": "ambiguous_multiple_top_hits",
                    }
                )
                continue

        resolved = selected["resolved"]
        review_flag, review_reason = _review_flag(
            pident=selected["pident"],
            qcovs=selected["qcovs"],
            evalue=selected["evalue"],
            config=config,
        )

        status = "assigned" if resolved["taxon_rank"] != "unassigned" else "unassigned"

        rows.append(
            {
                "feature_id": feature_id,
                "marker": marker,
                "status": status,
                "status_reason": resolved["status_reason"],
                "taxon_rank": resolved["taxon_rank"],
                "taxon_name": resolved["taxon_name"],
                "resolution_flag": resolved["resolution_flag"],
                "ambiguity_flag": ambiguity_flag,
                "annotation_scope_flag": _annotation_scope_flag(str(selected["raw_taxonomy_text"])),
                "low_confidence_flag": review_flag,
                "low_confidence_reason": review_reason,
                "confidence": selected["bitscore"],
                "confidence_metric": "bitscore",
                "pident": selected["pident"],
                "evalue": selected["evalue"],
                "qcovs": selected["qcovs"],
                "sseqid": str(selected["sseqid"]),
                "raw_taxonomy_text": str(selected["raw_taxonomy_text"]),
                "stitle": str(selected["stitle"]),
            }
        )

    ann_df = pd.DataFrame(rows)
    if ann_df.empty:
        raise Step3Error(f"No annotation rows generated for {marker}")

    return ann_df


def join_annotation_to_asv(
    marker: str,
    asv_df: pd.DataFrame,
    feature_col: str,
    sample_cols: List[str],
    seq_map: Dict[str, str],
    ann_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    base = asv_df.rename(columns={feature_col: "feature_id"}).copy()
    base["feature_id"] = base["feature_id"].astype(str).str.strip()
    base["sequence"] = base["feature_id"].map(seq_map)

    if base["sequence"].isna().any():
        missing = base.loc[base["sequence"].isna(), "feature_id"].tolist()[:5]
        raise Step3Error(f"Missing sequence for feature_id in {marker}: {missing}")

    ann_key = ann_df.drop_duplicates(subset=["feature_id"]).copy()
    merged = base.merge(ann_key, on="feature_id", how="left")

    if merged["marker"].isna().any():
        missing = merged.loc[merged["marker"].isna(), "feature_id"].tolist()[:5]
        raise Step3Error(f"Missing annotation rows after join for {marker}: {missing}")

    for col in sample_cols:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)
        if (merged[col] < 0).any():
            raise Step3Error(f"Negative reads detected in ASV table for {marker}, column {col}")

    long_df = merged.melt(
        id_vars=[
            "feature_id",
            "sequence",
            "marker",
            "status",
            "status_reason",
            "taxon_rank",
            "taxon_name",
            "resolution_flag",
            "ambiguity_flag",
            "annotation_scope_flag",
            "low_confidence_flag",
            "low_confidence_reason",
            "confidence",
            "confidence_metric",
            "pident",
            "evalue",
            "qcovs",
            "sseqid",
            "raw_taxonomy_text",
            "stitle",
        ],
        value_vars=sample_cols,
        var_name="sample_id",
        value_name="reads",
    )

    long_df["reads"] = pd.to_numeric(long_df["reads"], errors="coerce").fillna(0)
    long_df = long_df[long_df["reads"] > 0].copy()
    long_df["reads"] = long_df["reads"].round().astype(int)

    candidate_cols = [
        "sample_id",
        "feature_id",
        "sequence",
        "marker",
        "reads",
        "taxon_rank",
        "taxon_name",
        "confidence",
        "confidence_metric",
        "status",
        "status_reason",
        "resolution_flag",
        "ambiguity_flag",
        "annotation_scope_flag",
        "low_confidence_flag",
        "low_confidence_reason",
        "pident",
        "evalue",
        "qcovs",
        "sseqid",
        "raw_taxonomy_text",
        "stitle",
    ]
    candidate_df = long_df[candidate_cols].sort_values(
        ["sample_id", "marker", "feature_id"],
        kind="stable",
    ).reset_index(drop=True)

    unassigned_df = ann_df[ann_df["status"] == "unassigned"].copy()
    unassigned_df["sequence"] = unassigned_df["feature_id"].map(seq_map)
    unassigned_cols = [
        "feature_id",
        "sequence",
        "marker",
        "status",
        "status_reason",
        "taxon_rank",
        "taxon_name",
        "resolution_flag",
        "ambiguity_flag",
        "annotation_scope_flag",
        "low_confidence_flag",
        "low_confidence_reason",
        "confidence",
        "confidence_metric",
        "pident",
        "evalue",
        "qcovs",
        "sseqid",
        "raw_taxonomy_text",
        "stitle",
    ]
    unassigned_df = unassigned_df[unassigned_cols].drop_duplicates().reset_index(drop=True)

    return candidate_df, unassigned_df


def build_annotation_rank_summary(ann_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        ann_df.groupby(
            [
                "marker",
                "taxon_rank",
                "status",
                "resolution_flag",
                "ambiguity_flag",
                "annotation_scope_flag",
                "low_confidence_flag",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="n_features")
        .sort_values(
            ["marker", "status", "taxon_rank", "resolution_flag"],
            kind="stable",
        )
        .reset_index(drop=True)
    )
    return summary


def write_step3_outputs(
    output_dir: Path,
    candidate_tables: Dict[str, pd.DataFrame],
    confidence_df: pd.DataFrame,
    rank_summary_df: pd.DataFrame,
    unassigned_df: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    candidate_tables["ITS1"].to_csv(
        output_dir / "ITS1_candidate_lineages.tsv",
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    candidate_tables["ITS2"].to_csv(
        output_dir / "ITS2_candidate_lineages.tsv",
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    confidence_df.to_csv(
        output_dir / "annotation_confidence.tsv",
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    rank_summary_df.to_csv(
        output_dir / "annotation_rank_summary.tsv",
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    unassigned_df.to_csv(
        output_dir / "unassigned_features.tsv",
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )

    for name in REQUIRED_OUTPUTS:
        if not (output_dir / name).exists():
            raise Step3Error(f"Failed to generate required output: {name}")


def main() -> int:
    args = parse_args()

    result_dir = Path(args.result_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else result_dir
    reference_root = Path(args.reference_root).resolve()

    annotation_cfg = AnnotationConfig(
        ambiguity_policy=args.ambiguity_policy,
        rank_resolution_policy=args.rank_resolution_policy,
        allow_unassigned=args.allow_unassigned,
        sample_id_regex=args.sample_id_regex,
        review_pident_below=args.review_pident_below,
        review_qcov_below=args.review_qcov_below,
        review_evalue_above=args.review_evalue_above,
    )

    print("[Step3] Path configuration")
    print(f"  Windows result dir: {args.windows_result_dir}")
    print(f"  WSL result dir:     {args.result_dir}")
    print(f"  Active result dir:  {result_dir}")
    print(f"  Reference root dir: {reference_root}")
    print(f"  Output dir:         {output_dir}")
    print(f"  Method stack:       {METHOD_STACK['name']}")

    try:
        _require_tool("blastn")
        _require_tool("makeblastdb")

        file_map = load_step2_outputs(result_dir)
        validated = validate_step2_input_schema(file_map, args.sample_id_regex)

        with tempfile.TemporaryDirectory(prefix="step03_refcheck_") as tmp:
            tmpdir = Path(tmp)

            layouts = locate_and_validate_reference_resources(reference_root, tmpdir)
            provenance_df = freeze_reference_provenance(layouts)

            index_root = tmpdir / "blast_index"
            db_its1 = build_reference_index("ITS1", layouts["ITS1"], index_root)
            db_its2 = build_reference_index("ITS2", layouts["ITS2"], index_root)

            ann_its1 = annotate_rep_seqs_by_marker(
                marker="ITS1",
                rep_fasta=file_map["ITS1_rep_seqs.fasta"],
                db_prefix=db_its1,
                layout=layouts["ITS1"],
                config=annotation_cfg,
                temp_dir=tmpdir,
            )
            ann_its2 = annotate_rep_seqs_by_marker(
                marker="ITS2",
                rep_fasta=file_map["ITS2_rep_seqs.fasta"],
                db_prefix=db_its2,
                layout=layouts["ITS2"],
                config=annotation_cfg,
                temp_dir=tmpdir,
            )

            cand_its1, unassigned_its1 = join_annotation_to_asv(
                marker="ITS1",
                asv_df=validated["ITS1_asv_df"],
                feature_col=validated["ITS1_feature_col"],
                sample_cols=validated["ITS1_sample_cols"],
                seq_map=validated["ITS1_seq_map"],
                ann_df=ann_its1,
            )
            cand_its2, unassigned_its2 = join_annotation_to_asv(
                marker="ITS2",
                asv_df=validated["ITS2_asv_df"],
                feature_col=validated["ITS2_feature_col"],
                sample_cols=validated["ITS2_sample_cols"],
                seq_map=validated["ITS2_seq_map"],
                ann_df=ann_its2,
            )

            confidence_df = pd.concat([ann_its1, ann_its2], ignore_index=True)
            confidence_df = confidence_df.merge(provenance_df, on="marker", how="left")
            confidence_df = confidence_df.sort_values(["marker", "feature_id"], kind="stable").reset_index(drop=True)

            rank_summary_df = build_annotation_rank_summary(pd.concat([ann_its1, ann_its2], ignore_index=True))

            unassigned_df = pd.concat([unassigned_its1, unassigned_its2], ignore_index=True)
            unassigned_df = unassigned_df.merge(provenance_df, on="marker", how="left")
            unassigned_df = unassigned_df.sort_values(["marker", "feature_id"], kind="stable").reset_index(drop=True)

            write_step3_outputs(
                output_dir=output_dir,
                candidate_tables={"ITS1": cand_its1, "ITS2": cand_its2},
                confidence_df=confidence_df,
                rank_summary_df=rank_summary_df,
                unassigned_df=unassigned_df,
            )

        print("[Step3] Completed successfully. Required outputs written.")
        return 0

    except Step3Error as e:
        print(f"[Step3] STOP: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
