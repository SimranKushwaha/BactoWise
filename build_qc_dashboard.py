#!/usr/bin/env python3
"""
Build a QC dashboard from an existing Bactowise v3 output folder.

Expected input:
  - Master_Table_Annotation.xlsx

Default output:
  - QC-Dashboard.html
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

LOCAL_CONFLICT_REVIEW_MAX_CONFIDENCE = 85.0
REVIEWMAYBE_PRIORITY_MAX_CONFIDENCE = 60.0

EXPORTED_TO_INTERNAL_COLUMNS = {
    "Contig_ID": "Sequence_ID",
    "Consensus-Gene-Name": "Best-Gene",
    "Likely-Gene-Synonyms": "Gene-Putative-Synonyms",
    "Bactowise-Confidence-Score": "Confidence-Score",
    "PGAP_Predicted_Pseudogene": "Pseudogene",
    "Gene_Names_from_Input_Gffs": "Gene",
    "Product_from_Input_Gffs": "Product",
    "Locus_Tag_from_Input_Gffs": "Original_Locus_Tag",
}


def normalise_input_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Accept either pipeline-internal headers or user-facing Excel headers."""
    df = df.rename(columns={k: v for k, v in EXPORTED_TO_INTERNAL_COLUMNS.items() if k in df.columns})
    if "Status" in df.columns:
        df["Status"] = df["Status"].replace({"ReviewMaybe": "ScreenMaybe"})
    if "Pseudogene" in df.columns:
        df["Pseudogene"] = (
            df["Pseudogene"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(["true", "1", "yes"])
        )
    return df


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


def series_table(title: str, series: pd.Series) -> str:
    rows = "".join(
        f"<tr><td>{label}</td><td>{int(value)}</td></tr>"
        for label, value in series.items()
    )
    return f"<section><h2>{title}</h2><table><tbody>{rows}</tbody></table></section>"


def build_dashboard_html(df: pd.DataFrame) -> str:
    df = blank_screenmaybe_hp_gene_display(df)
    priority_screenmaybe_mask = (
        (df["Status"] == "ScreenMaybe") &
        (df["Confidence-Score"] <= REVIEWMAYBE_PRIORITY_MAX_CONFIDENCE)
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

    summary = {
        "Total features": int(len(df)),
        "High confidence (>=85)": int((df["Confidence-Score"] >= 85).sum()),
        "Medium confidence (60-84)": int(
            ((df["Confidence-Score"] >= 60) & (df["Confidence-Score"] < 85)).sum()
        ),
        "Low confidence (<60)": int((df["Confidence-Score"] < 60).sum()),
        "Mean confidence": round(float(df["Confidence-Score"].mean()), 2) if len(df) else 0.0,
        "Worth reviewing annotations": primary_manual_check_count,
    }
    summary_html = "".join(f"<li><strong>{k}:</strong> {v}</li>" for k, v in summary.items())
    local_conflict_tags = set()
    grouped = {
        key: group.sort_values(["Start", "End", "Locus_Tag"])
        for key, group in df.groupby(["Sequence_ID", "Direction", "Type"], dropna=False)
    }
    for _, group in grouped.items():
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
            ("ReviewMaybe (<=60)", grouped_counts["ReviewMaybe (<=60)"], f"Lower-confidence ReviewMaybe rows prioritised for review at confidence <= {REVIEWMAYBE_PRIORITY_MAX_CONFIDENCE:.0f}."),
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
        manual_df["check_group"] = "Local Conflict"
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

    return f"""<!doctype html>
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Build QC-Dashboard.html from an existing output table.")
    parser.add_argument(
        "--output-folder",
        required=True,
        help="Folder containing Master_Table_Annotation.xlsx",
    )
    parser.add_argument(
        "--input-table",
        default=None,
        help="Optional explicit path to Master_Table_Annotation.xlsx",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional explicit output HTML path (default: <output-folder>/QC-Dashboard.html)",
    )
    args = parser.parse_args()

    output_folder = Path(args.output_folder)
    input_table = Path(args.input_table) if args.input_table else output_folder / "Master_Table_Annotation.xlsx"
    out_path = Path(args.out) if args.out else output_folder / "QC-Dashboard.html"

    if not input_table.exists():
        raise FileNotFoundError(f"Input table not found: {input_table}")

    df = pd.read_excel(input_table, sheet_name="Master_Table")
    df = normalise_input_columns(df)
    required = ["Type", "Status", "Source-of-Gene-Name", "Confidence-Score"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input table missing required columns: {', '.join(missing)}")

    html_doc = build_dashboard_html(df)
    out_path.write_text(html_doc, encoding="utf-8")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
