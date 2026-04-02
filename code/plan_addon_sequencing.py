#!/usr/bin/env python3
"""
Plan add-on UG sequencing: select below-median samples per omic, convert effective
deficit to raw reads using demux stats, scale to 5B SE reads total.

ATAC-seq: median effective reads excludes sample 31_ZI_4_3 (contaminated original);
that sample is always included with raw need for a full median worth of effective reads
(prompt.md rule 6). Outputs default to add_on_sequencing_plan/data/output/:
  - CSV read plan: effective + raw read counts as integers; scaling/proportions 4 decimals (rule 8)
  - One PNG per omic: <prefix>_*_<OMIC>_effective_reads_before_after.png (before vs after bars + medians)
  - *_pooling.xlsx with four sheets: ATAC-seq_1, ATAC-seq_2, RNA-seq, sRNA-seq
    (column A = proportion_reads_within_pool; columns F–G = demux molarity correction per prompt;
    undiluted volume uses corrected nM in column G).

Run from repo (or this folder) with uv:
  cd add_on_sequencing_plan && uv run python code/plan_addon_sequencing.py

Defaults read inputs from ``data/input/`` with paths mirroring the parent multiomics repo
(see README). Optional local ``prompt.md`` is not required to run the script.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from effective_reads_plots import write_before_after_effective_plots
from pooling_workbook import build_pooling_workbook, sorted_sample_ids_from_atac

TARGET_SE_READS = 5_000_000_000

# ATAC: original data discarded (contamination); median excludes this sample; it is always
# scheduled with deficit_effective = median (replace with ~median effective reads total).
ATAC_SPECIAL_SAMPLE = "31_ZI_4_3"
ATAC_SPECIAL_CONC_NG_UL = 4.06
ATAC_SPECIAL_SIZE_BP = 663.0


def bundle_root() -> Path:
    """add_on_sequencing_plan/ (parent of ``code/``)."""
    return Path(__file__).resolve().parent.parent


def _input_relpath(*parts: str) -> Path:
    return bundle_root().joinpath("data", "input", *parts)


def load_demux(path: Path) -> tuple[dict[str, int], float]:
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    total = float(d["total_reads"])
    assigned = float(d["assigned_reads"])
    p = assigned / total if total else 0.0
    counts = {str(k): int(v) for k, v in d["sample_counts"].items()}
    return counts, p


def load_atac_effective(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"Sample": "sample_raw", "Aligned Reads": "effective"})
    df["sample_id"] = df["sample_raw"].astype(str).str.replace(
        r"_REP1$", "", regex=True
    )
    df["effective"] = pd.to_numeric(df["effective"], errors="coerce")
    return df[["sample_id", "effective"]]


def load_rna_effective(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in ("uniquely_mapped", "multimapped"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["effective"] = df["uniquely_mapped"] + df["multimapped"]
    df["sample_id"] = df["Sample"].astype(str)
    return df[["sample_id", "effective"]]


def load_srna_effective(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["effective"] = pd.to_numeric(df["miRNA"], errors="coerce")
    df["sample_id"] = df["Sample"].astype(str)
    return df[["sample_id", "effective"]]


def sort_plan_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Within each omic: largest deficit first; ATAC special sample first (its own pool)."""
    blocks: list[pd.DataFrame] = []
    for omic in ("ATAC-seq", "RNA-seq", "sRNA-seq"):
        sub = df[df["omics_type"] == omic]
        if sub.empty:
            continue
        if omic == "ATAC-seq":
            rest = sub[sub["sample_id"] != ATAC_SPECIAL_SAMPLE].copy()
            sp = sub[sub["sample_id"] == ATAC_SPECIAL_SAMPLE].copy()
            rest = rest.sort_values("deficit_effective", ascending=False)
            sub = pd.concat([sp, rest], ignore_index=True)
        else:
            sub = sub.sort_values("deficit_effective", ascending=False)
        blocks.append(sub)
    return pd.concat(blocks, ignore_index=True)


def allocate_integer_reads(raw_need: np.ndarray, target: int) -> np.ndarray:
    """Largest-remainder allocation so integer counts sum exactly to target."""
    raw_need = np.asarray(raw_need, dtype=np.float64)
    total = float(raw_need.sum())
    if total <= 0:
        return np.zeros_like(raw_need, dtype=np.int64)
    quotas = raw_need / total * target
    floors = np.floor(quotas).astype(np.int64)
    remainder = int(target - floors.sum())
    if remainder > 0:
        frac = quotas - floors
        order = np.argsort(-frac)
        for i in range(remainder):
            floors[order[i % len(order)]] += 1
    return floors


def select_and_compute(
    omics_label: str,
    eff: pd.DataFrame,
    sample_counts: dict[str, int],
    p: float,
    max_samples: int | None,
) -> tuple[pd.DataFrame, list[str]]:
    """Median, below-median, rank by deficit, optional cap, raw_need."""
    n = len(eff)
    med = float(eff["effective"].median())
    if p <= 0:
        return pd.DataFrame(), [f"{omics_label}: demux assigned_reads/total_reads is 0"]

    below = eff[eff["effective"] < med].copy()
    below["deficit_effective"] = med - below["effective"]
    below = below.sort_values("deficit_effective", ascending=False)
    if max_samples is not None:
        below = below.head(int(max_samples))

    warnings: list[str] = []
    rows: list[dict] = []
    for _, r in below.iterrows():
        sid = r["sample_id"]
        e = float(r["effective"])
        delta_e = float(r["deficit_effective"])
        if e <= 0:
            warnings.append(f"{omics_label} {sid}: effective reads <= 0, skipped")
            continue
        if sid not in sample_counts:
            warnings.append(
                f"{omics_label} {sid}: missing from demux sample_counts, skipped"
            )
            continue
        R = float(sample_counts[sid])
        raw_need = delta_e * R / (e * p)
        rows.append(
            {
                "omics_type": omics_label,
                "sample_id": sid,
                "effective_reads": e,
                "median_effective": med,
                "deficit_effective": delta_e,
                "raw_reads_needed": raw_need,
            }
        )

    out = pd.DataFrame(rows)
    return out, warnings


def select_and_compute_atac(
    eff: pd.DataFrame,
    sample_counts: dict[str, int],
    p: float,
    max_samples: int | None,
) -> tuple[pd.DataFrame, list[str]]:
    """ATAC: median excludes ATAC_SPECIAL_SAMPLE; that sample always gets raw_need for a full median worth of effective reads."""
    warnings: list[str] = []
    med_series = eff.loc[eff["sample_id"] != ATAC_SPECIAL_SAMPLE, "effective"]
    if med_series.empty:
        return pd.DataFrame(), [f"ATAC-seq: no samples left to compute median (exclude {ATAC_SPECIAL_SAMPLE})"]
    med = float(med_series.median())

    if p <= 0:
        return pd.DataFrame(), ["ATAC-seq: demux assigned_reads/total_reads is 0"]

    below = eff[
        (eff["effective"] < med) & (eff["sample_id"] != ATAC_SPECIAL_SAMPLE)
    ].copy()
    below["deficit_effective"] = med - below["effective"]
    below = below.sort_values("deficit_effective", ascending=False)
    if max_samples is not None:
        below = below.head(int(max_samples))

    rows: list[dict] = []
    for _, r in below.iterrows():
        sid = r["sample_id"]
        e = float(r["effective"])
        delta_e = float(r["deficit_effective"])
        if e <= 0:
            warnings.append(f"ATAC-seq {sid}: effective reads <= 0, skipped")
            continue
        if sid not in sample_counts:
            warnings.append(f"ATAC-seq {sid}: missing from demux sample_counts, skipped")
            continue
        R = float(sample_counts[sid])
        raw_need = delta_e * R / (e * p)
        rows.append(
            {
                "omics_type": "ATAC-seq",
                "sample_id": sid,
                "effective_reads": e,
                "median_effective": med,
                "deficit_effective": delta_e,
                "raw_reads_needed": raw_need,
            }
        )

    spec = eff[eff["sample_id"] == ATAC_SPECIAL_SAMPLE]
    if spec.empty:
        warnings.append(
            f"ATAC-seq: {ATAC_SPECIAL_SAMPLE} not in MultiQC table; special row omitted"
        )
    else:
        e = float(spec["effective"].iloc[0])
        if e <= 0:
            warnings.append(
                f"ATAC-seq {ATAC_SPECIAL_SAMPLE}: effective reads <= 0; special row omitted"
            )
        elif ATAC_SPECIAL_SAMPLE not in sample_counts:
            warnings.append(
                f"ATAC-seq {ATAC_SPECIAL_SAMPLE}: missing from demux sample_counts; special row omitted"
            )
        else:
            R = float(sample_counts[ATAC_SPECIAL_SAMPLE])
            delta_e = med  # median - 0
            raw_need = delta_e * R / (e * p)
            rows.append(
                {
                    "omics_type": "ATAC-seq",
                    "sample_id": ATAC_SPECIAL_SAMPLE,
                    "effective_reads": e,
                    "median_effective": med,
                    "deficit_effective": delta_e,
                    "raw_reads_needed": raw_need,
                }
            )

    return pd.DataFrame(rows), warnings


def main() -> None:
    ap = argparse.ArgumentParser(description="Plan add-on UG sequencing (5B SE cap).")
    ap.add_argument(
        "--atac-multiqc",
        type=Path,
        default=_input_relpath(
            "atac-seq",
            "data",
            "output",
            "multiqc",
            "broad_peak",
            "ATAC-seq_36_samples_UG_SE_multiqc_report_data",
            "mqc_picard_alignment_summary__name_Aligned_Reads_ylab_Reads_cpswitch_counts_label_Number_of_Reads_.txt",
        ),
    )
    ap.add_argument(
        "--rna-multiqc",
        type=Path,
        default=_input_relpath(
            "rna-seq",
            "data",
            "output",
            "multiqc",
            "star_rsem",
            "multiqc_report_data",
            "multiqc_star.txt",
        ),
    )
    ap.add_argument(
        "--srna-multiqc",
        type=Path,
        default=_input_relpath(
            "srna-seq",
            "data",
            "output",
            "multiqc",
            "smallRNA-seq-analysis-of-Dmel-cold-adaptation-samples_multiqc_report_data",
            "mirtrace_rna_categories_plot.txt",
        ),
    )
    ap.add_argument(
        "--atac-demux",
        type=Path,
        default=_input_relpath(
            "data_UG_seq",
            "data",
            "de-multiplexed_reads",
            "atacseq",
            "demux_stats.json",
        ),
    )
    ap.add_argument(
        "--rna-demux",
        type=Path,
        default=_input_relpath(
            "data_UG_seq",
            "data",
            "de-multiplexed_reads",
            "rnaseq",
            "demux_stats.json",
        ),
    )
    ap.add_argument(
        "--srna-demux",
        type=Path,
        default=_input_relpath(
            "data_UG_seq",
            "data",
            "de-multiplexed_reads",
            "smallrna",
            "demux_stats.json",
        ),
    )
    ap.add_argument(
        "--n-atac",
        type=int,
        default=None,
        metavar="N",
        help="Max ATAC samples to include (below median, by priority). Default: all.",
    )
    ap.add_argument(
        "--n-rna",
        type=int,
        default=None,
        metavar="N",
        help="Max RNA samples to include. Default: all below median.",
    )
    ap.add_argument(
        "--n-srna",
        type=int,
        default=None,
        metavar="N",
        help="Max sRNA samples to include. Default: all below median.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=bundle_root() / "data" / "output",
        help="Directory for output CSV.",
    )
    ap.add_argument(
        "--prefix",
        default="addon_sequencing",
        help="Output filename prefix (before _<na>_<nr>_<ns>.csv).",
    )
    ap.add_argument(
        "--pooling-rna-xlsx",
        type=Path,
        default=bundle_root() / "data" / "input" / "2025-12-03_pooling_WS.xlsx",
        help="RNA + small RNA pooling template (.xlsx).",
    )
    ap.add_argument(
        "--pooling-atac-xlsx",
        type=Path,
        default=bundle_root() / "data" / "input" / "pooling_for_ngs_all_samples.xlsx",
        help="ATAC pooling template (.xlsx), sheet Pool1_2.",
    )
    ap.add_argument(
        "--skip-pooling-xlsx",
        action="store_true",
        help="Do not write the four-sheet pooling workbook.",
    )
    args = ap.parse_args()

    atac_eff = load_atac_effective(args.atac_multiqc)
    rna_eff = load_rna_effective(args.rna_multiqc)
    srna_eff = load_srna_effective(args.srna_multiqc)

    expected_n = 36
    for name, df in (
        ("ATAC-seq", atac_eff),
        ("RNA-seq", rna_eff),
        ("sRNA-seq", srna_eff),
    ):
        n = len(df)
        print(f"{name}: {n} samples in MultiQC table.")
        if n != expected_n:
            print(
                f"  WARNING: expected {expected_n} samples, got {n}.",
                file=sys.stderr,
            )

    atac_counts, atac_p = load_demux(args.atac_demux)
    rna_counts, rna_p = load_demux(args.rna_demux)
    srna_counts, srna_p = load_demux(args.srna_demux)

    parts: list[pd.DataFrame] = []
    all_warnings: list[str] = []

    sub_atac, w_atac = select_and_compute_atac(
        atac_eff, atac_counts, atac_p, args.n_atac
    )
    all_warnings.extend(w_atac)
    if not sub_atac.empty:
        parts.append(sub_atac)

    for label, eff, counts, p, n_cap in (
        ("RNA-seq", rna_eff, rna_counts, rna_p, args.n_rna),
        ("sRNA-seq", srna_eff, srna_counts, srna_p, args.n_srna),
    ):
        sub, w = select_and_compute(label, eff, counts, p, n_cap)
        all_warnings.extend(w)
        if not sub.empty:
            parts.append(sub)

    for w in all_warnings:
        print(w, file=sys.stderr)

    if not parts:
        print("No samples selected; nothing to write.", file=sys.stderr)
        sys.exit(1)

    out = pd.concat(parts, ignore_index=True)
    out = sort_plan_rows(out)

    pool_name = np.where(
        out["omics_type"].to_numpy() == "ATAC-seq",
        np.where(
            out["sample_id"].to_numpy().astype(str) == ATAC_SPECIAL_SAMPLE,
            "ATAC-seq_1",
            "ATAC-seq_2",
        ),
        out["omics_type"].to_numpy().astype(str),
    )
    out["pool_name"] = pool_name

    total_raw = float(out["raw_reads_needed"].sum())
    if total_raw <= 0:
        print("Total raw_reads_needed is 0.", file=sys.stderr)
        sys.exit(1)

    scaling = TARGET_SE_READS / total_raw
    out["scaling_factor"] = scaling

    # Read share within each sequencing pool (UG pool), from pre-scale needs
    pool_sums = out.groupby("pool_name")["raw_reads_needed"].transform("sum")
    out["proportion_within_pool"] = out["raw_reads_needed"] / pool_sums
    out["proportion_across_omics"] = out["raw_reads_needed"] / total_raw

    raw_arr = out["raw_reads_needed"].to_numpy()
    int_reads = allocate_integer_reads(raw_arr, TARGET_SE_READS)
    out["raw_reads_to_sequence"] = int_reads

    assert int(out["raw_reads_to_sequence"].sum()) == TARGET_SE_READS

    rn = out["raw_reads_needed"].to_numpy(dtype=np.float64)
    rts = out["raw_reads_to_sequence"].to_numpy(dtype=np.float64)
    de = out["deficit_effective"].to_numpy(dtype=np.float64)
    e0 = out["effective_reads"].to_numpy(dtype=np.float64)
    eff_gain = np.zeros(len(out), dtype=np.float64)
    need_pos = rn > 0
    eff_gain[need_pos] = rts[need_pos] / rn[need_pos] * de[need_pos]
    total_eff_after = np.rint(e0 + eff_gain).astype(np.int64)

    spec_mask = (out["omics_type"].to_numpy() == "ATAC-seq") & (
        out["sample_id"].to_numpy().astype(str) == ATAC_SPECIAL_SAMPLE
    )
    eff_export = np.rint(out["effective_reads"].to_numpy()).astype(np.int64)
    eff_export[spec_mask] = 0
    total_eff_export = total_eff_after.copy()
    total_eff_export[spec_mask] = np.rint(eff_gain[spec_mask]).astype(np.int64)

    na = int((out["omics_type"] == "ATAC-seq").sum())
    nr = int((out["omics_type"] == "RNA-seq").sum())
    ns = int((out["omics_type"] == "sRNA-seq").sum())

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / f"{args.prefix}_{na}_{nr}_{ns}.csv"

    export = pd.DataFrame(
        {
            "omics_type": out["omics_type"],
            "pool_name": out["pool_name"],
            "sample_id": out["sample_id"],
            "effective_reads": eff_export,
            "median_effective_reads": np.rint(out["median_effective"].to_numpy()).astype(
                np.int64
            ),
            "effective_reads_to_reach_median": np.rint(
                out["deficit_effective"].to_numpy()
            ).astype(np.int64),
            "raw_reads_needed": np.rint(out["raw_reads_needed"].to_numpy()).astype(
                np.int64
            ),
            "raw_reads_to_sequence": out["raw_reads_to_sequence"].astype(np.int64),
            "scaling_factor": out["scaling_factor"].round(4),
            "total_effective_reads_expected_after_resequencing": total_eff_export,
            "proportion_within_pool": out["proportion_within_pool"].round(4),
            "proportion_across_omics": out["proportion_across_omics"].round(4),
        }
    )
    export.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(export)} rows, sum raw_reads_to_sequence = {TARGET_SE_READS}).")

    plot_paths = write_before_after_effective_plots(
        export,
        out_path.with_suffix(""),
        full_effective_by_omic={
            "ATAC-seq": atac_eff[["sample_id", "effective"]].copy(),
            "RNA-seq": rna_eff[["sample_id", "effective"]].copy(),
            "sRNA-seq": srna_eff[["sample_id", "effective"]].copy(),
        },
    )
    for pp in plot_paths:
        print(f"Wrote {pp}")

    if not args.skip_pooling_xlsx:
        if not args.pooling_rna_xlsx.is_file() or not args.pooling_atac_xlsx.is_file():
            print(
                "WARNING: pooling template xlsx missing; skip pooling workbook.",
                file=sys.stderr,
            )
        else:
            sorted_ids = sorted_sample_ids_from_atac(atac_eff)

            at1 = out[out["pool_name"] == "ATAC-seq_1"]
            at2 = out[out["pool_name"] == "ATAC-seq_2"]
            rna_sub = out[out["pool_name"] == "RNA-seq"]
            srna_sub = out[out["pool_name"] == "sRNA-seq"]
            pools = {
                "ATAC-seq_1": {
                    "sample_ids": at1["sample_id"].astype(str).tolist(),
                    "proportion_reads_within_pool": at1["proportion_within_pool"]
                    .astype(float)
                    .tolist(),
                    "overrides": {
                        ATAC_SPECIAL_SAMPLE: {
                            "conc_ng_ul": ATAC_SPECIAL_CONC_NG_UL,
                            "size_bp": ATAC_SPECIAL_SIZE_BP,
                        }
                    },
                },
                "ATAC-seq_2": {
                    "sample_ids": at2["sample_id"].astype(str).tolist(),
                    "proportion_reads_within_pool": at2["proportion_within_pool"]
                    .astype(float)
                    .tolist(),
                    "overrides": {},
                },
                "RNA-seq": {
                    "sample_ids": rna_sub["sample_id"].astype(str).tolist(),
                    "proportion_reads_within_pool": rna_sub["proportion_within_pool"]
                    .astype(float)
                    .tolist(),
                    "overrides": {},
                },
                "sRNA-seq": {
                    "sample_ids": srna_sub["sample_id"].astype(str).tolist(),
                    "proportion_reads_within_pool": srna_sub["proportion_within_pool"]
                    .astype(float)
                    .tolist(),
                    "overrides": {},
                },
            }
            xlsx_path = args.output_dir / f"{args.prefix}_{na}_{nr}_{ns}_pooling.xlsx"
            build_pooling_workbook(
                xlsx_path,
                args.pooling_atac_xlsx,
                args.pooling_rna_xlsx,
                sorted_ids,
                pools,
                {
                    "ATAC-seq": atac_counts,
                    "RNA-seq": rna_counts,
                    "sRNA-seq": srna_counts,
                },
            )
            print(f"Wrote {xlsx_path}")


if __name__ == "__main__":
    main()
