"""Bar plots: effective reads before vs expected after re-sequencing (per omic)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

BAR_BEFORE = "#4C72B0"
BAR_AFTER = "#55A868"


def _lighter_line_color(rgb_hex: str, mix_white: float = 0.5) -> str:
    """Tint bar color toward white for median reference lines."""
    r, g, b = mcolors.to_rgb(rgb_hex)
    r = r + (1.0 - r) * mix_white
    g = g + (1.0 - g) * mix_white
    b = b + (1.0 - b) * mix_white
    return mcolors.to_hex((r, g, b))


def write_before_after_effective_plots(
    df: pd.DataFrame,
    plot_stem: Path,
    *,
    full_effective_by_omic: dict[str, pd.DataFrame] | None = None,
) -> list[Path]:
    """
    One PNG per omic in df: grouped bars + median reference lines.

    When ``full_effective_by_omic`` is provided (``sample_id``, ``effective`` per omic),
    samples not in the plan are drawn **before-only** (no after bar). Otherwise only
    plan rows are shown (legacy behavior).
    """
    required = (
        "omics_type",
        "sample_id",
        "effective_reads",
        "median_effective_reads",
        "total_effective_reads_expected_after_resequencing",
    )
    for c in required:
        if c not in df.columns:
            raise ValueError(f"plot export missing column {c!r}")

    line_before = _lighter_line_color(BAR_BEFORE)
    line_after = _lighter_line_color(BAR_AFTER)

    written: list[Path] = []
    for omic in ("ATAC-seq", "RNA-seq", "sRNA-seq"):
        sub = df[df["omics_type"] == omic].copy()
        if sub.empty:
            continue

        m_before = float(sub["median_effective_reads"].iloc[0])
        after_vals = sub["total_effective_reads_expected_after_resequencing"].to_numpy(
            dtype=np.float64
        )
        m_after = float(np.median(after_vals))

        after_map = dict(
            zip(
                sub["sample_id"].astype(str),
                sub["total_effective_reads_expected_after_resequencing"].astype(float),
            )
        )
        plan_ids = set(after_map)
        plan_before_map = dict(
            zip(sub["sample_id"].astype(str), sub["effective_reads"].astype(float))
        )

        if full_effective_by_omic and omic in full_effective_by_omic:
            full = full_effective_by_omic[omic]
            if "sample_id" not in full.columns or "effective" not in full.columns:
                raise ValueError(f"full_effective_by_omic[{omic!r}] needs sample_id, effective")
            eff_full = pd.to_numeric(full["effective"], errors="coerce").fillna(0.0)
            bef_map = dict(
                zip(full["sample_id"].astype(str), eff_full.to_numpy(dtype=np.float64))
            )
            plan_order = sub["sample_id"].astype(str).tolist()
            rest = [
                s
                for s in full["sample_id"].astype(str).tolist()
                if s not in plan_ids
            ]
            rest.sort(key=lambda s: (float(bef_map[s]), str(s)))
            samples = plan_order + rest
        else:
            bef_map = dict(
                zip(
                    sub["sample_id"].astype(str),
                    sub["effective_reads"].astype(float),
                )
            )
            samples = sub["sample_id"].astype(str).tolist()

        n = len(samples)
        x = np.arange(n)
        w = 0.36
        fig_w = max(7.5, 0.22 * n + 2.5)
        fig, ax = plt.subplots(figsize=(fig_w, 5.5))

        for i, sid in enumerate(samples):
            if sid in plan_ids:
                b = float(plan_before_map[sid])
            else:
                b = float(bef_map[sid])
            if sid in after_map:
                a = float(after_map[sid])
                ax.bar(
                    i - w / 2,
                    b,
                    width=w,
                    color=BAR_BEFORE,
                )
                ax.bar(
                    i + w / 2,
                    a,
                    width=w,
                    color=BAR_AFTER,
                )
            else:
                ax.bar(i, b, width=w, color=BAR_BEFORE)

        ax.axhline(
            m_before,
            color=line_before,
            linestyle="--",
            linewidth=1.4,
        )
        ax.axhline(
            m_after,
            color=line_after,
            linestyle="--",
            linewidth=1.4,
        )

        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="None",
                markersize=9,
                markerfacecolor=BAR_BEFORE,
                markeredgecolor=BAR_BEFORE,
                label="Before (current)",
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="None",
                markersize=9,
                markerfacecolor=BAR_AFTER,
                markeredgecolor=BAR_AFTER,
                label="After (expected, plan only)",
            ),
            Line2D(
                [0],
                [0],
                color=line_before,
                linestyle="--",
                linewidth=1.4,
                label=f"Median before sequencing ({m_before:,.0f})",
            ),
            Line2D(
                [0],
                [0],
                color=line_after,
                linestyle="--",
                linewidth=1.4,
                label=f"Median after (plan cohort) ({m_after:,.0f})",
            ),
        ]
        ax.legend(handles=legend_elements, loc="upper right", fontsize=8)

        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_ylabel("Effective reads")
        ax.set_title(f"{omic}: effective reads before vs. expected after add-on")
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}")
        )
        fig.tight_layout()

        safe = omic.replace("-", "_")
        outp = plot_stem.parent / f"{plot_stem.name}_{safe}_effective_reads_before_after.png"
        fig.savefig(outp, dpi=150)
        plt.close(fig)
        written.append(outp)
    return written
