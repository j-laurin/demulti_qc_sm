#!/usr/bin/env python3

"""Extract cutadapt statistics from log files and generate plots. Snakemake version."""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sns

TOTAL_RE = re.compile(r"Total reads processed:\s*([\d,]+)")
WRITTEN_RE = re.compile(r"Reads written \(passing filters\):\s*([\d,]+)")
ADAPTER_BLOCK_RE = re.compile(
    r"=== Adapter\s+'?([^=']+?)'?\s*===.*?Trimmed:\s*([\d,]+)\s*times",
    re.S,
)


def parse_cutadapt_summary(path: Path) -> tuple[int, int]:
    """Return (total_reads_processed, reads_written) from a cutadapt log."""
    text = path.read_text()
    total_m = TOTAL_RE.search(text)
    written_m = WRITTEN_RE.search(text)
    if not total_m or not written_m:
        raise ValueError(f"Could not parse summary stats from {path}")
    total = int(total_m.group(1).replace(",", ""))
    written = int(written_m.group(1).replace(",", ""))
    return total, written


def parse_adapter_blocks(path: Path) -> dict:
    """Return {adapter_name: trimmed_count} for each '=== Adapter ===' block."""
    text = path.read_text()
    blocks = ADAPTER_BLOCK_RE.findall(text)
    if not blocks:
        raise ValueError(f"Could not find any '=== Adapter ===' blocks in {path}")
    return {name.strip(): int(count.replace(",", "")) for name, count in blocks}


def collect_sample_stats(sample: str, step1_path: Path, demux_path: Path,
                          step3_dir: Path, n_libraries: int) -> dict:
    for p in (step1_path, demux_path, step3_dir):
        if not p.exists():
            raise FileNotFoundError(f"Expected path not found: {p}")

    raw_total, step1_written = parse_cutadapt_summary(step1_path)

    demux_total, _ = parse_cutadapt_summary(demux_path)
    adapter_counts = parse_adapter_blocks(demux_path)
    demux_assigned = sum(adapter_counts.values())
    demux_unassigned = demux_total - demux_assigned

    step3_logs = sorted(Path(step3_dir).glob("*.log"))
    if len(step3_logs) != n_libraries:
        raise ValueError(
            f"Expected {n_libraries} step3 logs for sample {sample} in "
            f"{step3_dir}, found {len(step3_logs)}"
        )

    library_written = {}
    for log in step3_logs:
        _, written_i = parse_cutadapt_summary(log)
        library_written[log.stem] = written_i
    after_step3 = sum(library_written.values())

    return {
        "sample": sample,
        "raw_reads": raw_total,
        "after_step1": step1_written,
        "demux_input": demux_total,
        "after_demux": demux_assigned,
        "demux_unassigned": demux_unassigned,
        "adapter_counts": adapter_counts,
        "after_step3": after_step3,
        "library_written": library_written,
    }


def make_sample_page(pdf: PdfPages, stats: dict):
    sample = stats["sample"]

    cum_df = pd.DataFrame({
        "stage": ["raw", "step1\n(adaptor trim)", "step2\n(demux, assigned)",
                  "step3\n(length filter)"],
        "reads": [stats["raw_reads"], stats["after_step1"],
                  stats["after_demux"], stats["after_step3"]],
    })

    demux_items = list(stats["adapter_counts"].items()) + \
        [("unassigned", stats["demux_unassigned"])]
    demux_df = pd.DataFrame(demux_items, columns=["barcode", "reads"])

    table_df = cum_df.copy()
    table_df["% of raw"] = [f"{r / stats['raw_reads'] * 100:.1f}%"
                             for r in cum_df["reads"]]
    table_df["reads"] = table_df["reads"].map(lambda x: f"{x:,}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.27, 11.69))  # A4 portrait
    fig.suptitle(f"Sample: {sample}", fontsize=14, fontweight="bold")

    sns.barplot(data=cum_df, x="stage", y="reads", ax=ax1, color="#4C72B0")
    ax1.set_title("Remaining reads per step (cumulative)")
    ax1.set_xlabel("")
    ax1.set_ylabel("reads")
    for p in ax1.patches:
        ax1.annotate(f"{int(p.get_height()):,}",
                      (p.get_x() + p.get_width() / 2, p.get_height()),
                      ha="center", va="bottom", fontsize=8)

    sns.barplot(data=demux_df, x="barcode", y="reads", ax=ax2, color="#DD8452")
    ax2.set_title("Demultiplexing outcome (step 2)")
    ax2.set_xlabel("")
    ax2.set_ylabel("reads")
    ax2.tick_params(axis="x", rotation=45)
    for p in ax2.patches:
        ax2.annotate(f"{int(p.get_height()):,}",
                      (p.get_x() + p.get_width() / 2, p.get_height()),
                      ha="center", va="bottom", fontsize=8)

    plt.tight_layout(rect=[0, 0.28, 1, 0.95])

    ax_table = fig.add_axes([0.08, 0.02, 0.84, 0.20])
    ax_table.axis("off")
    tbl = ax_table.table(
        cellText=table_df.values,
        colLabels=table_df.columns,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.3)

    pdf.savefig(fig)
    plt.close(fig)


def main():
    samples = snakemake.params.samples
    step1_paths = snakemake.input.step1
    demux_paths = snakemake.input.demux
    step3_dirs = snakemake.params.step3_dirs
    n_libraries = snakemake.params.get("n_libraries", 8)

    pdf_out = Path(snakemake.output.pdf)
    csv_out = Path(snakemake.output.csv)

    all_rows = []
    with PdfPages(pdf_out) as pdf:
        for sample, step1_path, demux_path, step3_dir in zip(
            samples, step1_paths, demux_paths, step3_dirs
        ):
            print(f"Processing {sample}...")
            stats = collect_sample_stats(
                sample, Path(step1_path), Path(demux_path), Path(step3_dir), n_libraries
            )
            make_sample_page(pdf, stats)

            row = {k: v for k, v in stats.items()
                   if k not in ("adapter_counts", "library_written")}
            row.update({f"adapter_{k}": v for k, v in stats["adapter_counts"].items()})
            row.update({f"lib_{k}": v for k, v in stats["library_written"].items()})
            all_rows.append(row)

    pd.DataFrame(all_rows).to_csv(csv_out, index=False)
    print(f"Done! Saved outputs to '{pdf_out}' and '{csv_out}'.")


if __name__ == "__main__":
    main()