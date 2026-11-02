#!/usr/bin/env python3

import argparse
import csv
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


COLOURS = {
    "frameshift": "#B24A50",
    "stop": "#B24A50",
    "splice": "#9B6C9E",
    "missense": "#C58B45",
    "synonymous": "#5A9A91",
    "intronic": "#2563EB",
    "utr": "#668C80",
    "other": "#A0A9B2",
    "sv": "#BD7952",
    "cnv": "#69947B",
    "known": "#C82E48",
}


def consequence_class(value):
    value = (value or "").lower()

    if "frameshift" in value:
        return "frameshift"
    if "stop_gained" in value or "stop_lost" in value:
        return "stop"
    if "splice_acceptor" in value or "splice_donor" in value:
        return "splice"
    if "missense" in value:
        return "missense"
    if "splice_region" in value:
        return "splice"
    if "synonymous" in value:
        return "synonymous"
    if "utr" in value:
        return "utr"
    if "intron" in value:
        return "intronic"

    return "other"


def read_bed(path):
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.split()
                if len(fields) >= 3:
                    yield fields


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--gene_bed", required=True)
    parser.add_argument("--exon_bed", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    matches = [
        r for r in read_bed(args.gene_bed)
        if len(r) >= 4 and r[3] == args.gene
    ]

    if not matches:
        raise ValueError(f"Gene {args.gene} not found in BED")

    chrom, start, end = matches[0][:3]
    start, end = int(start), int(end)
    strand = matches[0][5] if len(matches[0]) > 5 else "."

    with open(args.csv, newline="") as handle:
        variants = list(csv.DictReader(handle))

    exons = []
    for fields in read_bed(args.exon_bed):
        if fields[0] != chrom:
            continue

        left = max(start, int(fields[1]))
        right = min(end, int(fields[2]))

        if left < right:
            exons.append((left, right))

    if not exons:
        raise ValueError(f"No exon intervals for {args.gene}")

    tracks = [
        ("small_variant", "Small variants", 3),
        ("structural_variant", "Structural variants", 2),
        ("cnv", "Copy-number variants", 1),
    ]

    fig, ax = plt.subplots(figsize=(14, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.axhspan(-0.35, 0.4, color="#F5F7F9", zorder=0)

    for y in (1, 2, 3):
        ax.axhline(y, color="#E6E9ED", linewidth=0.8, zorder=0)

    plotted_counts = {}
    consequence_counts = Counter()
    known_plotted = False

    for variant_type, label, y in tracks:
        selected = [
            v for v in variants
            if v.get("variant_type") == variant_type
        ]

        plotted = 0

        for v in selected:
            try:
                position = int(v["position"])
                variant_end = int(v.get("variant_end") or position)
            except (ValueError, TypeError, KeyError):
                continue

            left = max(start, position - 1)
            right = min(end, max(position, variant_end))

            if left >= right:
                continue

            known = v.get("known_variant", "").lower() == "yes"
            known_plotted |= known

            if variant_type == "small_variant":
                category = consequence_class(v.get("consequence"))
                consequence_counts[category] += 1

                colour = (
                    COLOURS["known"] if known
                    else COLOURS[category]
                )

                ax.vlines(
                    position, y - 0.16, y + 0.08,
                    color=colour, linewidth=1.2, alpha=0.85
                )
                ax.scatter(
                    position, y + 0.08,
                    s=100 if known else 25,
                    marker="*" if known else "o",
                    color=colour,
                    edgecolor="white",
                    linewidth=0.4,
                    zorder=4
                )

            else:
                colour = (
                    COLOURS["known"] if known
                    else COLOURS["sv" if variant_type == "structural_variant" else "cnv"]
                )

                if variant_end > position:
                    ax.hlines(
                        y, left, right,
                        color=colour, linewidth=5, zorder=3
                    )
                else:
                    ax.scatter(
                        position, y,
                        s=65,
                        marker="*" if known else "D",
                        color=colour,
                        zorder=4
                    )

            plotted += 1

        plotted_counts[variant_type] = plotted

        if plotted == 0:
            ax.text(
                (start + end) / 2, y,
                "No variants reported",
                ha="center", va="center",
                fontsize=9, color="#98A2AC"
            )

    # Gene backbone and annotated exon union.
    ax.hlines(
        0, start, end,
        color="#8995A3", linewidth=2, zorder=1
    )

    for left, right in exons:
        ax.add_patch(
            Rectangle(
                (left, -0.10), right - left, 0.20,
                facecolor="#34465B",
                edgecolor="none",
                zorder=3
            )
        )

    # Strand direction.
    direction = -1 if strand == "-" else 1
    arrow_start = start + 0.50 * (end - start)
    arrow_end = arrow_start + direction * 0.06 * (end - start)

    ax.annotate(
        "",
        xy=(arrow_end, -0.23),
        xytext=(arrow_start, -0.23),
        arrowprops=dict(
            arrowstyle="->",
            color="#566578",
            linewidth=1.5
        )
    )

    ax.set_xlim(start, end)
    ax.set_ylim(-0.5, 3.65)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels([
        f"Exon union (n={len(exons)})",
        f"CNVs (n={plotted_counts['cnv']})",
        f"SVs (n={plotted_counts['structural_variant']})",
        f"Small variants (n={plotted_counts['small_variant']})",
    ])

    ax.tick_params(axis="y", length=0, pad=12, labelsize=9)
    ax.tick_params(axis="x", labelsize=9)

    ax.ticklabel_format(
        axis="x", style="plain", useOffset=False
    )

    ax.set_xlabel(
        f"Genomic position on {chrom} (GRCh38)",
        fontsize=10, labelpad=12
    )

    ax.set_title(
        f"{args.gene} | Genomic variant landscape",
        loc="left", fontsize=15, fontweight="bold", pad=20
    )

    ax.text(
        0, 1.025,
        f"Gene interval: {chrom}:{start + 1:,}-{end:,}  |  Gene orientation: {'Reverse strand (−)' if strand == '-' else 'Forward strand (+)' if strand == '+' else 'Unknown'}",
        transform=ax.transAxes,
        fontsize=9, color="#667487"
    )

    legend_labels = {
        "frameshift": "Frameshift",
        "stop": "Stop gained/lost",
        "splice": "Splice-related",
        "missense": "Missense",
        "synonymous": "Synonymous",
        "intronic": "Intronic",
        "utr": "UTR",
        "other": "Other",
    }

    handles = [
        Line2D(
            [0], [0],
            marker="o", linestyle="none",
            markerfacecolor=COLOURS[category],
            markeredgecolor="none",
            markersize=9,
            label=f"{legend_labels[category]} ({count})"
        )
        for category, count in consequence_counts.items()
    ]

    if known_plotted:
        handles.append(
            Line2D(
                [0], [0],
                marker="*", linestyle="none",
                color=COLOURS["known"],
                markersize=11,
                label="Previously identified variant"
            )
        )

    if handles:
        fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.54, 0.035),
            ncol=len(handles),
            frameon=False,
            fontsize=10
        )

    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    ax.spines["bottom"].set_color("#AAB3BC")

    fig.subplots_adjust(
        left=0.20, right=0.97,
        top=0.82, bottom=0.27
    )

    for extension in ("png", "pdf"):
        fig.savefig(
            f"{args.output_prefix}.{extension}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white"
        )

    plt.close(fig)


if __name__ == "__main__":
    main()
