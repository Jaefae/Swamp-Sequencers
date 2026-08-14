#!/usr/bin/env python3
"""Render SwampBench CSV output into README charts.

One SVG per metric per theme, so each chart can be linked on its own. A README
is read on a white canvas *and* a near-black one, so the pair is joined with
<picture>:

    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="...-dark.svg">
      <img alt="..." src="...-light.svg">
    </picture>

Requires matplotlib (`pip install matplotlib`).

Usage:
    python bench/render_charts.py [bench/results/data.csv]
"""

import csv
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

matplotlib.rcParams["svg.fonttype"] = "path"

# For rendering in both dark and light mode. The two series colors are a
# validated categorical pair; they clear colorblind separation and 3:1 contrast
# against both canvases. Re-validate if you change them rather than eyeballing.
THEMES = {
    "light": {
        "array": "#2a78d6",
        "tree": "#eb6834",
        "grid": "#e7e6e1",
        "axis": "#c3c2b7",
        "primary": "#0b0b0b",
        "secondary": "#52514e",
        "muted": "#898781",
    },
    "dark": {
        "array": "#3987e5",
        "tree": "#d95926",
        "grid": "#272d36",
        "axis": "#444c56",
        "primary": "#ffffff",
        "secondary": "#c3c2b7",
        "muted": "#8b949e",
    },
}

# Text is embedded as outlines, so the viewer never needs these installed.
FONT_STACK = ["Segoe UI", "DejaVu Sans", "sans-serif"]

SUBTITLE = "E. coli reference genome (5M+ base pairs) · clang 22.1.4 · -O3"


def read_csv(path):
    """CSV -> {(experiment, structure, metric): {x: value}}."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            key = (row["experiment"], row["structure"], row["metric"])
            out.setdefault(key, {})[int(row["x"])] = float(row["value"])
    return out


def series(data, experiment, structure, metric, xs):
    return [data[(experiment, structure, metric)][x] for x in xs]


def draw_ratio_bracket(ax, x, lo, hi, span, t):
    """Vertical bracket spanning the two endpoints, labelled with their ratio.

    This is the comparison the chart exists to make, so it gets drawn rather
    than left for the reader to compute from two axis positions.
    """
    bx = x + 0.07 * span
    cap = 0.02 * span
    ax.plot([bx, bx], [lo, hi], color=t["muted"], linewidth=1,
            clip_on=False, zorder=5)
    for y in (lo, hi):
        ax.plot([bx - cap, bx + cap], [y, y], color=t["muted"], linewidth=1,
                clip_on=False, zorder=5)
    ax.annotate(f"{hi / lo:.1f}×", xy=(bx + cap * 1.8, (lo + hi) / 2),
                ha="left", va="center", fontsize=12, fontweight="bold",
                color=t["primary"], clip_on=False, zorder=6)


def draw_chart(fig, ax, spec, t):
    """One line chart: 2px lines, a soft wash beneath, endpoint labels only."""
    xs = spec["xs"]
    xmin, xmax = min(xs), max(xs)
    span = (xmax - xmin) or 1

    for name, ckey, vals in spec["lines"]:
        color = t[ckey]
        ax.fill_between(xs, vals, color=color, alpha=0.09, linewidth=0,
                        zorder=2)
        # A single marker on the last point anchors the value label. Marking
        # every point turns a 50-step sweep into a bead necklace.
        ax.plot(xs, vals, color=color, linewidth=2, marker="o",
                markevery=[-1], markersize=5,
                solid_capstyle="round", solid_joinstyle="round",
                label=name, clip_on=False, zorder=3)
        ax.annotate(spec["yfmt"](vals[-1]), xy=(xs[-1], vals[-1]),
                    xytext=(0, 11), textcoords="offset points",
                    ha="center", va="bottom", fontsize=10.5,
                    fontweight="bold", color=t["primary"], clip_on=False,
                    zorder=4)

    ax.set_xlabel(spec["xlabel"], color=t["secondary"], fontsize=10, labelpad=7)

    ymax = max(max(vals) for _, _, vals in spec["lines"])
    ax.set_ylim(0, ymax * 1.16)
    # Right margin holds the ratio bracket and its label.
    ax.set_xlim(xmin, xmax + 0.22 * span)

    # Let the locator choose a handful of round ticks. Ticking every measured
    # point is unreadable once the sweep has more than a few steps.
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: spec["xfmt"](v)))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: spec.get("ytickfmt", spec["yfmt"])(v)))

    lo, hi = sorted(vals[-1] for _, _, vals in spec["lines"])
    draw_ratio_bracket(ax, xmax, lo, hi, span, t)

    # Recessive chrome: hairline solid gridlines, baseline only.
    ax.grid(axis="y", color=t["grid"], linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(t["axis"])
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(colors=t["muted"], labelsize=9, length=0, pad=5)


def build_specs(data):
    scaling_xs = sorted(data[("scaling", "array", "build_ms")])

    as_int = lambda v: f"{v:,.0f}"
    as_2dp = lambda v: f"{v:.2f}"
    as_bases = lambda v: f"{v / 1_000_000:.1f}M"

    # Short patterns match so often that the timing measures building the
    # result list, not searching, and their spike flattens everything else.
    # Which lengths those are depends on the text, so derive it instead of
    # hardcoding: at large m a sampled pattern matches about once, so the
    # smallest total_hits stands in for the query count.
    hits = data[("query", "array", "total_hits")]
    all_m = sorted(data[("query", "array", "query_us")])
    baseline = min(hits[m] for m in all_m)
    query_xs = [m for m in all_m if hits[m] <= 3 * baseline]
    dropped = [m for m in all_m if m not in query_xs]

    note = None
    if dropped:
        note = ("Omits m = " + ", ".join(str(m) for m in dropped) +
                ": patterns this short match too many times, so the timing is "
                "dominated by building the result list rather than searching.")

    # Stay in MB even when the tree runs into gigabytes: in GB the array reads
    # "0.04", which throws away the only precision that side of the comparison
    # has. The bracket carries the ratio either way.
    scaled = lambda vs: [v / 1e6 for v in vs]

    return [
        {
            "slug": "memory",
            "title": "Index memory · MB",
            "xs": scaling_xs, "xfmt": as_bases,
            "xlabel": "text length n", "yfmt": as_int,
            "lines": [
                ("suffix array", "array",
                 scaled(series(data, "scaling", "array", "memory_bytes", scaling_xs))),
                ("suffix tree", "tree",
                 scaled(series(data, "scaling", "tree", "memory_bytes", scaling_xs))),
            ],
        },
        {
            "slug": "build",
            "title": "Build time · ms",
            "xs": scaling_xs, "xfmt": as_bases,
            "xlabel": "text length n", "yfmt": as_int,
            "lines": [
                ("suffix array", "array",
                 series(data, "scaling", "array", "build_ms", scaling_xs)),
                ("suffix tree", "tree",
                 series(data, "scaling", "tree", "build_ms", scaling_xs)),
            ],
        },
        {
            "slug": "query",
            "title": "Query latency · µs",
            "note": note,
            "xs": query_xs, "xfmt": as_int,
            "xlabel": "pattern length m",
            "yfmt": as_2dp, "ytickfmt": as_int,
            "lines": [
                ("suffix array", "array",
                 series(data, "query", "array", "query_us", query_xs)),
                ("suffix tree", "tree",
                 series(data, "query", "tree", "query_us", query_xs)),
            ],
        },
    ]


def render(spec, theme, out_path):
    t = THEMES[theme]
    plt.rcParams["font.family"] = FONT_STACK

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.72,
                        bottom=0.235 if spec.get("note") else 0.14)

    draw_chart(fig, ax, spec, t)

    fig.text(0.105, 0.965, spec["title"], fontsize=14, fontweight="bold",
             color=t["primary"], ha="left", va="top")
    fig.text(0.105, 0.885, SUBTITLE, fontsize=9.5, color=t["muted"],
             ha="left", va="top")
    if spec.get("note"):
        # Wrap explicitly: matplotlib's own wrap= does not respect the figure
        # margins, so a long note runs off the canvas.
        fig.text(0.105, 0.025, textwrap.fill(spec["note"], 84), fontsize=8.5,
                 color=t["muted"], ha="left", va="bottom", linespacing=1.5)

    # Legend always present for two or more series, so identity never rests on
    # color alone.
    legend = fig.legend(loc="upper right", bbox_to_anchor=(0.985, 0.99),
                        ncol=2, frameon=False, fontsize=10, handlelength=1.6,
                        columnspacing=1.6)
    for txt in legend.get_texts():
        txt.set_color(t["secondary"])

    fig.savefig(out_path, format="svg", transparent=True)
    plt.close(fig)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        root, "bench", "results", "ecoli-50.csv")
    data = read_csv(src)

    outdir = os.path.dirname(src)
    for spec in build_specs(data):
        for theme in THEMES:
            path = os.path.join(outdir, f"{spec['slug']}-{theme}.svg")
            render(spec, theme, path)
            print(f"wrote {os.path.relpath(path, root)}")


if __name__ == "__main__":
    main()
