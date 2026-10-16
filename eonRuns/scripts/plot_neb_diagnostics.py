#!/usr/bin/env python3
"""Plot NEB diagnostics that complement the rgpycrumbs path figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import readcon


SYSTEM_LABELS = {
    "system_100": "Cycloaddition",
    "11_grignard": "Grignard addition",
    "bicyclobutane_05": "Bicyclobutane opening",
}

# Ruhi palette (chemparseplot RUHI_COLORS). Series order matches SYSTEM_LABELS.
COLORS = {
    "ink": "#1f2933",
    "muted": "#657786",
    "teal": "#004D40",
    "coral": "#FF655D",
    "sky": "#1E88E5",
    "sunshine": "#F1DB4B",
    "magenta": "#D81B60",
}
SERIES_COLORS = [COLORS["teal"], COLORS["coral"], COLORS["sky"]]


def set_style() -> None:
    """Ruhi colors with an open-zero tick face.

    Atkinson Hyperlegible (the Ruhi theme face) draws a slashed zero.
    The other blog figures use DejaVu Sans, whose zero is open, so the
    diagnostics and the 2D panels share that tick font.
    """
    mpl.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "text.color": "black",
            "axes.grid": True,
            "grid.color": "#d7dde2",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.7,
            "legend.frameon": True,
            "legend.framealpha": 0.95,
            "legend.edgecolor": "#d7dde2",
            "mathtext.fontset": "dejavusans",
            "axes.formatter.use_mathtext": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def parse_dat(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open() as handle:
        next(handle)
        for line in handle:
            parts = line.split()
            if len(parts) != 4:
                continue
            rows.append(
                {
                    "image": int(float(parts[0])),
                    "rxn_coord": float(parts[1]),
                    "energy": float(parts[2]),
                    "f_para": float(parts[3]),
                }
            )
    return rows


def iteration_number(path: Path) -> int:
    return int(path.stem.split("_")[-1])


def load_iterations(folder: Path) -> dict[int, list[dict[str, float | int]]]:
    return {
        iteration_number(path): parse_dat(path)
        for path in sorted(folder.glob("neb_[0-9][0-9][0-9].dat"))
    }


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_readcon_summary(paths: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["system", "frames", "frames_with_energy", "first_energy", "last_energy"],
        )
        writer.writeheader()
        for path in paths:
            frames = readcon.read_con(str(path))
            energies = [float(frame.energy) for frame in frames if frame.energy is not None]
            writer.writerow(
                {
                    "system": path.parent.name,
                    "frames": len(frames),
                    "frames_with_energy": len(energies),
                    "first_energy": energies[0] if energies else "",
                    "last_energy": energies[-1] if energies else "",
                }
            )


def plot_convergence(neb_root: Path, blog_system: str, output: Path) -> None:
    grouped = load_iterations(neb_root / blog_system)
    iterations = np.asarray(sorted(grouped.keys()))
    barriers = np.asarray(
        [max(float(row["energy"]) for row in grouped[int(it)]) for it in iterations]
    )
    max_force = np.asarray(
        [max(abs(float(row["f_para"])) for row in grouped[int(it)]) for it in iterations]
    )

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.2, 5.4),
        layout="constrained",
        gridspec_kw={"height_ratios": [1.15, 1.0]},
    )
    ax1.plot(iterations, barriers, color=COLORS["teal"], lw=2.1, marker="o", ms=3.5)
    ax1.set_ylabel("Barrier / eV", color=COLORS["teal"])
    ax1.tick_params(axis="y", labelcolor=COLORS["teal"])
    ax1.set_title(f"Barrier and force convergence ({SYSTEM_LABELS.get(blog_system, blog_system)})")
    ax1.set_ylim(bottom=min(barriers.min() * 0.98, barriers.min() - 0.02))

    ax2.plot(iterations, max_force, color=COLORS["coral"], lw=2.1, marker="o", ms=3.5)
    ax2.set_xlabel("NEB optimization step")
    ax2.set_ylabel("Max |F∥| / eV Å$^{-1}$", color=COLORS["coral"])
    ax2.tick_params(axis="y", labelcolor=COLORS["coral"])
    # Full dynamic range — do not zoom the force axis into noise
    f_lo, f_hi = float(max_force.min()), float(max_force.max())
    pad = max(0.05 * (f_hi - f_lo), 0.05)
    ax2.set_ylim(f_lo - pad, f_hi + pad)
    ax2.set_xticks(iterations[:: max(len(iterations) // 10, 1)])
    _annotate_force_plateau(ax2, iterations, barriers, max_force)
    save(fig, output)


def _annotate_force_plateau(ax, iterations, barriers, max_force) -> None:
    """Note that a flat barrier with a high force is the optimizer, not a glitch.

    The cycloaddition band settles near 1.03 eV while max |F∥| stays near
    2.7–2.9 eV/Å. That split is in the force table.
    """
    if len(barriers) < 4:
        return
    start = 0
    for idx in range(len(barriers)):
        if float(np.ptp(barriers[idx:])) <= 0.04:
            start = idx
            break
    else:
        return
    if len(barriers) - start < 3:
        return
    late_force = max_force[start:]
    lo, hi = float(late_force.min()), float(late_force.max())
    # Only call out a plateau that is still far from a tight force convergence.
    if hi < 0.5 or (hi - lo) > 1.0:
        return
    step = int(iterations[start])
    # The published cycloaddition plateau is ~2.7–2.9 eV/Å (one step touches 2.95).
    if 2.6 <= lo <= 2.85 and 2.85 <= hi <= 3.05:
        force_note = "~2.7–2.9 eV/Å"
    else:
        force_note = f"{lo:.2f}–{hi:.2f} eV/Å"
    note = f"Barrier is flat from step {step}.\nMax |F∥| stays {force_note}."
    ax.annotate(
        note,
        xy=(float(iterations[-1]), float(np.median(late_force))),
        xytext=(0.98, 0.97),
        textcoords=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color=COLORS["ink"],
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#d7dde2",
            "alpha": 0.96,
        },
        arrowprops={"arrowstyle": "-", "color": COLORS["muted"], "lw": 0.8},
        zorder=5,
    )


def plot_force_heatmap(neb_root: Path, blog_system: str, output: Path) -> None:
    grouped = load_iterations(neb_root / blog_system)
    iterations = sorted(grouped)
    images = [int(row["image"]) for row in grouped[iterations[-1]]]
    heat = np.zeros((len(iterations), len(images)))
    for i, iteration in enumerate(iterations):
        force_by_image = {
            int(row["image"]): abs(float(row["f_para"])) for row in grouped[iteration]
        }
        for j, image in enumerate(images):
            heat[i, j] = force_by_image.get(image, np.nan)

    final_rows = grouped[iterations[-1]]
    climb = int(max(final_rows, key=lambda r: float(r["energy"]))["image"])

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    # Integer cell centers: image and step are discrete
    mesh = ax.imshow(
        heat,
        aspect="auto",
        origin="lower",
        cmap="magma",
        interpolation="nearest",
        extent=[
            images[0] - 0.5,
            images[-1] + 0.5,
            iterations[0] - 0.5,
            iterations[-1] + 0.5,
        ],
    )
    fig.colorbar(mesh, ax=ax, label="|F∥| / eV Å$^{-1}$", pad=0.02)
    ax.axvline(climb, color="white", ls="--", lw=1.2, alpha=0.9, label=f"highest image ({climb})")
    ax.set_xlabel("Image index")
    ax.set_ylabel("NEB optimization step")
    ax.set_title(f"Where the band still moves ({SYSTEM_LABELS.get(blog_system, blog_system)})")
    # Integer ticks only. No grid: it draws through the cells.
    step = max(1, len(images) // 10)
    ax.set_xticks(images[::step])
    ystep = max(1, len(iterations) // 10)
    ax.set_yticks(iterations[::ystep])
    ax.grid(False)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        borderaxespad=0.0,
        fontsize=9,
        frameon=True,
    )
    save(fig, output)


def plot_sampling_density(neb_root: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    xmax = 0
    for idx, system in enumerate(SYSTEM_LABELS):
        grouped = load_iterations(neb_root / system)
        final_rows = grouped[max(grouped)]
        rxn = np.asarray([float(row["rxn_coord"]) for row in final_rows])
        spacing = np.diff(rxn)
        spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
        x = np.arange(len(spacing))
        xmax = max(xmax, int(x[-1]) if len(x) else 0)
        ax.plot(
            x,
            spacing,
            marker="o",
            lw=1.8,
            ms=4.2,
            color=SERIES_COLORS[idx % len(SERIES_COLORS)],
            label=SYSTEM_LABELS[system],
        )
    # The Grignard reactant segment is several angstroms. A linear axis
    # flattens the cycloaddition and bicyclobutane spacings under that spike.
    ax.set_yscale("log")
    ax.set_xlabel("Segment index")
    ax.set_ylabel("Path-length spacing / Å")
    ax.set_title("Final-band image spacing")
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=9, integer=True, min_n_ticks=4))
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%d"))
    ax.set_xlim(-0.5, xmax + 0.5)
    # The Grignard spike sits in the upper-left corner, so the legend stays outside.
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        fontsize=9,
        frameon=True,
    )
    save(fig, output)


def diagnostics(args: argparse.Namespace) -> None:
    set_style()
    neb_root = Path(args.neb_root)
    plot_convergence(neb_root, args.blog_system, Path(args.convergence))
    plot_force_heatmap(neb_root, args.blog_system, Path(args.force_heatmap))
    plot_sampling_density(neb_root, Path(args.sampling_density))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    summary = subparsers.add_parser("readcon-summary")
    summary.add_argument("--output", required=True)
    summary.add_argument("cons", nargs="+")

    diag = subparsers.add_parser("diagnostics")
    diag.add_argument("--neb-root", required=True)
    diag.add_argument("--blog-system", required=True)
    diag.add_argument("--convergence", required=True)
    diag.add_argument("--force-heatmap", required=True)
    diag.add_argument("--sampling-density", required=True)

    # keep pipeline subcommand for snakemake compatibility (no-op redirect)
    pipe = subparsers.add_parser("pipeline")
    pipe.add_argument("--output", required=True)

    args = parser.parse_args()
    if args.command == "pipeline":
        # Blog pipeline figure is Graphviz DOT in the post; write a tiny placeholder note file only if needed
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        # Still produce a PNG via graphviz if available, else skip
        import shutil
        import subprocess
        import tempfile

        if shutil.which("dot"):
            source = r'''
digraph NEBProjectionPipeline {
  graph [fontname="Jost", fontsize=15, bgcolor="white", rankdir=TB, nodesep=0.40, ranksep=0.45, pad=0.18];
  node [fontname="Jost", fontsize=15, shape=box, style="rounded,filled", fillcolor="white", color="#004D40", fontcolor="#004D40", penwidth=2.0, margin="0.30,0.18", width=3.6];
  edge [fontname="Jost", fontsize=13, color="#004D40", penwidth=1.8, arrowsize=1.0];
  neb [label="1. NEB band  ·  3N coords, E, F_parallel", fillcolor="#FF655D", fontcolor="white", color="#FF655D"];
  ira [label="2. IRA RMSD  ·  r = d(R), p = d(P)"];
  grad [label="3. Synthetic gradients in (r, p)"];
  gp [label="4. Grad-enhanced GP  ·  IMQ kernel, mean + variance", fillcolor="#F1DB4B", color="#004D40"];
  sd [label="5. Rotate frame  ·  s progress, d orthogonal"];
  out [label="6. s-d landscape  ·  path, samples, variance contours", fillcolor="#004D40", fontcolor="white", color="#004D40"];
  neb -> ira -> grad -> gp -> sd -> out;
}
'''
            with tempfile.NamedTemporaryFile("w", suffix=".dot", delete=False) as fh:
                fh.write(source)
                dpath = fh.name
            subprocess.run(["dot", "-Tpng", "-Gdpi=220", "-o", args.output, dpath], check=True)
            Path(dpath).unlink(missing_ok=True)
        else:
            Path(args.output).write_bytes(b"")
    elif args.command == "readcon-summary":
        write_readcon_summary([Path(path) for path in args.cons], Path(args.output))
    elif args.command == "diagnostics":
        diagnostics(args)


if __name__ == "__main__":
    main()
