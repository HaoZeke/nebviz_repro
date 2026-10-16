#!/usr/bin/env python3
"""Force-versus-step line for the blog cycloaddition, plus readcon helpers."""

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

# Ruhi palette (chemparseplot RUHI_COLORS).
COLORS = {
    "ink": "#1f2933",
    "teal": "#004D40",
    "coral": "#FF655D",
}


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
            "axes.grid": False,
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


def _plateau_start(barriers: np.ndarray) -> int | None:
    """First index whose tail stays within 0.04 eV (the barrier has flattened)."""
    if len(barriers) < 4:
        return None
    for idx in range(len(barriers)):
        if float(np.ptp(barriers[idx:])) <= 0.04 and len(barriers) - idx >= 3:
            return idx
    return None


def plot_force_vs_step(neb_root: Path, blog_system: str, output: Path) -> None:
    """One line plot: barrier and max parallel force against optimization step.

    The cycloaddition barrier settles while max |F∥| stays near 2.7–2.9 eV/Å.
    Both series share the step axis so that split is visible without a second panel.
    """
    grouped = load_iterations(neb_root / blog_system)
    iterations = np.asarray(sorted(grouped.keys()), dtype=float)
    barriers = np.asarray(
        [max(float(row["energy"]) for row in grouped[int(it)]) for it in iterations]
    )
    max_force = np.asarray(
        [max(abs(float(row["f_para"])) for row in grouped[int(it)]) for it in iterations]
    )
    label = SYSTEM_LABELS.get(blog_system, blog_system)

    fig, ax = plt.subplots(figsize=(6.8, 3.9), layout="constrained")
    ax.plot(iterations, barriers, color=COLORS["teal"], lw=2.3, label="Barrier")
    ax.set_xlabel("NEB optimization step")
    ax.set_ylabel("Barrier / eV", color=COLORS["teal"])
    ax.tick_params(axis="y", labelcolor=COLORS["teal"])
    b_lo, b_hi = float(barriers.min()), float(barriers.max())
    b_pad = max(0.08 * (b_hi - b_lo), 0.02)
    ax.set_ylim(b_lo - b_pad, b_hi + b_pad)

    ax_f = ax.twinx()
    ax_f.plot(iterations, max_force, color=COLORS["coral"], lw=2.3, label="Max |F∥|")
    ax_f.set_ylabel("Max |F∥| / eV Å$^{-1}$", color=COLORS["coral"])
    ax_f.tick_params(axis="y", labelcolor=COLORS["coral"])
    f_lo, f_hi = float(max_force.min()), float(max_force.max())
    f_pad = max(0.12 * (f_hi - f_lo), 0.04)
    ax_f.set_ylim(f_lo - f_pad, f_hi + 1.6 * f_pad)
    ax_f.grid(False)
    ax_f.spines["top"].set_visible(False)

    ax.set_xlim(float(iterations[0]) - 0.4, float(iterations[-1]) + 0.4)
    step0 = int(iterations[0])
    step1 = int(iterations[-1])
    stride = 4 if step1 - step0 > 16 else 2
    ax.set_xticks(list(range(step0, step1 + 1, stride)))
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%d"))
    ax.set_title(label)
    ax.grid(False)

    handles = [
        mpl.lines.Line2D([0], [0], color=COLORS["teal"], lw=2.3, label="Barrier"),
        mpl.lines.Line2D([0], [0], color=COLORS["coral"], lw=2.3, label="Max |F∥|"),
    ]
    # Both series leave the upper-right empty after the first few steps.
    ax.legend(handles=handles, loc="upper right", frameon=False, borderaxespad=0.15)

    start = _plateau_start(barriers)
    if start is not None:
        late_force = max_force[start:]
        lo, hi = float(late_force.min()), float(late_force.max())
        if hi >= 0.5 and (hi - lo) <= 1.0:
            step = int(iterations[start])
            if 2.6 <= lo <= 2.85 and 2.85 <= hi <= 3.05:
                force_note = "~2.7–2.9 eV/Å"
            else:
                force_note = f"{lo:.2f}–{hi:.2f} eV/Å"
            # Sit in the open band above the plateau, clear of the right-hand ticks.
            ax.annotate(
                f"Barrier flat from step {step}\nMax |F∥| stays {force_note}",
                xy=(0.50, 0.70),
                xycoords="axes fraction",
                ha="center",
                va="top",
                fontsize=10,
                color=COLORS["ink"],
                zorder=5,
            )
    save(fig, output)


def diagnostics(args: argparse.Namespace) -> None:
    set_style()
    plot_force_vs_step(Path(args.neb_root), args.blog_system, Path(args.force))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    summary = subparsers.add_parser("readcon-summary")
    summary.add_argument("--output", required=True)
    summary.add_argument("cons", nargs="+")

    diag = subparsers.add_parser("diagnostics")
    diag.add_argument("--neb-root", required=True)
    diag.add_argument("--blog-system", required=True)
    diag.add_argument("--force", required=True)

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
