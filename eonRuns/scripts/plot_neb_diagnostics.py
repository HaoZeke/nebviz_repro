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

COLORS = {
    "ink": "#1f2933",
    "muted": "#657786",
    "teal": "#007c78",
    "coral": "#d95f4f",
    "gold": "#b7791f",
    "plum": "#6b4c9a",
}


def set_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#d7dde2",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.7,
            "legend.frameon": False,
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
                    "image": int(parts[0]),
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
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_pipeline(path: Path) -> None:
    set_style()
    fig, ax = plt.subplots(figsize=(12.2, 2.45))
    ax.set_axis_off()
    labels = [
        ("Inputs", "reactant\nproduct\nsaddle"),
        ("eOn", "PET-MAD\nreadcon v2 con\nNEB movies"),
        ("rgpycrumbs", "path profile\nimage profile\nRMSD profile"),
        ("2D RMSD", "sampled images\nGP surface\nuncertainty"),
        ("Diagnostics", "barrier trace\nforce heatmap\nspacing"),
    ]
    xs = np.linspace(0.08, 0.92, len(labels))
    for idx, (x, (title, body)) in enumerate(zip(xs, labels)):
        rect = mpl.patches.FancyBboxPatch(
            (x - 0.085, 0.33),
            0.17,
            0.42,
            boxstyle="round,pad=0.018,rounding_size=0.015",
            linewidth=1.2,
            edgecolor=COLORS["teal"] if idx % 2 == 0 else COLORS["plum"],
            facecolor="#f7fafc",
        )
        ax.add_patch(rect)
        ax.text(x, 0.62, title, ha="center", va="center", weight="bold", color=COLORS["ink"])
        ax.text(x, 0.47, body, ha="center", va="center", color=COLORS["muted"], linespacing=1.25)
        if idx < len(labels) - 1:
            ax.annotate(
                "",
                xy=(xs[idx + 1] - 0.105, 0.54),
                xytext=(x + 0.105, 0.54),
                arrowprops={"arrowstyle": "->", "lw": 1.6, "color": COLORS["gold"]},
            )
    ax.text(0.5, 0.12, "Snakemake target graph for the blog figures", ha="center", color=COLORS["ink"])
    save(fig, path)


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
    iterations = np.asarray(list(grouped.keys()))
    barriers = np.asarray([max(float(row["energy"]) for row in grouped[it]) for it in iterations])
    max_force = np.asarray([max(abs(float(row["f_para"])) for row in grouped[it]) for it in iterations])

    fig, ax1 = plt.subplots(figsize=(7.1, 4.1))
    ax1.plot(iterations, barriers, color=COLORS["teal"], lw=2.0, label="barrier")
    ax1.set_xlabel("NEB optimization step")
    ax1.set_ylabel("Barrier / eV", color=COLORS["teal"])
    ax1.tick_params(axis="y", labelcolor=COLORS["teal"])
    ax2 = ax1.twinx()
    ax2.plot(iterations, max_force, color=COLORS["coral"], lw=2.0, label="max |parallel force|")
    ax2.set_ylabel("Max |parallel force| / eV A$^{-1}$", color=COLORS["coral"])
    ax2.tick_params(axis="y", labelcolor=COLORS["coral"])
    ax1.set_title("Barrier and force convergence")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="upper right")
    save(fig, output)


def plot_force_heatmap(neb_root: Path, blog_system: str, output: Path) -> None:
    grouped = load_iterations(neb_root / blog_system)
    iterations = list(grouped)
    images = [int(row["image"]) for row in grouped[iterations[-1]]]
    heat = np.zeros((len(iterations), len(images)))
    for i, iteration in enumerate(iterations):
        force_by_image = {int(row["image"]): abs(float(row["f_para"])) for row in grouped[iteration]}
        for j, image in enumerate(images):
            heat[i, j] = force_by_image.get(image, np.nan)

    fig, ax = plt.subplots(figsize=(7.1, 4.2))
    mesh = ax.imshow(
        heat,
        aspect="auto",
        origin="lower",
        cmap="magma",
        extent=[min(images), max(images), min(iterations), max(iterations)],
    )
    fig.colorbar(mesh, ax=ax, label="|parallel force| / eV A$^{-1}$")
    ax.set_xlabel("Image index")
    ax.set_ylabel("NEB optimization step")
    ax.set_title("Where the band still moves")
    save(fig, output)


def plot_sampling_density(neb_root: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.1))
    colors = [COLORS["teal"], COLORS["gold"], COLORS["plum"]]
    for idx, system in enumerate(SYSTEM_LABELS):
        grouped = load_iterations(neb_root / system)
        final_rows = grouped[max(grouped)]
        rxn = np.asarray([float(row["rxn_coord"]) for row in final_rows])
        spacing = np.diff(rxn)
        ax.plot(
            np.arange(len(spacing)) + idx * 0.08,
            spacing,
            marker="o",
            lw=1.7,
            ms=4.0,
            color=colors[idx],
            label=SYSTEM_LABELS[system],
        )
    ax.set_xlabel("Segment index")
    ax.set_ylabel("Path-length spacing / A")
    ax.set_title("Final-band image spacing")
    ax.legend(loc="best")
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

    pipeline = subparsers.add_parser("pipeline")
    pipeline.add_argument("--output", required=True)

    summary = subparsers.add_parser("readcon-summary")
    summary.add_argument("--output", required=True)
    summary.add_argument("cons", nargs="+")

    diag = subparsers.add_parser("diagnostics")
    diag.add_argument("--neb-root", required=True)
    diag.add_argument("--blog-system", required=True)
    diag.add_argument("--convergence", required=True)
    diag.add_argument("--force-heatmap", required=True)
    diag.add_argument("--sampling-density", required=True)

    args = parser.parse_args()
    if args.command == "pipeline":
        plot_pipeline(Path(args.output))
    elif args.command == "readcon-summary":
        write_readcon_summary([Path(path) for path in args.cons], Path(args.output))
    elif args.command == "diagnostics":
        diagnostics(args)


if __name__ == "__main__":
    main()
