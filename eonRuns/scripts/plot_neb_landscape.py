#!/usr/bin/env python3
"""Blog 2D NEB landscapes.

Wraps ``rgpycrumbs`` ``plt-neb`` and then:

- leaves the GP surface as fitted (no white wash over unsampled regions)
- keeps σ² labels fully inside the axes
- moves the on-plot SP / reference-saddle labels off the star markers and the path
- uses DejaVu Sans for ticks so the zero is open
  (Atkinson Hyperlegible, the Ruhi face, draws a slashed zero)
- renders every structure strip with the same xyzrender ``paton`` preset

The s–d window stays the library square (1 Å of progress matches 1 Å of
deviation). Grignard and bicyclobutane spans differ by more than 2×, so the
shared map size is the matching scale; stretching either window to a common
Å-per-inch would crop the longer path or pad the shorter one with empty GP.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Bbox

log = logging.getLogger("plot_neb_landscape")

_ON_PLOT_LABELS = {"R", "SP", "P", "Birkholz", "mlFSM", "ORCA"}
_NUDGE_LABELS = {"SP", "Birkholz", "mlFSM", "ORCA"}
_ORIG_SAVEFIG = mpl.figure.Figure.savefig
_POLISHED: set[int] = set()


def _use_open_zero_font() -> None:
    """Tick labels use DejaVu Sans (open zero), not a slashed-zero face."""
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["mathtext.fontset"] = "dejavusans"
    mpl.rcParams["axes.formatter.use_mathtext"] = False


def _install_font_hook() -> None:
    import chemparseplot.plot.theme as theme_mod
    import rgpycrumbs.eon.plt_neb as plt_neb

    original = theme_mod.setup_global_theme

    def setup_global_theme(theme):
        original(theme)
        _use_open_zero_font()

    theme_mod.setup_global_theme = setup_global_theme
    plt_neb.setup_global_theme = setup_global_theme


def _is_sigma_label(text: str) -> bool:
    compact = text.replace(" ", "").lower()
    return "σ" in text or "sigma" in compact or r"\sigma" in text


def _pt(fig) -> float:
    """Pixels per typographic point at this figure's dpi."""
    return fig.dpi / 72.0


def _contain_sigma_labels(fig, ax) -> None:
    """Move σ² labels fully inside the axes, or drop ones that cannot fit."""
    renderer = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(renderer)
    margin = 8.0 * _pt(fig)
    inner = Bbox.from_extents(
        ax_box.x0 + margin,
        ax_box.y0 + margin,
        ax_box.x1 - margin,
        ax_box.y1 - margin,
    )
    for text in list(ax.texts):
        if not _is_sigma_label(text.get_text()):
            continue
        text.set_clip_on(False)
        for _ in range(3):
            box = text.get_window_extent(renderer).padded(1.5 * _pt(fig))
            shift_x = 0.0
            shift_y = 0.0
            if box.x0 < inner.x0:
                shift_x = inner.x0 - box.x0
            elif box.x1 > inner.x1:
                shift_x = inner.x1 - box.x1
            if box.y0 < inner.y0:
                shift_y = inner.y0 - box.y0
            elif box.y1 > inner.y1:
                shift_y = inner.y1 - box.y1
            if not shift_x and not shift_y:
                break
            display = ax.transData.transform(text.get_position())
            moved = ax.transData.inverted().transform(
                (display[0] + shift_x, display[1] + shift_y)
            )
            text.set_position((float(moved[0]), float(moved[1])))
        box = text.get_window_extent(renderer)
        if box.x0 < ax_box.x0 or box.x1 > ax_box.x1 or box.y0 < ax_box.y0 or box.y1 > ax_box.y1:
            text.set_visible(False)


def _star_markers(fig, ax) -> list[tuple[str, float, float, float]]:
    """Saddle stars as (label, display x, display y, radius in pixels)."""
    px = _pt(fig)
    markers: list[tuple[str, float, float, float]] = []
    for coll in ax.collections:
        if not isinstance(coll, PathCollection):
            continue
        offsets = np.asarray(coll.get_offsets(), dtype=float)
        if offsets.ndim != 2 or offsets.shape[1] != 2 or not (1 <= len(offsets) <= 4):
            continue
        label = coll.get_label() or ""
        if label.startswith("_"):
            continue
        sizes = np.atleast_1d(np.asarray(coll.get_sizes(), dtype=float))
        area = float(sizes[0]) if len(sizes) else 280.0
        # Matplotlib scales a radius-0.5 star by sqrt(area). Pad past the points.
        radius = (0.62 * np.sqrt(max(area, 1.0)) + 3.0) * px
        for x, y in offsets:
            if not np.isfinite(x) or not np.isfinite(y):
                continue
            display = ax.transData.transform((float(x), float(y)))
            markers.append((label, float(display[0]), float(display[1]), radius))
    for line in ax.lines:
        if line.get_marker() != "*":
            continue
        label = line.get_label() or ""
        markersize = float(line.get_markersize() or 12.0)
        radius = (0.55 * markersize + 3.0) * px
        for x, y in zip(line.get_xdata(), line.get_ydata(), strict=False):
            display = ax.transData.transform((float(x), float(y)))
            markers.append((label, float(display[0]), float(display[1]), radius))
    return markers


def _boxes_overlap(left: Bbox, right: Bbox, pad: float = 0.0) -> bool:
    return not (
        left.x1 < right.x0 - pad
        or left.x0 > right.x1 + pad
        or left.y1 < right.y0 - pad
        or left.y0 > right.y1 + pad
    )


def _clearance(box: Bbox, cx: float, cy: float, radius: float) -> float:
    """Gap between a text box and a circular marker. Negative when they meet."""
    nearest_x = min(max(cx, box.x0), box.x1)
    nearest_y = min(max(cy, box.y0), box.y1)
    return float(np.hypot(nearest_x - cx, nearest_y - cy) - radius)


# Offsets in points from the marker the label belongs to. The text anchor is
# the bottom-center (va=bottom, ha=center), so a negative y has to clear the
# whole word, not just the star. Kept in points so 300 dpi matches a draft.
_LABEL_SHIFT_PT = (
    (26, -34),
    (16, -38),
    (36, -26),
    (-30, -34),
    (22, 16),
    (-26, 16),
    (30, 8),
    (-32, 6),
    (0, 22),
    (18, 26),
    (-20, 26),
    (34, -8),
    (-36, -6),
    (8, 34),
    (-12, 34),
    (0, -40),
    (42, 18),
    (-44, 18),
)


def _path_samples(fig, ax) -> list[tuple[float, float, float]]:
    """Thick path stroke as (display x, display y, radius in pixels)."""
    px = _pt(fig)
    samples: list[tuple[float, float, float]] = []
    for coll in ax.collections:
        if not isinstance(coll, LineCollection):
            continue
        widths = np.atleast_1d(np.asarray(coll.get_linewidth(), dtype=float))
        width = float(np.nanmax(widths)) if len(widths) else 0.0
        if not np.isfinite(width) or width < 2.0:
            continue
        radius = (0.5 * max(width, 6.5) + 7.0) * px
        for seg in coll.get_segments():
            seg = np.asarray(seg, dtype=float)
            if seg.ndim != 2 or len(seg) < 2:
                continue
            for start, end in zip(seg[:-1], seg[1:], strict=False):
                length = float(np.hypot(*(end - start)))
                steps = max(int(length / 0.05), 1)
                for frac in np.linspace(0.0, 1.0, steps + 1):
                    data = start * (1.0 - frac) + end * frac
                    if not np.isfinite(data).all():
                        continue
                    display = ax.transData.transform(data)
                    samples.append((float(display[0]), float(display[1]), radius))
    return samples


def _separate_marker_labels(fig, ax) -> None:
    """Put SP and reference labels clear of both saddle stars."""
    for patch in list(ax.patches):
        if isinstance(patch, FancyArrowPatch):
            patch.remove()
    renderer = fig.canvas.get_renderer()
    stars = _star_markers(fig, ax)
    path_samples = _path_samples(fig, ax)
    labels = [text for text in ax.texts if text.get_text() in _ON_PLOT_LABELS and text.get_visible()]
    if not labels or not stars:
        log.debug("Label pass skipped: %d labels, %d stars", len(labels), len(stars))
        return
    px = _pt(fig)
    ax_box = ax.get_window_extent(renderer)
    margin = 4.0 * px
    inner = Bbox.from_extents(
        ax_box.x0 + margin,
        ax_box.y0 + margin,
        ax_box.x1 - margin,
        ax_box.y1 - margin,
    )
    placed: list = []
    for text in labels:
        if text.get_text() not in _NUDGE_LABELS:
            placed.append(text)
            continue
        original = tuple(float(v) for v in text.get_position())
        anchor = ax.transData.transform(original)
        named = [star for star in stars if star[0] == text.get_text()]
        pool = named or stars
        _label, sx, sy, _radius = min(
            pool, key=lambda star: (star[1] - anchor[0]) ** 2 + (star[2] - anchor[1]) ** 2
        )
        best_pos = None
        best_score = -1e9
        for shift_x_pt, shift_y_pt in _LABEL_SHIFT_PT:
            display = (sx + shift_x_pt * px, sy + shift_y_pt * px)
            text.set_position(tuple(ax.transData.inverted().transform(display)))
            box = text.get_window_extent(renderer).padded(1.5 * px)
            if not inner.contains(box.x0, box.y0) or not inner.contains(box.x1, box.y1):
                continue
            gaps = [_clearance(box, mx, my, radius) for _name, mx, my, radius in stars]
            if path_samples:
                gaps.append(min(_clearance(box, px_, py_, radius) for px_, py_, radius in path_samples))
            for other in placed:
                if not other.get_visible():
                    continue
                other_box = other.get_window_extent(renderer)
                if _boxes_overlap(box, other_box, pad=3.0 * px):
                    gaps.append(-8.0 * px)
                else:
                    gaps.append(4.0 * px)
            for other in ax.texts:
                if other is text or not other.get_visible() or not _is_sigma_label(other.get_text()):
                    continue
                if _boxes_overlap(box, other.get_window_extent(renderer), pad=2.0 * px):
                    gaps.append(-8.0 * px)
            gap = min(gaps) if gaps else 0.0
            if gap < 6.0 * px:
                score = gap - 100.0 * px
            else:
                score = min(gap, 12.0 * px) - 0.08 * px * float(np.hypot(shift_x_pt, shift_y_pt))
            if score > best_score:
                best_score = score
                pos = text.get_position()
                best_pos = (float(pos[0]), float(pos[1]))
        text.set_position(best_pos if best_pos is not None else original)
        if best_pos is None:
            log.info("No clear spot for %s; left the original placement", text.get_text())
        if best_pos is not None:
            placed_disp = ax.transData.transform(best_pos)
            log.info(
                "Placed %s at offset (%.0f, %.0f) pt from its star",
                text.get_text(),
                (placed_disp[0] - sx) / px,
                (placed_disp[1] - sy) / px,
            )
        placed.append(text)


def _polish_landscape(fig) -> None:
    if id(fig) in _POLISHED:
        return
    axes = [
        ax
        for ax in fig.axes
        if any(_is_sigma_label(text.get_text()) for text in ax.texts)
        or any(text.get_text() in _NUDGE_LABELS for text in ax.texts)
    ]
    if not axes:
        return
    _POLISHED.add(id(fig))
    fig.canvas.draw()
    for ax in axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname("DejaVu Sans")
        for text in ax.texts:
            if text.get_zorder() < 20:
                text.set_zorder(20)
        _contain_sigma_labels(fig, ax)
        _separate_marker_labels(fig, ax)
    fig.canvas.draw()


def _install_savefig_hook() -> None:
    def savefig(self, *args, **kwargs):
        try:
            _polish_landscape(self)
        except Exception as exc:
            log.warning("Landscape label pass skipped: %s", exc)
        return _ORIG_SAVEFIG(self, *args, **kwargs)

    mpl.figure.Figure.savefig = savefig


def _parse_figsize(values: list[str]) -> tuple[float, float]:
    if len(values) != 2:
        msg = "--figsize expects width and height"
        raise SystemExit(msg)
    return float(values[0]), float(values[1])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--con-file", type=Path, required=True)
    parser.add_argument("--output-file", type=Path, required=True)
    parser.add_argument("--dat-pattern", required=True)
    parser.add_argument("--path-pattern", required=True)
    parser.add_argument("--cache-file", type=Path, required=True)
    parser.add_argument("--saddle", type=Path, required=True)
    parser.add_argument("--saddle-label", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--surface-type", required=True)
    parser.add_argument("--ira-kmax", type=float, default=14.0)
    parser.add_argument("--figsize", nargs=2, default=["5.37", "5.37"])
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--zoom-ratio", type=float, default=0.35)
    parser.add_argument("--fontsize-base", type=int, default=12)
    parser.add_argument("--facecolor", default="white")
    parser.add_argument("--rotation", default="15x,55y,0z")
    parser.add_argument("--strip-renderer", default="xyzrender")
    parser.add_argument("--xyzrender-config", default="paton")
    parser.add_argument("--perspective-tilt", type=float, default=8.0)
    parser.add_argument("--cmap-landscape", default="ruhi_diverging")
    parser.add_argument("--theme", default="ruhi")
    parser.add_argument("--plot-structures", default="crit_points")
    return parser


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    args = build_parser().parse_args()
    _install_font_hook()
    _install_savefig_hook()
    _use_open_zero_font()

    from rgpycrumbs.eon.plot_config import run_plot
    from rgpycrumbs.eon.plt_neb import plot_neb_from_settings

    output = args.output_file
    output.parent.mkdir(parents=True, exist_ok=True)
    run_plot(
        "neb",
        plot_neb_from_settings,
        con_file=args.con_file,
        output_file=output,
        input_dat_pattern=args.dat_pattern,
        input_path_pattern=args.path_pattern,
        cache_file=args.cache_file,
        additional_con=[(args.saddle, args.saddle_label)],
        title=args.title,
        surface_type=args.surface_type,
        ira_kmax=args.ira_kmax,
        figsize=_parse_figsize(args.figsize),
        dpi=args.dpi,
        zoom_ratio=args.zoom_ratio,
        fontsize_base=args.fontsize_base,
        facecolor=args.facecolor,
        rotation=args.rotation,
        strip_renderer=args.strip_renderer,
        xyzrender_config=args.xyzrender_config,
        perspective_tilt=args.perspective_tilt,
        cmap_landscape=args.cmap_landscape,
        theme=args.theme,
        plot_structures=args.plot_structures,
        plot_type="landscape",
        project_path=True,
        show_pts=True,
        landscape_path="all",
        show_legend=True,
        force_recompute=False,
    )
    if not output.is_file() or output.stat().st_size < 1000:
        msg = f"Landscape figure was not written: {output}"
        raise SystemExit(msg)
    log.info("Wrote %s", output)


if __name__ == "__main__":
    main()
