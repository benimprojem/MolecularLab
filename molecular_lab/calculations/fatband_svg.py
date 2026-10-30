"""Dependency-free SVG renderer for prepared fat-band plot data."""
from __future__ import annotations

from html import escape
from math import isfinite
from typing import Sequence

from molecular_lab.calculations.fatband_plot import FatBandPlotData


_DEFAULT_PALETTE = (
    "#2563eb", "#dc2626", "#059669", "#9333ea", "#d97706",
    "#0891b2", "#be185d", "#4d7c0f", "#7c3aed", "#0f766e",
)


def render_fatband_svg(
    data: FatBandPlotData,
    *,
    width: int = 1280,
    height: int = 800,
    energy_range_ev: tuple[float, float] | None = None,
    title: str = "Fat-band structure",
    palette: Sequence[str] = _DEFAULT_PALETTE,
) -> str:
    """Render plot-ready fat-band data as a self-contained SVG string.

    Thin neutral curves show the underlying band energies once per band/spin.
    Colored circles show the selected atom/orbital contribution; their radii
    scale with display_weight. No external plotting or GUI dependency is used.
    """
    if isinstance(width, bool) or isinstance(height, bool) or width < 500 or height < 350:
        raise ValueError("SVG dimensions must be at least 500x350 pixels")
    if not palette or any(not isinstance(color, str) or not color.strip() for color in palette):
        raise ValueError("palette must contain at least one non-empty color")
    if not data.kpoint_distances_inverse_angstrom:
        raise ValueError("fat-band plot data has no k-point path")
    all_energies = [
        point.energy_ev for series in data.series for point in series.points
    ]
    if not all_energies:
        raise ValueError("fat-band plot data contains no band energies")
    if energy_range_ev is None:
        low, high = min(all_energies), max(all_energies)
        if high - low < 1e-9:
            low, high = low - 1.0, high + 1.0
        else:
            padding = (high - low) * 0.06
            low, high = low - padding, high + padding
    else:
        low, high = energy_range_ev
        if not isfinite(low) or not isfinite(high) or high <= low:
            raise ValueError("energy_range_ev must be a finite increasing pair")

    left, right, top, bottom = 88.0, float(width - 230), 76.0, float(height - 92)
    x_values = data.kpoint_distances_inverse_angstrom
    x_min, x_max = x_values[0], x_values[-1]
    if x_max - x_min < 1e-12:
        x_max = x_min + 1.0

    def sx(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * (right - left)

    def sy(value: float) -> float:
        return bottom - (value - low) / (high - low) * (bottom - top)

    group_keys = sorted({series.group_key for series in data.series})
    group_labels = {
        series.group_key: series.label for series in data.series
    }
    colors = {key: palette[index % len(palette)] for index, key in enumerate(group_keys)}
    max_weight = max(
        (point.display_weight for series in data.series for point in series.points),
        default=0.0,
    )

    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<title>{escape(title)}</title>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2:.2f}" y="34" text-anchor="middle" font-size="22" font-family="sans-serif" font-weight="600">{escape(title)}</text>',
    ]

    # Horizontal energy grid and labels.
    for tick_index in range(6):
        energy = low + (high - low) * tick_index / 5
        y = sy(energy)
        chunks.append(
            f'<line x1="{left:.2f}" y1="{y:.2f}" x2="{right:.2f}" y2="{y:.2f}" stroke="#e5e7eb" stroke-width="1"/>'
        )
        chunks.append(
            f'<text x="{left - 12:.2f}" y="{y + 4:.2f}" text-anchor="end" font-size="12" fill="#374151" font-family="sans-serif">{energy:.3f}</text>'
        )

    # High-symmetry boundaries and labels.
    for tick in data.symmetry_ticks:
        x = sx(tick.path_distance_inverse_angstrom)
        chunks.append(
            f'<line x1="{x:.2f}" y1="{top:.2f}" x2="{x:.2f}" y2="{bottom:.2f}" stroke="#9ca3af" stroke-width="1" stroke-dasharray="4 4"/>'
        )
        chunks.append(
            f'<text x="{x:.2f}" y="{bottom + 25:.2f}" text-anchor="middle" font-size="12" fill="#111827" font-family="sans-serif">{escape(tick.label)}</text>'
        )

    # Draw the band-energy curves once per band and spin channel.
    drawn_bands: set[tuple[int, str]] = set()
    for series in data.series:
        band_key = (series.band_index, series.spin_channel)
        if band_key in drawn_bands:
            continue
        drawn_bands.add(band_key)
        coordinates = " ".join(
            f"{sx(point.path_distance_inverse_angstrom):.2f},{sy(point.energy_ev):.2f}"
            for point in series.points
        )
        chunks.append(
            f'<polyline points="{coordinates}" fill="none" stroke="#6b7280" stroke-width="1.25" opacity="0.75"/>'
        )

    # Overlay atom/orbital weights as contribution-sized circles.
    for series in data.series:
        color = colors[series.group_key]
        for point in series.points:
            if point.display_weight <= 0 or max_weight <= 0:
                continue
            ratio = min(1.0, point.display_weight / max_weight)
            radius = 1.5 + 6.5 * ratio ** 0.5
            x = sx(point.path_distance_inverse_angstrom)
            y = sy(point.energy_ev)
            chunks.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{escape(color)}" fill-opacity="0.72" stroke="#ffffff" stroke-width="0.45">'
                f'<title>{escape(series.label)} · band {series.band_index} · {series.spin_channel} · raw weight {point.raw_weight:.6g}</title></circle>'
            )

    # Axes.
    chunks.extend([
        f'<line x1="{left:.2f}" y1="{top:.2f}" x2="{left:.2f}" y2="{bottom:.2f}" stroke="#111827" stroke-width="1.4"/>',
        f'<line x1="{left:.2f}" y1="{bottom:.2f}" x2="{right:.2f}" y2="{bottom:.2f}" stroke="#111827" stroke-width="1.4"/>',
        f'<text x="{(left + right) / 2:.2f}" y="{height - 22}" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#111827">{escape(data.x_axis_label)}</text>',
        f'<text x="22" y="{(top + bottom) / 2:.2f}" text-anchor="middle" transform="rotate(-90 22 {(top + bottom) / 2:.2f})" font-size="14" font-family="sans-serif" fill="#111827">{escape(data.y_axis_label)}</text>',
    ])

    # Group legend.
    legend_x = right + 24
    legend_y = top + 12
    chunks.append(
        f'<text x="{legend_x:.2f}" y="{legend_y:.2f}" font-size="13" font-weight="600" font-family="sans-serif" fill="#111827">Orbital contribution</text>'
    )
    legend_y += 24
    for key in group_keys:
        chunks.append(
            f'<circle cx="{legend_x + 5:.2f}" cy="{legend_y - 4:.2f}" r="5" fill="{escape(colors[key])}"/>'
        )
        chunks.append(
            f'<text x="{legend_x + 17:.2f}" y="{legend_y:.2f}" font-size="11" font-family="sans-serif" fill="#374151">{escape(group_labels[key])}</text>'
        )
        legend_y += 20
        if legend_y > bottom - 8:
            break

    chunks.append("</svg>")
    return "\n".join(chunks)
