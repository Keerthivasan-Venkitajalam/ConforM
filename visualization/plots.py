"""Self-contained SVG plots (no matplotlib dependency required).

Produces the conformational state-space explorer plot and simple bar charts
directly as SVG strings so they can be embedded in the HTML report and the
Streamlit dashboard without extra rendering dependencies.
"""
from __future__ import annotations

PALETTE = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2"]


def _scale(values, lo_out, hi_out):
    lo, hi = min(values), max(values)
    span = (hi - lo) or 1.0
    return [lo_out + (v - lo) / span * (hi_out - lo_out) for v in values]


def state_space_svg(pca_coords: list[list[float]], labels: list[str],
                    selected: str | None = None, pocket_states: list[str] | None = None,
                    width: int = 640, height: int = 420) -> str:
    """2D PCA state-space map. Marks which states contain the selected pocket."""
    if not pca_coords:
        return "<svg xmlns='http://www.w3.org/2000/svg'></svg>"
    pocket_states = pocket_states or []
    xs = [c[0] for c in pca_coords]
    ys = [c[1] if len(c) > 1 else 0.0 for c in pca_coords]
    pad = 70
    px = _scale(xs, pad, width - pad) if len(set(xs)) > 1 else [width / 2] * len(xs)
    py = _scale(ys, height - pad, pad) if len(set(ys)) > 1 else [height / 2] * len(ys)

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' "
        f"width='100%' style='max-width:{width}px;font-family:system-ui,sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='none'/>",
        f"<line x1='{pad-20}' y1='{height-pad+20}' x2='{width-pad+20}' y2='{height-pad+20}' "
        "stroke='currentColor' stroke-opacity='0.3'/>",
        f"<line x1='{pad-20}' y1='{pad-20}' x2='{pad-20}' y2='{height-pad+20}' "
        "stroke='currentColor' stroke-opacity='0.3'/>",
        f"<text x='{width/2}' y='{height-12}' text-anchor='middle' font-size='13' "
        "fill='currentColor' opacity='0.75'>PC1 (conformational state space)</text>",
        f"<text x='16' y='{height/2}' text-anchor='middle' font-size='13' fill='currentColor' "
        f"opacity='0.75' transform='rotate(-90 16 {height/2})'>PC2</text>",
    ]
    for i, label in enumerate(labels):
        has_pocket = label in pocket_states
        is_selected = label == selected
        color = "#E45756" if is_selected else ("#F58518" if has_pocket else "#4C78A8")
        r = 13 if is_selected else 9
        parts.append(f"<circle cx='{px[i]:.1f}' cy='{py[i]:.1f}' r='{r}' fill='{color}' "
                     f"fill-opacity='0.85' stroke='white' stroke-width='2'/>")
        if is_selected:
            parts.append(f"<circle cx='{px[i]:.1f}' cy='{py[i]:.1f}' r='{r+7}' fill='none' "
                         "stroke='#E45756' stroke-width='2' stroke-dasharray='4 3'/>")
        parts.append(f"<text x='{px[i]:.1f}' y='{py[i]-r-7:.1f}' text-anchor='middle' "
                     f"font-size='12' font-weight='600' fill='currentColor'>{label}</text>")

    legend = [("#E45756", "selected state (cryptic pocket)"), ("#F58518", "pocket-containing"),
              ("#4C78A8", "other sampled state")]
    for i, (color, text) in enumerate(legend):
        y = 22 + i * 19
        parts.append(f"<circle cx='{width-215}' cy='{y-4}' r='6' fill='{color}'/>")
        parts.append(f"<text x='{width-202}' y='{y}' font-size='11.5' fill='currentColor' "
                     f"opacity='0.85'>{text}</text>")
    parts.append("</svg>")
    return "".join(parts)


def bar_chart_svg(labels: list[str], values: list[float], title: str = "",
                  unit: str = "", width: int = 640, bar_h: int = 26,
                  lower_is_better: bool = False) -> str:
    if not labels:
        return "<svg xmlns='http://www.w3.org/2000/svg'></svg>"
    height = 44 + len(labels) * (bar_h + 8)
    label_w = 260
    max_abs = max(abs(v) for v in values) or 1.0
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' "
        f"width='100%' style='max-width:{width}px;font-family:system-ui,sans-serif'>"
    ]
    if title:
        parts.append(f"<text x='0' y='16' font-size='13' font-weight='600' "
                     f"fill='currentColor'>{title}</text>")
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 34 + i * (bar_h + 8)
        w = abs(value) / max_abs * (width - label_w - 70)
        best = (value == min(values)) if lower_is_better else (value == max(values))
        color = "#54A24B" if best else "#4C78A8"
        short = label if len(label) <= 34 else label[:31] + "..."
        parts.append(f"<text x='0' y='{y+bar_h*0.68:.0f}' font-size='11.5' "
                     f"fill='currentColor' opacity='0.9'>{short}</text>")
        parts.append(f"<rect x='{label_w}' y='{y}' width='{w:.1f}' height='{bar_h}' "
                     f"fill='{color}' rx='3'/>")
        parts.append(f"<text x='{label_w + w + 8:.1f}' y='{y+bar_h*0.68:.0f}' font-size='11.5' "
                     f"fill='currentColor' opacity='0.85'>{value:.2f}{unit}</text>")
    parts.append("</svg>")
    return "".join(parts)
