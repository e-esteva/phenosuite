"""
viewer_utils.py – In-browser channel compositing for the Masquerade Shiny viewer.

Loaded alongside masquerade.py via reticulate. Maintains a module-level
channel cache (safe: reticulate gives each R session its own Python process).
"""

from __future__ import annotations

import base64
import colorsys
import io

import matplotlib
matplotlib.use("Agg")                  # headless – no display needed
import matplotlib.image as mpimg
import numpy as np
import tifffile
from tifffile import TiffFile


# ── Default colours per channel type ─────────────────────────────────

# 20 perceptually-distinct cluster colours
_CLUSTER_RGB: list[tuple[float, float, float]] = [
    (0.902, 0.098, 0.294),  # red
    (0.235, 0.706, 0.294),  # green
    (1.000, 0.882, 0.098),  # yellow
    (0.263, 0.388, 0.847),  # blue
    (0.961, 0.510, 0.192),  # orange
    (0.569, 0.118, 0.706),  # purple
    (0.259, 0.831, 0.957),  # cyan
    (0.941, 0.196, 0.902),  # magenta
    (0.749, 0.937, 0.271),  # lime
    (0.980, 0.745, 0.745),  # pink
    (0.502, 0.306, 0.165),  # brown
    (0.200, 0.800, 0.800),  # teal
    (0.902, 0.502, 0.000),  # amber
    (0.502, 0.000, 0.502),  # dark purple
    (0.000, 0.502, 0.502),  # dark teal
    (0.800, 0.000, 0.200),  # crimson
    (0.200, 0.600, 0.000),  # forest green
    (0.000, 0.200, 0.800),  # navy
    (0.600, 0.400, 0.000),  # sienna
    (0.400, 0.800, 0.400),  # sage
]

# Named colours for well-known markers (prefix-matched, case-insensitive)
_MARKER_RGB: list[tuple[str, tuple[float, float, float]]] = [
    ("dapi",    (0.20, 0.40, 1.00)),
    ("hoechst", (0.20, 0.40, 1.00)),
    ("cd3",     (0.20, 0.90, 0.20)),
    ("cd8",     (0.20, 0.90, 0.90)),
    ("cd4",     (0.90, 0.20, 0.90)),
    ("cd20",    (1.00, 1.00, 0.20)),
    ("cd68",    (1.00, 0.20, 0.20)),
    ("foxp3",   (1.00, 0.60, 0.20)),
    ("ki67",    (0.70, 0.20, 0.90)),
    ("pd1",     (0.20, 0.50, 0.90)),
    ("pdl1",    (0.20, 0.50, 0.90)),
    ("epcam",   (1.00, 0.50, 0.00)),
    ("panck",   (1.00, 0.50, 0.00)),
    ("ck",      (1.00, 0.50, 0.00)),
    ("cd45",    (0.00, 0.80, 0.40)),
    ("cd56",    (0.60, 0.00, 0.80)),
    ("cd31",    (0.80, 0.40, 0.00)),
    ("cd163",   (0.80, 0.20, 0.20)),
    ("cd11b",   (0.60, 0.80, 0.20)),
    ("sma",     (0.20, 0.80, 0.80)),
]

# Fallback: 36 evenly-spaced hues at high saturation & value, for any
# unrecognised marker (so every channel gets a distinct non-grey colour)
_AUTO_MARKER_RGB: list[tuple[float, float, float]] = [
    colorsys.hsv_to_rgb(i / 36, 0.85, 0.95) for i in range(36)
]


def _default_rgb(
    name: str,
    mask_index: int,
    auto_marker_index: int,
) -> tuple[float, float, float]:
    if name.endswith("_mask-expanded"):
        return _CLUSTER_RGB[mask_index % len(_CLUSTER_RGB)]
    lname = name.lower()
    for prefix, rgb in _MARKER_RGB:
        if lname.startswith(prefix):
            return rgb
    # Auto-assign a distinct hue for any panel marker not in the list above
    return _AUTO_MARKER_RGB[auto_marker_index % len(_AUTO_MARKER_RGB)]


# ── Module-level channel cache ────────────────────────────────────────

_channels: dict[str, np.ndarray] = {}    # name → float32 (H, W)
_channel_meta: list[dict] = []           # ordered metadata returned to R


def load_tiff_channels(path: str) -> list[dict]:
    """Read a Masquerade output TIFF into the cache.

    Returns a list of dicts — one per channel — with keys:
        name, is_mask, color_r, color_g, color_b
    R receives this as a list of named lists.
    """
    global _channels, _channel_meta
    _channels = {}
    _channel_meta = []

    with TiffFile(path) as tif:
        labels: list[str] = []
        if tif.imagej_metadata and "Labels" in tif.imagej_metadata:
            raw = tif.imagej_metadata["Labels"]
            labels = list(raw) if not isinstance(raw, str) else [raw]

        mask_idx = 0
        auto_marker_idx = 0
        for i, page in enumerate(tif.series[0].pages):
            name = labels[i] if i < len(labels) else f"layer_{i:03d}"
            arr  = page.asarray().astype(np.float32)
            _channels[name] = arr

            is_mask = name.endswith("_mask-expanded")
            r, g, b = _default_rgb(name, mask_idx, auto_marker_idx)
            if is_mask:
                mask_idx += 1
            else:
                auto_marker_idx += 1

            _channel_meta.append({
                "name":    name,
                "is_mask": bool(is_mask),
                "color_r": float(r),
                "color_g": float(g),
                "color_b": float(b),
            })

    return _channel_meta


def get_image_dims() -> list[int]:
    """Return [height, width] of the cached image."""
    if not _channels:
        return [0, 0]
    arr = next(iter(_channels.values()))
    return [int(arr.shape[0]), int(arr.shape[1])]


# ── Compositing ───────────────────────────────────────────────────────

def composite_and_encode(
    visible_names: list[str],
    color_r:       list[float],
    color_g:       list[float],
    color_b:       list[float],
    brightnesses:  list[float],
) -> str:
    """Additively composite visible channels and return a base64 PNG string.

    Parameters are parallel lists indexed by position in visible_names.
    Returns "" if there are no visible channels or no cached data.
    """
    if not _channels or not visible_names:
        return ""

    ref = next(iter(_channels.values()))
    H, W = ref.shape
    out = np.zeros((H, W, 3), dtype=np.float32)

    for idx, name in enumerate(visible_names):
        arr = _channels.get(name)
        if arr is None:
            continue
        vmax = float(arr.max())
        if vmax <= 0:
            continue

        layer = (arr / vmax) * float(brightnesses[idx])
        out[:, :, 0] += layer * float(color_r[idx])
        out[:, :, 1] += layer * float(color_g[idx])
        out[:, :, 2] += layer * float(color_b[idx])

    np.clip(out, 0.0, 1.0, out=out)
    rgba = np.dstack([out, np.ones((H, W), dtype=np.float32)])   # H×W×4, 0-1

    buf = io.BytesIO()
    mpimg.imsave(buf, rgba, format="png")
    buf.seek(0)
    return base64.b64encode(buf.getvalue()).decode("utf-8")
