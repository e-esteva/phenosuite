"""
viewer.py – Masquerade napari viewer with interactive magicgui control panel.

Launch with no arguments — a napari window opens immediately with a
docked panel for loading files and controlling layer visibility:

    python viewer.py

Requirements:
    pip install "napari[all]" tifffile numpy
    pip install zarr dask pandas   # optional features
"""

from __future__ import annotations

import sys
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import tifffile
from tifffile import TiffFile

try:
    import napari
    from magicgui import widgets as mw
    from napari.qt.threading import thread_worker
except ImportError as exc:
    sys.exit(
        f"Missing dependency: {exc}\n"
        "Install with:  pip install 'napari[all]'"
    )


# ── Palette + colormaps ───────────────────────────────────────────────

CLUSTER_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe",
    "#469990", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9",
]

_MARKER_CMAPS: list[tuple[str, str]] = [
    ("dapi",    "blue"),
    ("hoechst", "blue"),
    ("cd3",     "green"),
    ("cd8",     "cyan"),
    ("cd4",     "magenta"),
    ("cd20",    "yellow"),
    ("cd68",    "red"),
    ("foxp3",   "bop orange"),
    ("ki67",    "bop purple"),
    ("pd1",     "bop blue"),
    ("pdl1",    "bop blue"),
]


def _pick_colormap(name: str) -> str:
    lname = name.lower()
    for prefix, cmap in _MARKER_CMAPS:
        if lname.startswith(prefix):
            return cmap
    return "gray"


def _hex_to_rgba(h: str) -> tuple[float, float, float, float]:
    h = h.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) / 255 for i in (0, 2, 4))
    return (r, g, b, 1.0)


# ── Data loading (runs in background thread) ──────────────────────────

def _parse_biomarker(desc: str | None) -> str | None:
    """Extract Biomarker name from PerkinElmer/Akoya per-page XML."""
    if not desc:
        return None
    try:
        el = ElementTree.fromstring(desc).find("Biomarker")
        if el is not None and el.text:
            return el.text.strip().replace(" ", "-")
    except ElementTree.ParseError:
        pass
    return None


def _load_qptiff(path: str, lazy: bool) -> list[tuple]:
    """Load QPTIFF channels. Returns list of (array, name, colormap, store_ref)."""
    results = []

    if lazy:
        try:
            import zarr
            import dask.array as da

            with TiffFile(path) as tif:
                names = [
                    _parse_biomarker(getattr(p, "description", "")) or f"ch_{i:03d}"
                    for i, p in enumerate(tif.series[0].pages)
                ]
            store = tifffile.ZarrTiffStore(path, series=0, level=0)
            z = zarr.open(store, mode="r")
            arr = da.from_zarr(z["0"] if isinstance(z, zarr.hierarchy.Group) else z)

            if arr.ndim == 3:
                for i in range(arr.shape[0]):
                    n = names[i] if i < len(names) else f"ch_{i:03d}"
                    results.append((arr[i], n, _pick_colormap(n), store))
            else:
                n = names[0] if names else "ch_000"
                results.append((arr, n, _pick_colormap(n), store))
            return results

        except Exception as exc:
            print(f"[viewer] Lazy load failed ({exc}) — falling back to eager.")

    with TiffFile(path) as tif:
        for i, page in enumerate(tif.series[0].pages):
            desc = getattr(page, "description", "") or ""
            name = _parse_biomarker(desc) or f"ch_{i:03d}"
            results.append((page.asarray(), name, _pick_colormap(name), None))
    return results


def _load_masks(path: str) -> tuple[dict, dict]:
    """Split a Masquerade output TIFF into cluster masks and marker channels."""
    cluster_masks: dict[str, np.ndarray] = {}
    marker_channels: dict[str, np.ndarray] = {}

    with TiffFile(path) as tif:
        labels: list[str] = []
        if tif.imagej_metadata and "Labels" in tif.imagej_metadata:
            raw = tif.imagej_metadata["Labels"]
            labels = list(raw) if not isinstance(raw, str) else [raw]

        for i, page in enumerate(tif.series[0].pages):
            name = labels[i] if i < len(labels) else f"layer_{i:03d}"
            arr = page.asarray()
            if name.endswith("_mask-expanded"):
                cluster_masks[name] = arr
            else:
                marker_channels[name] = arr

    return cluster_masks, marker_channels


def _build_label_array(cluster_masks: dict) -> tuple[np.ndarray, list[str]]:
    names = list(cluster_masks.keys())
    ref = next(iter(cluster_masks.values()))
    label_img = np.zeros(ref.shape, dtype=np.int32)
    for idx, (_, mask) in enumerate(cluster_masks.items(), start=1):
        label_img[mask > 0] = idx
    return label_img, names


@thread_worker
def _load_worker(img_path: str, mask_path: str | None, csv_path: str | None, lazy: bool) -> dict:
    """Background worker — loads all data and returns a result dict."""
    result: dict = {
        "image_channels": [],   # [(arr, name, cmap, store_ref), ...]
        "label_img":      None,
        "cluster_names":  [],
        "color_dict":     {},
        "marker_channels": {},  # compressed markers from mask TIFF
        "point_data":     [],   # [(coords_yx, cluster_id, color), ...]
        "error":          None,
    }
    try:
        result["image_channels"] = _load_qptiff(img_path, lazy)

        if mask_path:
            cluster_masks, marker_channels = _load_masks(mask_path)
            if cluster_masks:
                label_img, cluster_names = _build_label_array(cluster_masks)
                color_dict = {0: (0.0, 0.0, 0.0, 0.0)}
                color_dict.update({
                    i + 1: _hex_to_rgba(CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)])
                    for i in range(len(cluster_names))
                })
                result["label_img"]     = label_img
                result["cluster_names"] = cluster_names
                result["color_dict"]    = color_dict
            result["marker_channels"] = marker_channels

        if csv_path:
            import pandas as pd
            meta = pd.read_csv(csv_path)
            if {"x", "y", "cluster"}.issubset(meta.columns):
                for i, cid in enumerate(sorted(meta["cluster"].unique())):
                    sub    = meta[meta["cluster"] == cid]
                    coords = sub[["y", "x"]].values   # napari: (row, col) = (y, x)
                    color  = CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
                    result["point_data"].append((coords, cid, color))

    except Exception as exc:
        result["error"] = str(exc)

    return result


# ── Control panel ─────────────────────────────────────────────────────

def create_masquerade_panel(viewer: napari.Viewer) -> mw.Container:
    """Build the docked magicgui control panel and wire up all callbacks."""

    # ── tracked layer names (so we can clear on reload) ──
    _state: dict = {
        "image_layers": [],
        "mask_layer":   None,
        "point_layers": [],
        "store_refs":   [],   # zarr stores: must stay alive while dask arrays are used
    }

    # ── file pickers ──────────────────────────────────────────────────
    img_edit  = mw.FileEdit(label="Source Image",     filter="*.qptiff *.tiff *.tif", mode="r")
    mask_edit = mw.FileEdit(label="Mask TIFF (opt.)", filter="*.tiff *.tif",          mode="r")
    csv_edit  = mw.FileEdit(label="Spatial CSV (opt.)", filter="*.csv",               mode="r")
    lazy_cb   = mw.CheckBox(text="Lazy / dask  (large files, requires zarr + dask)")
    load_btn  = mw.PushButton(text="Load")
    status_lbl = mw.Label(value="Select an image and click Load.")

    # ── image channel controls (shown after load) ─────────────────────
    ch_show_btn = mw.PushButton(text="Show All")
    ch_hide_btn = mw.PushButton(text="Hide All")
    ch_row = mw.Container(widgets=[ch_show_btn, ch_hide_btn], layout="horizontal")
    ch_section = mw.Container(widgets=[
        mw.Label(value="Image Channels"),
        ch_row,
    ])
    ch_section.visible = False

    # ── cluster / mask controls (shown after mask or CSV load) ────────
    opacity_sl = mw.FloatSlider(label="Mask opacity", min=0.0, max=1.0, value=0.55, step=0.05)
    pt_show_btn = mw.PushButton(text="Show Points")
    pt_hide_btn = mw.PushButton(text="Hide Points")
    pt_row = mw.Container(widgets=[pt_show_btn, pt_hide_btn], layout="horizontal")
    cluster_section = mw.Container(widgets=[
        mw.Label(value="Cluster Controls"),
        opacity_sl,
        pt_row,
    ])
    cluster_section.visible = False

    # ── assemble panel ────────────────────────────────────────────────
    panel = mw.Container(widgets=[
        mw.Label(value="MASQUERADE VIEWER"),
        img_edit,
        mask_edit,
        csv_edit,
        lazy_cb,
        load_btn,
        status_lbl,
        ch_section,
        cluster_section,
    ])

    # ── callbacks ─────────────────────────────────────────────────────

    def _remove_tracked_layers() -> None:
        names = (
            _state["image_layers"]
            + ([_state["mask_layer"]] if _state["mask_layer"] else [])
            + _state["point_layers"]
        )
        for name in names:
            try:
                viewer.layers.remove(name)
            except (KeyError, ValueError):
                pass
        _state["image_layers"].clear()
        _state["point_layers"].clear()
        _state["mask_layer"] = None
        _state["store_refs"].clear()

    def _on_result(result: dict) -> None:
        """Called on the main Qt thread when the worker finishes."""
        if result["error"]:
            status_lbl.value = f"Error: {result['error']}"
            load_btn.enabled = True
            return

        _remove_tracked_layers()

        # Image channels
        for i, (arr, name, cmap, store_ref) in enumerate(result["image_channels"]):
            viewer.add_image(
                arr,
                name=name,
                colormap=cmap,
                blending="additive",
                visible=(i == 0),   # show only first channel (usually DAPI)
            )
            _state["image_layers"].append(name)
            if store_ref is not None:
                _state["store_refs"].append(store_ref)

        # Cluster labels layer
        if result["label_img"] is not None:
            lname = "Cluster Masks"
            viewer.add_labels(
                result["label_img"],
                name=lname,
                color=result["color_dict"],
                opacity=opacity_sl.value,
            )
            _state["mask_layer"] = lname

        # Compressed marker channels from mask TIFF (hidden by default)
        for ch_name, ch_arr in result["marker_channels"].items():
            full_name = f"[mask] {ch_name}"
            viewer.add_image(
                ch_arr,
                name=full_name,
                colormap=_pick_colormap(ch_name),
                blending="additive",
                visible=False,
            )
            _state["image_layers"].append(full_name)

        # Cell scatter points (hidden by default, toggle per cluster)
        for coords, cid, color in result["point_data"]:
            pname = f"Cluster {cid}"
            viewer.add_points(
                coords,
                name=pname,
                size=6,
                face_color=color,
                edge_color="transparent",
                visible=False,
            )
            _state["point_layers"].append(pname)

        # Show/hide control sections
        ch_section.visible      = bool(_state["image_layers"])
        cluster_section.visible = bool(_state["mask_layer"] or _state["point_layers"])

        n_ch = len(result["image_channels"])
        n_cl = len(result["cluster_names"])
        n_pt = len(result["point_data"])
        parts = [f"{n_ch} channel(s)"]
        if n_cl:
            parts.append(f"{n_cl} cluster mask(s)")
        if n_pt:
            parts.append(f"{n_pt} point layer(s)")
        status_lbl.value = "Loaded: " + ", ".join(parts)
        load_btn.enabled = True

    def _on_error(exc: Exception) -> None:
        status_lbl.value = f"Error: {exc}"
        load_btn.enabled = True

    def _on_load() -> None:
        img_path = Path(str(img_edit.value))
        if not img_path.is_file():
            status_lbl.value = "Please select a valid source image."
            return

        mask_path = str(mask_edit.value)
        mask_path = mask_path if Path(mask_path).is_file() else None

        csv_path = str(csv_edit.value)
        csv_path = csv_path if Path(csv_path).is_file() else None

        status_lbl.value = "Loading..."
        load_btn.enabled = False

        worker = _load_worker(str(img_path), mask_path, csv_path, bool(lazy_cb.value))
        worker.returned.connect(_on_result)
        worker.errored.connect(_on_error)
        worker.start()

    def _set_image_visibility(visible: bool) -> None:
        for name in _state["image_layers"]:
            if name in viewer.layers:
                viewer.layers[name].visible = visible

    def _on_opacity(val: float) -> None:
        if _state["mask_layer"] and _state["mask_layer"] in viewer.layers:
            viewer.layers[_state["mask_layer"]].opacity = val

    def _set_point_visibility(visible: bool) -> None:
        for name in _state["point_layers"]:
            if name in viewer.layers:
                viewer.layers[name].visible = visible

    load_btn.clicked.connect(_on_load)
    ch_show_btn.clicked.connect(lambda: _set_image_visibility(True))
    ch_hide_btn.clicked.connect(lambda: _set_image_visibility(False))
    opacity_sl.changed.connect(_on_opacity)
    pt_show_btn.clicked.connect(lambda: _set_point_visibility(True))
    pt_hide_btn.clicked.connect(lambda: _set_point_visibility(False))

    return panel


# ── Entry point ───────────────────────────────────────────────────────

def main() -> None:
    viewer = napari.Viewer(title="Masquerade Viewer")
    panel  = create_masquerade_panel(viewer)
    viewer.window.add_dock_widget(panel, name="Masquerade", area="left")
    napari.run()


if __name__ == "__main__":
    main()
