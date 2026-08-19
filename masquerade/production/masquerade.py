"""
masquerade.py – Cell-mask overlay builder for multiplex TIFF images.

Memory-optimised rewrite
  • Coordinates are computed once and reused everywhere.
  • Full image is never materialised when only individual pages are needed.
  • Morphological dilation replaces manual pixel expansion.
  • Intermediate arrays are deleted and gc.collect()'d eagerly.
  • Deprecated scipy.ndimage.interpolation calls replaced with ndimage.zoom.
"""

import gc
import re

import numpy as np
import tifffile
from scipy import ndimage
from tifffile import TiffFile
from xml.etree import ElementTree


# ── helpers ──────────────────────────────────────────────────────────

def _compute_bounds(spatial_metadata):
    """Return (x_min, x_max, y_min, y_max) from spatial metadata."""
    xs = spatial_metadata["x"].astype(int)
    ys = spatial_metadata["y"].astype(int)
    return int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max())


def _apply_crop(arr_2d, bounds, adjust):
    """Crop a 2-D array to the bounding box when *adjust* is True."""
    if not adjust:
        return arr_2d
    x_min, x_max, y_min, y_max = bounds
    return arr_2d[y_min:y_max, x_min:x_max]


def _disk_structure(radius):
    """Filled circular structuring element of the given pixel radius."""
    y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return x**2 + y**2 <= radius**2


MASK_RADIUS = 4  # pixel radius of each per-cell mask circle


# ── public API ───────────────────────────────────────────────────────

def PreProcessImage(image_source, spatial_metadata, adjust_coords=True):
    """Load the TIFF, optionally crop to coordinate bounds.

    Returns
    -------
    image : ndarray  (C, H, W)
    raw_img_size : float  – approximate size in GB (uint8)
    bounds : tuple  – (x_min, x_max, y_min, y_max)
    """
    image = tifffile.imread(str(image_source))
    bounds = _compute_bounds(spatial_metadata)

    if adjust_coords:
        x_min, x_max, y_min, y_max = bounds
        image = image[:, y_min:y_max, x_min:x_max]

    raw_img_size = (np.prod(image.shape) * np.dtype("uint8").itemsize) / 1e9
    return image, raw_img_size, bounds


def get_mask_channels(
    image, spatial_metadata, raw_img_size, bounds, adjust_coords=True
):
    """Build per-cluster binary masks, dilate, compress.

    Returns
    -------
    channels : dict[str, ndarray]
    compression_factor : float
    """
    x_min, x_max, y_min, y_max = bounds
    cluster_ids = spatial_metadata["cluster"].unique()

    # Sum across channels → 2-D, then release the full image immediately.
    # image.sum() on uint16 promotes to float64 (8× per pixel), so freeing
    # the original C×H×W array before proceeding is critical.
    summed_image = image.sum(axis=0)  # (H, W)
    del image
    gc.collect()

    channels = {}
    compression_factor = 1.0

    for idx, cid in enumerate(cluster_ids):
        cluster = spatial_metadata[spatial_metadata["cluster"] == cid]
        if cluster.shape[0] <= 1:
            continue

        # Cluster pixel coordinates (optionally shifted to crop space)
        cx = cluster["x"].values.astype(int)
        cy = cluster["y"].values.astype(int)
        if adjust_coords:
            cx = cx - x_min
            cy = cy - y_min

        # Build boolean mask and dilate by 1 px (replaces manual ±1 expansion)
        mask = np.zeros(summed_image.shape, dtype=bool)
        # Clip to valid range
        cy_safe = np.clip(cy, 0, mask.shape[0] - 1)
        cx_safe = np.clip(cx, 0, mask.shape[1] - 1)
        mask[cy_safe, cx_safe] = True
        mask = ndimage.binary_dilation(mask, structure=_disk_structure(MASK_RADIUS))

        # Apply mask to the summed image
        masked = np.where(mask, summed_image, 0).astype(np.float64)

        # Compute compression factor once (same mask footprint assumed)
        if idx == 0:
            single_mask_gb = np.zeros_like(masked, dtype="uint8").nbytes / 1e9
            total_gb = single_mask_gb * len(cluster_ids) + raw_img_size
            compression_factor = min(1.0, np.sqrt(4.0 / total_gb))

        # Compress
        masked = ndimage.zoom(masked, compression_factor, order=3)
        channels[f"{cid}_mask-expanded"] = masked

    # Free the full image early – caller should not rely on it after this
    del summed_image
    gc.collect()

    return channels, compression_factor


def compress_marker_channels(
    image_source,
    channels,
    compression_factor,
    spatial_metadata,
    bounds,
    relevant_markers=None,
    adjust_coords=True,
):
    """Read individual TIFF pages (lazy), crop, compress, and add to *channels*.

    Reads one page at a time so memory never exceeds ~1 page + output.
    """
    x_min, x_max, y_min, y_max = bounds

    # Build the marker whitelist set (if provided)
    marker_set = None
    if relevant_markers is not None:
        raw = list(relevant_markers["x"])
        variants = set(raw)
        for m in raw:
            variants.add(m.replace("_", ""))
            variants.add(m.replace("_", "-"))
        marker_set = variants

    with TiffFile(str(image_source)) as tif:
        for page in tif.series[0].pages:
            desc = page.description
            if not desc:
                continue
            el = ElementTree.fromstring(desc).find("Biomarker")
            if el is None or el.text is None:
                continue

            name = el.text.replace(" ", "-")

            if name in channels:
                raise ValueError(
                    f"Duplicate channel name '{name}' – check TIFF metadata."
                )

            # If whitelist supplied, skip markers not in it
            if marker_set is not None and name not in marker_set:
                continue

            # Lazy read of single page
            arr = page.asarray()
            arr = _apply_crop(arr, (x_min, x_max, y_min, y_max), adjust_coords)
            arr = ndimage.zoom(arr.astype(np.float64), compression_factor, order=1)
            channels[name] = arr

            del arr
            gc.collect()

    return channels


def writeMaskTiff(channels, outPath):
    """Write channels page-by-page to an ImageJ-compatible TIFF.

    Avoids np.stack which would duplicate all channels into a single
    contiguous array, doubling peak memory at write time.
    """
    labels = list(channels.keys())

    with tifffile.TiffWriter(str(outPath), imagej=True) as tw:
        for k in labels:
            page = channels[k].astype("uint8")
            tw.write(page, contiguous=True, metadata={"Labels": labels} if k == labels[0] else None)
            del page
        gc.collect()
