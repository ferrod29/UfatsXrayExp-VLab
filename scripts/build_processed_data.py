"""
Turn the raw VLab exports in ``data/raw/`` into the ``.npy`` arrays the
notebooks load from ``data/processed/``.

Run from the repository root::

    python scripts/build_processed_data.py

The script is idempotent: it only reads ``data/raw/`` and rewrites
``data/processed/``, so it can be re-run after adding more exports.

Raw formats produced by the VLab
--------------------------------
Every export is a text file with ``#``-comment headers, then either

* a 420x420 comma-separated pixel matrix (imagers: BIU2, BIU3, X-Ray Eye), or
* three whitespace-separated columns ``pixel  x  y`` (the Gotthard strip
  detectors at SpA1 and the von Hamos spectrometer, and the LPD jet scan),
  where ``x`` is photon energy in eV for the spectrometers and jet position
  for the LPD, or
* the IPM's own ``i=<n> x=<pos>`` rows followed by the four diode readings.

Two facts about this data set that the folder names do not tell you
-------------------------------------------------------------------
1. ``data/raw/pink_X-ray/`` and ``data/raw/monochromated_x-ray/`` are
   **swapped**. Measured from the files themselves:

   =========================  ====================  =========================
   quantity                   ``pink_X-ray/``       ``monochromated_x-ray/``
   =========================  ====================  =========================
   bandwidth dE/E             3.2e-4 (narrow)       2.4-4.5e-3 (broad)
   shot-to-shot correlation   1.000 (identical)     0.52 (spiky, varies)
   integrated flux            1.0e8                 4.0e9  (39x brighter)
   =========================  ====================  =========================

   A monochromator can neither broaden its input nor add flux, so the broad,
   spiky, brighter set is the pink/SASE beam and the narrow, reproducible,
   dimmer one is the monochromatised beam. This script names the outputs by
   what the data are, not by which folder they sit in.

2. The delay folder names under ``data/raw/laser_on/`` (-360 ... +360, step
   15) carry no unit. They are treated as **femtoseconds**, matching the
   repository convention and the shape of the rise (noise before about
   -180 fs, half maximum near +80 fs, plateau by about +300 fs). Set
   ``DELAY_UNIT_FS`` below if the VLab scan was actually in picoseconds.
"""

from __future__ import annotations

import glob
import os
import re

import numpy as np

# Scale applied to the raw delay-folder names to get femtoseconds.
# 1.0 if the VLab scanned in fs, 1000.0 if it scanned in ps.
DELAY_UNIT_FS = 1.0

#: An acquisition whose integrated signal exceeds this multiple of the set's
#: median was taken with the attenuators out; it is saved separately rather
#: than averaged in, where it would swamp every other shot.
OUTLIER_FACTOR = 10.0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "data", "processed")


# ----------------------------------------------------------------------
# Raw readers
# ----------------------------------------------------------------------
def read_image(path):
    """420x420 comma-separated pixel matrix from an imager export."""
    return np.loadtxt(path, comments="#", delimiter=",")


def read_spectrum(path):
    """``(x, y)`` columns from a Gotthard/LPD export (column 0 is the pixel index)."""
    a = np.loadtxt(path, comments="#")
    return a[:, 1], a[:, 2]


def read_ipm(path):
    """``(scan_position, diodes)`` from the IPM export's ``i=<n> x=<pos>`` rows."""
    positions, diodes = [], []
    for line in open(path, encoding="utf-8"):
        if not line.startswith("i="):
            continue
        positions.append(float(re.search(r"x=(-?\d+\.?\d*)", line).group(1)))
        diodes.append([float(v) for v in line.split("\t")[1:] if v.strip()])
    return np.array(positions), np.array(diodes)


def read_ipm_from_image_headers(paths):
    """The four diode readings the BIU2 exports carry in their own header."""
    rows = []
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if "Readings from IPM" in line:
                rows.append([float(v) for v in line.split(":", 1)[1].split()])
                break
    return np.array(rows)


def files(*parts):
    return sorted(glob.glob(os.path.join(RAW, *parts)))


def save(name, array):
    path = os.path.join(OUT, name + ".npy")
    np.save(path, array)
    print(f"  {name + '.npy':34s} {str(array.shape):16s} {array.dtype}")
    return path


def split_outliers(paths):
    """Partition acquisitions into (typical, unattenuated) by integrated signal."""
    areas = np.array([np.trapezoid(*read_spectrum(p)[::-1]) for p in paths])
    bright = areas > OUTLIER_FACTOR * np.median(areas)
    return [p for p, b in zip(paths, bright) if not b], [p for p, b in zip(paths, bright) if b]


def stack_of(paths):
    return np.array([read_image(p) for p in paths], dtype=np.float32)


def save_image(name, array):
    """Imager arrays are stored as float32: detector counts span 1e-4 to 5e9,
    which float32 holds to a relative error of 6e-8, and the frames are large
    enough that the halved file size matters."""
    return save(name, np.asarray(array, dtype=np.float32))


def spectrum_pair(x, y):
    """Two columns ``(x, y)``, the layout ``vlab_utils.load_1d`` returns."""
    return np.column_stack([x, y])


def mean_spectrum(paths):
    xs = [read_spectrum(p) for p in paths]
    return spectrum_pair(xs[0][0], np.mean([y for _, y in xs], axis=0))


# ----------------------------------------------------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"reading  {RAW}\nwriting  {OUT}\n")

    # --- Imagers -------------------------------------------------------
    print("beam images (Ex. 2, 3, 9):")
    full = files("full_beam", "*.txt")
    save_image("biu2_stack", stack_of(full))                    # Ex. 2: jitter statistics
    save_image("biu2_pink", stack_of(full).mean(axis=0))        # Ex. 3: full (pink) beam
    # The single BIU2 frame filed under pink_X-ray/ belongs with the narrow,
    # dim spectra in that folder -- it carries ~5000x less flux than the full
    # beam above, which is what inserting the four-bounce mono does.
    save_image("biu2_mono", read_image(files("pink_X-ray", "BeamImagingUnit2*.txt")[0]))

    save_image("biu3_focus_tight_stack", stack_of(files("focused_beam", "*.txt")))
    save_image("biu3_focus_wide_stack", stack_of(files("focused_beam", "New folder", "*.txt")))

    # Ex. 9 wants the focused spot at the sample; the x-ray-eye/ set is the
    # same detector before focusing (~85 px FWHM vs ~2.9 px here).
    save_image("xray_eye", read_image(files("laser_off", "X-RayEye*.txt")[0]))
    save_image("xray_eye_unfocused_stack", stack_of(files("x-ray-eye", "*.txt")))

    # --- Incident spectra (SpA1) ---------------------------------------
    # Named by measured character, not by folder -- see the module docstring.
    print("\nincident spectra (Ex. 4, 5, 6)  [folder labels corrected]:")
    pink_paths, pink_bright = split_outliers(files("monochromated_x-ray", "Gotthard SpA1*.txt"))
    mono_paths, mono_bright = split_outliers(files("pink_X-ray", "Gotthard SpA1*.txt"))

    for i, path in enumerate(pink_paths[:3]):
        save(f"spa1_pink_shot{i}", spectrum_pair(*read_spectrum(path)))
    save("spa1_pink_avg", mean_spectrum(pink_paths))
    save("spa1_mono", mean_spectrum(mono_paths))
    if pink_bright:
        save("spa1_pink_unattenuated", spectrum_pair(*read_spectrum(pink_bright[0])))
    if mono_bright:
        save("spa1_mono_unattenuated", spectrum_pair(*read_spectrum(mono_bright[0])))
    print(f"    ({len(pink_paths)} pink + {len(mono_paths)} mono shots averaged; "
          f"{len(pink_bright) + len(mono_bright)} unattenuated shot(s) kept aside)")

    # --- IPM -----------------------------------------------------------
    print("\nintensity/position monitor (Ex. 7):")
    positions, diodes = read_ipm(files("IPM*.txt")[0])
    save("ipm", diodes)                     # (n_pulses, 4): I0 I1 I2 I3
    save("ipm_scan_position", positions)
    save("ipm_biu2_headers", read_ipm_from_image_headers(full))

    # --- Jet scan and ground-state emission ----------------------------
    print("\njet scan and ground state (Ex. 22, 23):")
    save("jet_scan", spectrum_pair(*read_spectrum(files("laser_off", "LPD*.txt")[0])))

    von_hamos_off = files("laser_off", "GotthardVonHamos*.txt")
    by_signal = sorted(von_hamos_off, key=lambda p: read_spectrum(p)[1].sum())
    save("kb_ground", spectrum_pair(*read_spectrum(by_signal[-1])))
    # The other laser-off acquisition integrates to ~0: the null measurement
    # that shows the difference channel is flat with the pump blocked.
    save("kb_laser_off_null", spectrum_pair(*read_spectrum(by_signal[0])))

    # --- Delay series --------------------------------------------------
    print("\ntransient delay series (Ex. 29, 32-35):")
    delay_dirs = sorted(os.listdir(os.path.join(RAW, "laser_on")), key=float)
    energy, rows, sems, counts = None, [], [], []
    for name in delay_dirs:
        paths = files("laser_on", name, "*.txt")
        spectra = [read_spectrum(p) for p in paths]
        energy = spectra[0][0]
        ys = np.array([y for _, y in spectra])
        rows.append(ys.mean(axis=0))
        # Repeat counts vary between delays (most have 5, one has 7), so the
        # per-repeat spectra cannot be stacked into one rectangular array;
        # keep the standard error and the count instead.
        sems.append(ys.std(axis=0, ddof=1) / np.sqrt(len(ys)))
        counts.append(len(ys))
    save("transients_energy", energy)
    save("transients_delays", np.array([float(d) for d in delay_dirs]) * DELAY_UNIT_FS)
    save("transients_matrix", np.array(rows))
    save("transients_sem", np.array(sems))
    save("transients_n_repeats", np.array(counts))
    print(f"    ({len(delay_dirs)} delays, {min(counts)}-{max(counts)} repeats each, averaged)")

    print(f"\ndone -- {len(glob.glob(os.path.join(OUT, '*.npy')))} arrays in data/processed/")


if __name__ == "__main__":
    main()
