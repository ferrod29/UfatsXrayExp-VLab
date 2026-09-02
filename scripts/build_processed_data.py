"""
Turn a VLab session's raw exports into the ``.npy`` arrays the notebooks load.

Each session lives in its own dated folder::

    data/raw/<session>/...      exports straight out of the VLab
    data/processed/<session>/   the arrays built from them

Run from the repository root::

    python scripts/build_processed_data.py              # newest session
    python scripts/build_processed_data.py 20260901     # a specific one
    python scripts/build_processed_data.py --all        # every session

The loaders in ``vlab_xray.vlab_utils`` read the newest built session by
default; ``V.use_session("20260226")`` or ``$VLAB_SESSION`` picks another.

Raw formats produced by the VLab
--------------------------------
Every export is a text file with ``#``-comment headers, then either

* a 420x420 comma-separated pixel matrix (imagers: BIU2, BIU3, X-Ray Eye), or
* three whitespace-separated columns ``pixel  x  y`` (the Gotthard strip
  detectors at SpA1 and the von Hamos spectrometer, the TAD timing trace and
  the LPD jet scan), where ``x`` is photon energy in eV for the spectrometers
  and a position/time axis for TAD and LPD, or
* the IPM's own ``i=<n> x=<pos>`` rows followed by the four diode readings.

Session layouts
---------------
Two have been seen so far and are detected automatically:

``by_condition`` (20260226)
    Folders named for what was being done -- ``full_beam/``,
    ``focused_beam/``, ``pink_X-ray/``, ``laser_off/``, ``laser_on/<delay>/``.

``by_detector`` (20260901)
    Folders named for the detector, split by beam condition --
    ``BIU1/``, ``SPA1/{pink,mono}``, ``IPM/{pink,mono}``, ``TAD/{pink,mono}``,
    ``TR-XES/<delay>/``, plus ``LPD/`` and ``VonHamos/``.

Which folder holds the pink beam is **verified from the data**, never trusted
from the folder name -- see :func:`classify_beam`. In the 20260226 session the
two labels are swapped, and building it prints a warning saying so.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np

# Scale applied to the raw delay-folder names to get femtoseconds.
# 1.0 if the VLab scanned in fs, 1000.0 if it scanned in ps.
DELAY_UNIT_FS = 1.0

#: An acquisition whose integrated signal exceeds this multiple of the set's
#: median was taken with the attenuators out; it is saved separately rather
#: than averaged in, where it would swamp every other shot.
OUTLIER_FACTOR = 10.0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_ROOT = os.path.join(ROOT, "data", "raw")
PROCESSED_ROOT = os.path.join(ROOT, "data", "processed")


# ----------------------------------------------------------------------
# Raw readers
# ----------------------------------------------------------------------
def read_image(path):
    """420x420 comma-separated pixel matrix from an imager export."""
    return np.loadtxt(path, comments="#", delimiter=",")


def read_spectrum(path):
    """``(x, y)`` columns from a Gotthard/TAD/LPD export (column 0 is the pixel index)."""
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
    """The four diode readings some imager exports carry in their own header."""
    rows = []
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if "Readings from IPM" in line:
                rows.append([float(v) for v in line.split(":", 1)[1].split()])
                break
    return np.array(rows)


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------
class Session:
    """One dated raw folder and the processed folder it builds into."""

    def __init__(self, name):
        self.name = name
        self.raw = os.path.join(RAW_ROOT, name)
        self.out = os.path.join(PROCESSED_ROOT, name)

    def files(self, *parts):
        return sorted(glob.glob(os.path.join(self.raw, *parts)))

    def save(self, name, array):
        os.makedirs(self.out, exist_ok=True)
        np.save(os.path.join(self.out, name + ".npy"), array)
        print(f"  {name + '.npy':34s} {str(array.shape):16s} {array.dtype}")

    def save_image(self, name, array):
        """Imager arrays are stored float32: counts span 1e-4 to 5e9, which
        float32 holds to a relative error of 6e-8, and the frames are large
        enough that the halved file size matters."""
        self.save(name, np.asarray(array, dtype=np.float32))


def spectrum_pair(x, y):
    """Two columns ``(x, y)``, the layout ``vlab_utils.load_1d`` returns."""
    return np.column_stack([x, y])


def stack_of(paths):
    return np.array([read_image(p) for p in paths], dtype=np.float32)


def mean_spectrum(paths):
    pairs = [read_spectrum(p) for p in paths]
    return spectrum_pair(pairs[0][0], np.mean([y for _, y in pairs], axis=0))


def split_outliers(paths):
    """Partition acquisitions into (typical, unattenuated) by integrated signal."""
    if len(paths) < 3:
        return list(paths), []
    areas = np.array([np.trapezoid(*read_spectrum(p)[::-1]) for p in paths])
    bright = areas > OUTLIER_FACTOR * np.median(areas)
    return ([p for p, b in zip(paths, bright) if not b],
            [p for p, b in zip(paths, bright) if b])


def classify_beam(set_a, set_b):
    """Decide which of two SpA1 sets is the pink beam, from the data alone.

    A monochromator can neither broaden its input nor add flux, so the pink
    (SASE) beam is the one that is spectrally wider, varies shot to shot, and
    carries more photons. Returns ``(pink_paths, mono_paths, report)``.
    """
    def describe(paths):
        pairs = [read_spectrum(p) for p in paths]
        x = pairs[0][0]
        ys = np.array([y for _, y in pairs])
        support = float(np.ptp(x[ys.mean(0) > 0])) if (ys.mean(0) > 0).any() else 0.0
        corr = (float(np.corrcoef(ys[0], ys[1])[0, 1]) if len(ys) > 1 else float("nan"))
        area = float(np.median([np.trapezoid(y, x) for y in ys]))
        return support, corr, area

    a, b = describe(set_a), describe(set_b)
    # Two of three signatures agreeing is enough; they rarely disagree.
    votes = (a[0] > b[0]) + (a[2] > b[2]) + ((a[1] < b[1]) if a[1] == a[1] and b[1] == b[1] else 0)
    a_is_pink = votes >= 2
    pink, mono = (set_a, set_b) if a_is_pink else (set_b, set_a)
    p, m = (a, b) if a_is_pink else (b, a)
    report = (f"pink: support {p[0]:.1f} eV, shot corr {p[1]:.3f}, flux {p[2]:.2e}  |  "
              f"mono: support {m[0]:.1f} eV, shot corr {m[1]:.3f}, flux {m[2]:.2e}")
    return pink, mono, report


def emit_spa1(session, pink_paths, mono_paths):
    pink, pink_bright = split_outliers(pink_paths)
    mono, mono_bright = split_outliers(mono_paths)
    for i, path in enumerate(pink[:3]):
        session.save(f"spa1_pink_shot{i}", spectrum_pair(*read_spectrum(path)))
    session.save("spa1_pink_avg", mean_spectrum(pink))
    session.save("spa1_mono", mean_spectrum(mono))
    if pink_bright:
        session.save("spa1_pink_unattenuated", spectrum_pair(*read_spectrum(pink_bright[0])))
    if mono_bright:
        session.save("spa1_mono_unattenuated", spectrum_pair(*read_spectrum(mono_bright[0])))
    kept = len(pink_bright) + len(mono_bright)
    print(f"    ({len(pink)} pink + {len(mono)} mono shots averaged"
          + (f"; {kept} unattenuated shot(s) kept aside)" if kept else ")"))


def emit_delay_series(session, delay_dirs, pattern="*.txt"):
    """Average the repeats at each delay onto one common energy grid.

    A session may switch spectrometer range part-way through, leaving the
    delays on different grids; everything is interpolated onto the overlap.
    """
    delays, grids, means, sems, counts = [], [], [], [], []
    for name, directory in delay_dirs:
        paths = sorted(glob.glob(os.path.join(directory, pattern)))
        if not paths:
            continue
        pairs = [read_spectrum(p) for p in paths]
        x = pairs[0][0]
        ys = np.array([y for _, y in pairs])
        delays.append(float(name) * DELAY_UNIT_FS)
        grids.append(x)
        means.append(ys.mean(axis=0))
        sems.append(ys.std(axis=0, ddof=1) / np.sqrt(len(ys)) if len(ys) > 1
                    else np.full(x.size, np.nan))
        counts.append(len(ys))

    order = np.argsort(delays)
    delays = np.array(delays)[order]
    grids = [grids[i] for i in order]
    means = [means[i] for i in order]
    sems = [sems[i] for i in order]
    counts = np.array(counts)[order]

    lo = max(g.min() for g in grids)
    hi = min(g.max() for g in grids)
    n = max(g.size for g in grids)
    energy = np.linspace(lo, hi, n)
    regridded = not all(np.array_equal(g, grids[0]) for g in grids)
    matrix = np.array([np.interp(energy, g, m) for g, m in zip(grids, means)])
    sem = np.array([np.interp(energy, g, s) for g, s in zip(grids, sems)])

    session.save("transients_energy", energy)
    session.save("transients_delays", delays)
    session.save("transients_matrix", matrix)
    session.save("transients_sem", sem)
    session.save("transients_n_repeats", counts)
    print(f"    ({len(delays)} delays, {counts.min()}-{counts.max()} repeats each, averaged"
          + (f"; grids differed, interpolated onto {lo:.1f}-{hi:.1f} eV)" if regridded else ")"))


# ----------------------------------------------------------------------
# Layout: 20260226 -- folders named for the measurement condition
# ----------------------------------------------------------------------
def build_by_condition(s: Session):
    print("beam images (Ex. 2, 3, 9):")
    full = s.files("full_beam", "*.txt")
    s.save_image("biu2_stack", stack_of(full))
    s.save_image("biu2_pink", stack_of(full).mean(axis=0))
    # The single BIU2 frame filed with the narrow, dim spectra carries ~5000x
    # less flux than the full beam, which is what the four-bounce mono does.
    mono_frames = s.files("pink_X-ray", "BeamImagingUnit2*.txt")
    if mono_frames:
        s.save_image("biu2_mono", read_image(mono_frames[0]))
    s.save_image("biu3_focus_tight_stack", stack_of(s.files("focused_beam", "*.txt")))
    s.save_image("biu3_focus_wide_stack", stack_of(s.files("focused_beam", "New folder", "*.txt")))
    s.save_image("xray_eye", read_image(s.files("laser_off", "X-RayEye*.txt")[0]))
    s.save_image("xray_eye_unfocused_stack", stack_of(s.files("x-ray-eye", "*.txt")))

    print("\nincident spectra (Ex. 4, 5, 6):")
    pink, mono, report = classify_beam(s.files("pink_X-ray", "Gotthard SpA1*.txt"),
                                       s.files("monochromated_x-ray", "Gotthard SpA1*.txt"))
    print(f"    {report}")
    print("    folder labels agree with the data"
          if os.path.basename(os.path.dirname(pink[0])) == "pink_X-ray"
          else "    NOTE: labels are SWAPPED here -- 'monochromated_x-ray/'"
               " holds the pink beam. Naming by content.")
    emit_spa1(s, pink, mono)

    print("\nintensity/position monitor (Ex. 7):")
    headers = read_ipm_from_image_headers(full)
    standalone = s.files("IPM*.txt") or s.files("*", "IPM*.txt")
    if standalone:
        positions, diodes = read_ipm(standalone[0])
        s.save("ipm", diodes)
        s.save("ipm_scan_position", positions)
        if headers.size:
            s.save("ipm_biu2_headers", headers)
    elif headers.size:
        # No dedicated IPM export in this session; the BIU2 frames carry the
        # same four diode readings in their own header, one row per frame.
        s.save("ipm", headers)
        print("    (no standalone IPM export -- using the BIU2 frame headers)")
    else:
        print("    (no IPM data in this session)")

    print("\njet scan and ground state (Ex. 22, 23):")
    s.save("jet_scan", spectrum_pair(*read_spectrum(s.files("laser_off", "LPD*.txt")[0])))
    off = sorted(s.files("laser_off", "GotthardVonHamos*.txt"),
                 key=lambda p: read_spectrum(p)[1].sum())
    s.save("kb_ground", spectrum_pair(*read_spectrum(off[-1])))
    # The other laser-off acquisition integrates to ~0: the null measurement
    # showing the difference channel is flat with the pump blocked.
    if len(off) > 1:
        s.save("kb_laser_off_null", spectrum_pair(*read_spectrum(off[0])))

    print("\ntransient delay series (Ex. 29, 32-35):")
    root = os.path.join(s.raw, "laser_on")
    emit_delay_series(s, [(d, os.path.join(root, d)) for d in os.listdir(root)])


# ----------------------------------------------------------------------
# Layout: 20260901 -- folders named for the detector
# ----------------------------------------------------------------------
def build_by_detector(s: Session):
    print("beam images (Ex. 2):")
    # The folder is named BIU1 but the exports identify the detector as BIU2.
    frames = s.files("BIU1", "*.txt") or s.files("BIU2", "*.txt")
    if frames:
        s.save_image("biu2_stack", stack_of(frames))
        s.save_image("biu2_pink", stack_of(frames).mean(axis=0))
        print("    (no separate monochromatic BIU2 frame in this session, so Ex. 3's\n"
              "     image comparison stays synthetic; its flux argument uses SpA1 below)")

    print("\nincident spectra (Ex. 4, 5, 6):")
    pink, mono, report = classify_beam(s.files("SPA1", "pink", "*.txt"),
                                       s.files("SPA1", "mono", "*.txt"))
    print(f"    {report}")
    print("    folder labels agree with the data"
          if os.path.basename(os.path.dirname(pink[0])) == "pink"
          else "    NOTE: folder labels are swapped; naming by content.")
    emit_spa1(s, pink, mono)

    print("\nintensity/position monitor (Ex. 7):")
    for label in ("pink", "mono"):
        paths = s.files("IPM", label, "*.txt")
        if not paths:
            continue
        positions, diodes = read_ipm(paths[0])
        # Repeat exports of the same scan; stack them into one pulse train.
        allrows = np.vstack([read_ipm(p)[1] for p in paths])
        s.save("ipm" if label == "pink" else "ipm_mono", allrows)
        if label == "pink":
            s.save("ipm_scan_position", positions)

    print("\narrival-time monitor (temporal overlap):")
    for label in ("pink", "mono"):
        paths = s.files("TAD", label, "*.txt")
        if paths:
            s.save(f"tad_{label}", mean_spectrum(paths))

    print("\njet scan and ground state (Ex. 22, 23):")
    lpd = s.files("LPD", "*.txt")
    if lpd:
        s.save("jet_scan", spectrum_pair(*read_spectrum(lpd[0])))
    else:
        print("    (LPD/ is empty -- no jet scan in this session)")
    von_hamos = s.files("VonHamos", "*.txt")
    if von_hamos:
        order = sorted(von_hamos, key=lambda p: read_spectrum(p)[1].sum())
        s.save("kb_ground", spectrum_pair(*read_spectrum(order[-1])))
    else:
        print("    (VonHamos/ is empty -- no ground-state Kbeta spectrum in this session)")

    print("\ntransient delay series (Ex. 29, 32-35):")
    root = os.path.join(s.raw, "TR-XES")
    if os.path.isdir(root):
        emit_delay_series(s, [(d, os.path.join(root, d)) for d in os.listdir(root)])


# ----------------------------------------------------------------------
# Layout: beam profiles only (20263108)
# ----------------------------------------------------------------------
def build_beam_profiles(s: Session):
    print("beam images:")
    biu2 = s.files("beam_profile", "*.txt")
    if biu2:
        s.save_image("biu2_stack", stack_of(biu2))
        s.save_image("biu2_pink", stack_of(biu2).mean(axis=0))
    biu3 = s.files("beam_profile", "BIU3", "*.txt")
    if biu3:
        s.save_image("biu3_focus_tight_stack", stack_of(biu3))


# ----------------------------------------------------------------------
LAYOUTS = [
    ("by_condition", lambda s: os.path.isdir(os.path.join(s.raw, "laser_on")), build_by_condition),
    ("by_detector", lambda s: os.path.isdir(os.path.join(s.raw, "SPA1")), build_by_detector),
    ("beam_profiles", lambda s: os.path.isdir(os.path.join(s.raw, "beam_profile")), build_beam_profiles),
]


def detect_layout(session: Session):
    for name, matches, builder in LAYOUTS:
        if matches(session):
            return name, builder
    return None, None


def available_sessions():
    if not os.path.isdir(RAW_ROOT):
        return []
    return sorted(d for d in os.listdir(RAW_ROOT)
                  if os.path.isdir(os.path.join(RAW_ROOT, d)))


def build(name):
    session = Session(name)
    layout, builder = detect_layout(session)
    print("=" * 70)
    if builder is None:
        print(f"{name}: unrecognised layout, skipped "
              f"(contents: {', '.join(sorted(os.listdir(session.raw))[:6])})")
        return
    print(f"session {name}   layout: {layout}")
    print(f"  raw {session.raw}\n  out {session.out}\n")
    builder(session)
    n = len(glob.glob(os.path.join(session.out, "*.npy")))
    print(f"\n-> {n} arrays in data/processed/{name}/")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("session", nargs="?", help="dated folder under data/raw/ (default: newest)")
    parser.add_argument("--all", action="store_true", help="build every session found")
    args = parser.parse_args(argv)

    sessions = available_sessions()
    if not sessions:
        sys.exit(f"no session folders under {RAW_ROOT}")
    if args.all:
        targets = sessions
    elif args.session:
        if args.session not in sessions:
            sys.exit(f"unknown session {args.session!r}; available: {', '.join(sessions)}")
        targets = [args.session]
    else:
        targets = [sessions[-1]]
        print(f"(no session given -- building the newest, {targets[0]}; "
              f"available: {', '.join(sessions)})\n")

    for name in targets:
        build(name)


if __name__ == "__main__":
    main()
