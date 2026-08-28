"""
vlab_utils.py
=============
Shared utilities for the *Ultrafast X-Ray Experiments* virtual-lab (VLab)
exercises (Versuchsmappe Exp21, FXE @ European XFEL).

The module gathers everything the four exercise notebooks reuse:

* physical constants and unit conversions (X-ray eV<->AA, optical nm<->eV);
* thin wrappers around ``xraydb`` that reproduce the CXRO transmission /
  attenuation data used throughout Section 5 (offline, no web access needed);
* small fitting helpers (Gaussian profile / FWHM, image centroid,
  area-normalisation, Bragg angles, IRF-convolved rise kinetics);
* graceful data loaders (``load_1d``, ``load_image``, ``load_delay_series``)
  that read a VLab export if present and otherwise fall back to the
  synthetic generators in :class:`Synthetic`, so every notebook is runnable
  end-to-end before any real data has been recorded.

Conventions (kept consistent with the ultrafast-pipeline skill)
---------------------------------------------------------------
* X-ray energies in eV, wavelengths in AA; optical wavelengths in nm.
* XES difference signal  dI = I(laser on) - I(laser off), area-normalised.
* Time delays in fs internally; a Gaussian instrument response (IRF) with
  FWHM ``irf_fwhm`` gates every kinetic model.

Author: generated for F. (UHH) -- teaching material, not a publication.
"""
from __future__ import annotations

import os
import glob
from dataclasses import dataclass

import numpy as np

from . import DATA_DIR as _PKG_DATA_DIR

# ----------------------------------------------------------------------
# Physical constants (CODATA, SI unless noted)
# ----------------------------------------------------------------------
H_PLANCK = 6.626_070_15e-34      # J s
HBAR = 1.054_571_817e-34         # J s
C_LIGHT = 2.997_924_58e8         # m/s
E_CHARGE = 1.602_176_634e-19     # C  (also: 1 eV in J)
N_A = 6.022_140_76e23            # 1/mol
K_B = 1.380_649e-23              # J/K
HC_EV_NM = 1239.841984          # eV nm  (h c / e, optical convenience)
HC_EV_ANG = 12398.41984         # eV AA  (X-ray convenience)

# Iron K-shell binding energy quoted in the script (eV)
FE_K_EDGE_EV = 7112.0
# Nominal Fe emission-line positions used for demo spectra (eV)
FE_KB13_EV = 7058.0             # K-beta_1,3 main line
FE_KBP_EV = 7045.0             # K-beta' spin-sensitive satellite


# ----------------------------------------------------------------------
# Unit conversions
# ----------------------------------------------------------------------
def xray_ev_to_angstrom(energy_ev):
    """X-ray photon energy [eV] -> wavelength [AA]."""
    return HC_EV_ANG / np.asarray(energy_ev, dtype=float)


def xray_angstrom_to_ev(lam_ang):
    """X-ray wavelength [AA] -> photon energy [eV]."""
    return HC_EV_ANG / np.asarray(lam_ang, dtype=float)


def optical_nm_to_ev(lam_nm):
    """Optical wavelength [nm] -> photon energy [eV]."""
    return HC_EV_NM / np.asarray(lam_nm, dtype=float)


def photon_energy_joule(lam_nm):
    """Optical wavelength [nm] -> single-photon energy [J]."""
    return H_PLANCK * C_LIGHT / (np.asarray(lam_nm, dtype=float) * 1e-9)


# ----------------------------------------------------------------------
# X-ray optical data (xraydb == offline CXRO)
# ----------------------------------------------------------------------
import xraydb  # noqa: E402


def mu_linear(formula, energy_ev, density=None):
    """Linear attenuation coefficient mu [1/cm] of *formula* at *energy_ev* [eV].

    Uses tabulated densities from xraydb when *density* is None (works for
    the elements / compounds registered there, e.g. 'Fe', 'H2O', 'C').
    """
    en = np.atleast_1d(np.asarray(energy_ev, dtype=float))
    mu = xraydb.material_mu(formula, en, density=density)
    if np.isscalar(energy_ev) or np.ndim(energy_ev) == 0:
        return float(np.asarray(mu).reshape(-1)[0])
    return mu


def transmission(formula, energy_ev, thickness_cm, density=None):
    """Beer-Lambert transmission T = exp(-mu * d) through *thickness_cm* [cm]."""
    mu = mu_linear(formula, energy_ev, density=density)
    return np.exp(-mu * thickness_cm)


def edge_energy(element, edge="K"):
    """Absorption-edge energy [eV] for *element* / *edge* (e.g. 'Fe','K')."""
    return xraydb.xray_edge(element, edge).energy


def delta_beta(formula, energy_ev, density=None):
    """Refractive-index decrement (delta, beta) with n = 1 - delta + i beta.

    Needed for the CRL focal-length estimate (Be lenses, Sec. 6.4.7).
    Unlike :func:`mu_linear`, *density* [g/cm^3] must be given explicitly --
    ``xraydb.xray_delta_beta`` has no density table to fall back on.
    """
    if density is None:
        raise ValueError(
            f"delta_beta({formula!r}, ...) needs an explicit density [g/cm^3] "
            "(e.g. density=1.848 for Be)."
        )
    en = np.atleast_1d(np.asarray(energy_ev, dtype=float))
    d, b, _ = xraydb.xray_delta_beta(formula, density, en)
    if np.isscalar(energy_ev) or np.ndim(energy_ev) == 0:
        return float(d[0]), float(b[0])
    return d, b


# ----------------------------------------------------------------------
# Crystal database + Bragg geometry (Sec. 4.5 / Ex. 21)
# ----------------------------------------------------------------------
# Cubic lattice constants [AA]
LATTICE_A = {"Si": 5.431020, "Ge": 5.657900, "C": 3.567}  # C = diamond


def d_spacing(material, hkl):
    """Interplanar spacing d_hkl [AA] for a cubic crystal."""
    a = LATTICE_A[material]
    h, k, l = hkl
    return a / np.sqrt(h * h + k * k + l * l)


def bragg_angle(energy_ev, material, hkl, order=1):
    """Bragg angle theta [deg] for reflection *hkl* at *energy_ev* [eV].

    Returns np.nan if the reflection cannot diffract the given energy
    (n*lambda / 2d > 1).
    """
    lam = xray_ev_to_angstrom(energy_ev)
    d = d_spacing(material, hkl)
    s = order * lam / (2.0 * d)
    with np.errstate(invalid="ignore"):
        return np.degrees(np.arcsin(np.where(s <= 1.0, s, np.nan)))


# ----------------------------------------------------------------------
# Chemistry helpers (Ex. 10-15, 19)
# ----------------------------------------------------------------------
from molmass import Formula  # noqa: E402


def molar_mass(formula):
    """Molar mass [g/mol] of a chemical *formula* string."""
    return Formula(formula).mass


# Building blocks for the sample
M_BIPY = molar_mass("C10H8N2")                  # 2,2'-bipyridine
M_FE_BIPY3 = molar_mass("Fe(C10H8N2)3")         # [Fe(bipy)3]2+ complex ion
M_FE_BIPY3_CL2 = molar_mass("FeCl2(C10H8N2)3")  # neutral salt
M_H2O = molar_mass("H2O")


# ----------------------------------------------------------------------
# Profile / peak fitting
# ----------------------------------------------------------------------
SIG2FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))     # sigma -> FWHM factor (~2.3548)


def gaussian(x, amp, cen, sigma, offset=0.0):
    """Gaussian with a constant offset."""
    return offset + amp * np.exp(-0.5 * ((x - cen) / sigma) ** 2)


def fit_gaussian(x, y):
    """Least-squares Gaussian fit. Returns dict with amp, cen, sigma, offset, fwhm."""
    from scipy.optimize import curve_fit

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    off0 = np.min(y)
    amp0 = np.max(y) - off0
    cen0 = x[np.argmax(y)]
    sig0 = max((x.max() - x.min()) / 6.0, np.finfo(float).eps)
    p0 = [amp0, cen0, sig0, off0]
    popt, pcov = curve_fit(gaussian, x, y, p0=p0, maxfev=20000)
    amp, cen, sigma, offset = popt
    return {"amp": amp, "cen": cen, "sigma": abs(sigma),
            "offset": offset, "fwhm": abs(sigma) * SIG2FWHM,
            "perr": np.sqrt(np.diag(pcov))}


def fwhm_from_data(x, y):
    """Model-free FWHM by linear interpolation of the half-maximum crossings."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    half = (y.max() + y.min()) / 2.0
    above = y >= half
    idx = np.where(np.diff(above.astype(int)) != 0)[0]
    if len(idx) < 2:
        return np.nan
    crossings = []
    for i in (idx[0], idx[-1]):
        x0, x1, y0, y1 = x[i], x[i + 1], y[i], y[i + 1]
        crossings.append(x0 + (half - y0) * (x1 - x0) / (y1 - y0))
    return abs(crossings[-1] - crossings[0])


def centroid_2d(image):
    """Intensity-weighted centroid (x_col, y_row) of a 2-D image, in pixels."""
    img = np.asarray(image, float)
    img = img - img.min()
    ys, xs = np.indices(img.shape)
    tot = img.sum()
    return xs.reshape(-1).dot(img.reshape(-1)) / tot, \
        ys.reshape(-1).dot(img.reshape(-1)) / tot


def profiles(image):
    """Return (col_profile, row_profile): sums along rows and columns."""
    img = np.asarray(image, float)
    return img.sum(axis=0), img.sum(axis=1)


# ----------------------------------------------------------------------
# Spectral processing
# ----------------------------------------------------------------------
def area_normalize(x, y):
    """Normalise y so that trapz(y, x) == 1 (area normalisation, Ex. 25)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    area = np.trapezoid(y, x)
    return y / area


# ----------------------------------------------------------------------
# Kinetics: IRF-convolved rise (Ex. 28, 34, 35)
# ----------------------------------------------------------------------
def irf_gaussian(t, t0, irf_fwhm):
    """Unit-area Gaussian instrument-response function centred at t0."""
    sigma = irf_fwhm / SIG2FWHM
    return np.exp(-0.5 * ((t - t0) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def _fine_grid(t, irf_fwhm, npts=4000):
    """Uniform fine grid spanning t padded by +-6 sigma of the IRF."""
    pad = 6 * irf_fwhm / SIG2FWHM
    return np.linspace(np.min(t) - pad, np.max(t) + pad, npts)


def _convolve_with_irf(signal, fine, irf_fwhm):
    """Convolve a signal defined on the (uniform) *fine* grid with the IRF."""
    kernel = irf_gaussian(fine, fine[fine.size // 2], irf_fwhm)
    kernel = kernel / kernel.sum()
    return np.convolve(signal, kernel, mode="same")


def irf_step_exp(t, t0, tau, irf_fwhm, amp=1.0):
    """Saturating rise ``amp*(1 - exp(-(t-t0)/tau))`` convolved with a Gaussian IRF.

    Numeric convolution (consistent with :func:`sequential_populations`),
    valid for any IRF/tau ratio. Used for the HS-population rise fit (Ex. 34).
    """
    t = np.asarray(t, float)
    fine = _fine_grid(t, irf_fwhm)
    step = np.where(fine >= 0.0, 1.0 - np.exp(-np.clip(fine, 0, None) / tau), 0.0)
    conv = _convolve_with_irf(step, fine, irf_fwhm)
    return amp * np.interp(t - t0, fine, conv)


def sequential_populations(t, t0, taus, irf_fwhm):
    """Populations of an A->B->...->Z sequential scheme, IRF-convolved.

    *taus* are the lifetimes of the successive precursor states; the final
    state is a stable product. State 0 receives unit population at t=0; the
    linear kinetic chain is integrated on a fine grid and each population is
    convolved with the Gaussian IRF. Returns array (n_states, n_times).
    """
    t = np.asarray(t, float)
    taus = np.asarray(taus, float)
    n = taus.size + 1
    k = 1.0 / taus
    fine = _fine_grid(t, irf_fwhm)
    dt = fine[1] - fine[0]
    i0 = int(np.searchsorted(fine, 0.0))        # first index with fine >= 0
    P = np.zeros((n, fine.size))
    P[0, i0] = 1.0                              # inject unit population at t=0
    for j in range(i0 + 1, fine.size):
        prev = P[:, j - 1]
        dP = np.zeros(n)
        for i in range(n - 1):
            rate = k[i] * prev[i]
            dP[i] -= rate
            dP[i + 1] += rate
        P[:, j] = prev + dP * dt
    out = np.zeros((n, t.size))
    for i in range(n):
        conv = _convolve_with_irf(P[i], fine, irf_fwhm)
        out[i] = np.clip(np.interp(t - t0, fine, conv), 0.0, None)
    return out


# ----------------------------------------------------------------------
# Data loading with synthetic fallback
# ----------------------------------------------------------------------
#: Repository ``data/`` directory, shared with the rest of the package so the
#: loaders below find a real VLab export regardless of the notebook's current
#: working directory.
DATA_DIR = str(_PKG_DATA_DIR)


#: Where ``scripts/build_processed_data.py`` writes the ``.npy`` arrays it
#: builds from ``data/raw/``. Searched before ``data/`` itself, so a rebuilt
#: export wins over anything left lying in the top-level data directory.
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")


def _exists(path):
    return path is not None and os.path.exists(path)


def _candidates(name, extensions):
    """Every path to try for *name*, processed arrays first, then ``data/``."""
    for directory in (PROCESSED_DIR, DATA_DIR):
        for ext in extensions:
            yield os.path.join(directory, name + ext)


def _read_any(path):
    """Read a ``.npy`` array or a text table, whichever *path* is."""
    if path.endswith(".npy"):
        return np.load(path)
    return _loadtxt(path, delimiter="," if path.endswith(".csv") else None)


def _loadtxt(path, delimiter=None):
    """``np.loadtxt`` that tolerates the text headers VLab exports carry.

    Lines are skipped from the top until the rest of the file parses as
    numbers, so a one- or two-line column caption (as in the CXRO exports)
    does not need to be stripped by hand.
    """
    for skip in range(5):
        try:
            return np.loadtxt(path, delimiter=delimiter, skiprows=skip)
        except ValueError:
            continue
    # Report the original, most informative error.
    return np.loadtxt(path, delimiter=delimiter)


def load_1d(name, gen=None, **gen_kw):
    """Load a 2-column (x, y) text file ``data/<name>`` or synthesise it.

    Looks in ``data/processed/`` first, then ``data/``. Recognised
    extensions: .npy (an (n, 2) array) and .txt/.csv/.dat (whitespace or
    comma separated). If no file is found and *gen* (a callable) is given, it
    is called with **gen_kw and its (x, y) return value is used instead.

    Returns ``(x, y, path)``, where *path* is the file actually read or
    ``None`` when the synthetic generator was used -- so a notebook can
    always say which of the two it is showing.
    """
    for path in _candidates(name, ("", ".npy", ".txt", ".csv", ".dat")):
        if _exists(path):
            arr = _read_any(path)
            return arr[:, 0], arr[:, 1], path
    if gen is not None:
        x, y = gen(**gen_kw)
        return x, y, None
    raise FileNotFoundError(f"{name} not found in {PROCESSED_DIR} or {DATA_DIR} "
                            "and no generator given")


def load_image(name, gen=None, **gen_kw):
    """Load a 2-D image (or a stack of them) or synthesise it.

    Looks for ``data/processed/<name>`` first, then ``data/<name>``, trying
    the .npy/.txt/.csv extensions. Returns ``(image, path)``, with *path*
    ``None`` when synthetic.
    """
    for path in _candidates(name, ("", ".npy", ".txt", ".csv")):
        if _exists(path):
            return _read_any(path), path
    if gen is not None:
        return gen(**gen_kw), None
    raise FileNotFoundError(f"{name} not found in {PROCESSED_DIR} or {DATA_DIR} "
                            "and no generator given")


def load_table(name, gen=None, **gen_kw):
    """Load a plain n-column table (e.g. the four IPM diode signals).

    Same search as :func:`load_1d` but returns the whole array rather than
    splitting off two columns: ``(array, path)``, *path* ``None`` when the
    generator was used.
    """
    for path in _candidates(name, ("", ".npy", ".txt", ".csv", ".dat")):
        if _exists(path):
            return _read_any(path), path
    if gen is not None:
        return gen(**gen_kw), None
    raise FileNotFoundError(f"{name} not found in {PROCESSED_DIR} or {DATA_DIR} "
                            "and no generator given")


def load_delay_series(subdir, gen=None, **gen_kw):
    """Load a delay series from ``data/<subdir>/`` or synthesise it.

    Two layouts are accepted, in this order:

    1. the three arrays ``data/processed/<subdir>_energy.npy``,
       ``_delays.npy`` and ``_matrix.npy``, as written by
       ``scripts/build_processed_data.py``;
    2. a directory ``data/<subdir>/`` holding one 2-column file per delay
       named ``<delay_fs>.txt`` (e.g. ``-150.txt``, ``100.txt``).

    Returns ``(energy, delays_fs, matrix, path)`` with matrix shape
    (n_delays, n_energy), *path* being the file/directory actually read or
    ``None`` when synthetic.
    """
    bundle = [os.path.join(PROCESSED_DIR, f"{subdir}_{part}.npy")
              for part in ("energy", "delays", "matrix")]
    if all(_exists(p) for p in bundle):
        energy, delays, matrix = (np.load(p) for p in bundle)
        order = np.argsort(delays)
        return energy, delays[order], matrix[order], bundle[2]

    d = os.path.join(DATA_DIR, subdir)
    files = sorted(glob.glob(os.path.join(d, "*.txt")),
                   key=lambda p: float(os.path.splitext(os.path.basename(p))[0]))
    if files:
        delays, rows, energy = [], [], None
        for p in files:
            arr = _loadtxt(p)
            energy = arr[:, 0]
            delays.append(float(os.path.splitext(os.path.basename(p))[0]))
            rows.append(arr[:, 1])
        return energy, np.array(delays), np.array(rows), d
    if gen is not None:
        e, t, m = gen(**gen_kw)
        return e, t, m, None
    raise FileNotFoundError(f"no delay files in {d} and no generator given")


# ----------------------------------------------------------------------
# Synthetic data generators (physically-motivated stand-ins for VLab output)
# ----------------------------------------------------------------------
class Synthetic:
    """Deterministic, physically-plausible stand-ins for VLab exports.

    Every generator seeds its own RNG so notebooks are reproducible. These
    are *teaching* stand-ins: peak positions and relative intensities follow
    Figs. 4.3/4.4 and 3.6/3.7, but absolute numbers are illustrative.
    """

    # --- Kb emission spectra -------------------------------------------------
    @staticmethod
    def kb_spectrum(energy=None, spin_fraction=0.0, noise=0.0, seed=0):
        """Fe K-beta emission spectrum vs spin.

        *spin_fraction* in [0,1] mixes a low-spin (0) and high-spin (1)
        line shape: the K-beta' satellite at ~7045 eV grows and the main
        K-beta_1,3 line shifts slightly to lower energy with spin.
        """
        if energy is None:
            energy = np.linspace(7020, 7080, 601)
        f = float(np.clip(spin_fraction, 0, 1))
        main_cen = FE_KB13_EV - 0.8 * f          # slight red shift with spin
        main_w = 2.6 + 0.9 * f                    # broadening with spin
        y = gaussian(energy, 1.0, main_cen, main_w)
        # spin-sensitive K-beta' satellite grows in with multiplicity
        y += gaussian(energy, 0.16 * f, FE_KBP_EV, 3.0)
        y += gaussian(energy, 0.05 * f, 7035.0, 4.0)  # low-energy tail
        if noise:
            rng = np.random.default_rng(seed)
            y = y + rng.normal(0, noise * y.max(), size=y.shape)
        return energy, y

    @staticmethod
    def reference_series(energy=None, noise=0.0, seed=1):
        """Kb reference spectra for multiplicities 1..5 (singlet..quintet)."""
        if energy is None:
            energy = np.linspace(7020, 7080, 601)
        mults = [1, 2, 3, 4, 5]
        specs = {}
        for m in mults:
            frac = (m - 1) / 4.0
            _, y = Synthetic.kb_spectrum(energy, spin_fraction=frac,
                                         noise=noise, seed=seed + m)
            specs[m] = area_normalize(energy, y)
        return energy, specs

    # --- SASE / mono incident spectra ---------------------------------------
    @staticmethod
    def sase_single_shot(energy=None, e0=7250.0, bw_rel=1e-3, n_spikes=110,
                         seed=0):
        """One spiky SASE single-shot spectrum around *e0* [eV]."""
        if energy is None:
            energy = np.linspace(e0 - 25, e0 + 25, 1000)
        rng = np.random.default_rng(seed)
        sigma = e0 * bw_rel / SIG2FWHM
        env = np.exp(-0.5 * ((energy - e0) / sigma) ** 2)
        y = np.zeros_like(energy)
        cens = rng.normal(e0, sigma, n_spikes)
        for c in cens:
            amp = rng.exponential(1.0)
            y += amp * gaussian(energy, 1.0, c, sigma / 12.0)
        y *= env
        return energy, y / y.max()

    @staticmethod
    def sase_average(energy=None, e0=7250.0, bw_rel=1e-3, n_shots=400):
        """Average of many SASE shots -> smooth (near-Gaussian) envelope."""
        if energy is None:
            energy = np.linspace(e0 - 25, e0 + 25, 1000)
        acc = np.zeros_like(energy)
        for s in range(n_shots):
            _, y = Synthetic.sase_single_shot(energy, e0, bw_rel, seed=s)
            acc += y
        return energy, acc / acc.max()

    @staticmethod
    def mono_spectrum(energy=None, e0=7250.0, darwin_rel=1.3e-4, seed=0):
        """Monochromatised beam: narrow (Si(111) Darwin width) profile."""
        if energy is None:
            energy = np.linspace(e0 - 4, e0 + 4, 800)
        sigma = e0 * darwin_rel / SIG2FWHM
        rng = np.random.default_rng(seed)
        y = gaussian(energy, 1.0, e0, sigma)
        y = y + rng.normal(0, 0.01, size=y.shape)
        return energy, np.clip(y, 0, None)

    # --- Beam images / jitter ------------------------------------------------
    @staticmethod
    def beam_image(shape=(120, 160), fwhm_px=(28, 40), center=None,
                   jitter=(0.0, 0.0), noise=0.02, seed=0):
        """A 2-D Gaussian X-ray beam image (BIU2) with optional pointing jitter."""
        rng = np.random.default_rng(seed)
        ny, nx = shape
        if center is None:
            center = (nx / 2, ny / 2)
        cx = center[0] + rng.normal(0, jitter[0])
        cy = center[1] + rng.normal(0, jitter[1])
        ys, xs = np.indices(shape)
        sx = fwhm_px[0] / SIG2FWHM
        sy = fwhm_px[1] / SIG2FWHM
        img = np.exp(-0.5 * (((xs - cx) / sx) ** 2 + ((ys - cy) / sy) ** 2))
        img = img + rng.normal(0, noise, size=shape)
        return np.clip(img, 0, None)

    @staticmethod
    def biu2_stack(n=20, jitter_px=(6.0, 0.6), seed=0):
        """Stack of *n* single-pulse BIU2 images with realistic jitter.

        Horizontal jitter ~one beam diameter, vertical ~10x smaller
        (Sec. 6.2.2). Returns array (n, ny, nx)."""
        return np.array([
            Synthetic.beam_image(jitter=jitter_px, seed=seed + i)
            for i in range(n)
        ])

    @staticmethod
    def ipm_signals(n=200, seed=0):
        """Four-diode IPM signals (I0..I3) for *n* pulses.

        Diodes arranged as a quadrant; horizontal position from (I1-I3),
        vertical from (I0-I2); sum ~ pulse energy. Encodes horizontal jitter
        larger than vertical, plus intensity fluctuations.
        """
        rng = np.random.default_rng(seed)
        # true beam position (mm) and intensity per pulse
        x = rng.normal(0.0, 0.05, n)     # horizontal jitter, mm (larger)
        y = rng.normal(0.0, 0.015, n)    # vertical jitter, mm (smaller)
        E = rng.normal(1.0, 0.08, n)     # relative pulse energy
        # geometry: foil-to-diode-plane L=40 mm, half-spacing a=20 mm
        L, a = 40.0, 20.0
        # simple linear response: fraction on each diode
        gx = 0.5 * x / a
        gy = 0.5 * y / a
        I0 = E * (0.25 + gy)   # top
        I2 = E * (0.25 - gy)   # bottom
        I1 = E * (0.25 + gx)   # right
        I3 = E * (0.25 - gx)   # left
        eps = rng.normal(0, 0.003, (n, 4))
        return np.clip(np.stack([I0, I1, I2, I3], axis=1) + eps, 0, None), \
            {"L": L, "a": a, "x_true": x, "y_true": y, "E_true": E}

    # --- Time-resolved XES ---------------------------------------------------
    @staticmethod
    def transient_map(energy=None, delays=None, irf_fwhm=80.0,
                      tau_mlct=25.0, tau_intermediate=130.0, seed=0, noise=0.001):
        """Time-resolved Fe K-beta *difference* map dI(E, t).

        Physical stand-in following Fig. 3.6/3.7: after excitation the system
        proceeds 1,3MLCT -> (intermediate 3T) -> 5T2(HS) via a short
        sequential cascade convolved with a Gaussian IRF. The difference
        spectrum of the *intermediate* differs slightly from that of the HS
        product, so the data support a two-component (SVD) description and a
        detectable ultrafast intermediate (Ex. 35).

        Returns (energy, delays_fs, matrix[n_delays, n_energy]).
        """
        if energy is None:
            energy = np.linspace(7020, 7080, 301)
        if delays is None:
            delays = np.concatenate([
                np.linspace(-200, 0, 9),
                np.array([25, 50, 75, 100, 150, 200, 300, 450, 600, 800])
            ])
        # difference spectra of the two transient species vs the LS ground state
        e, ls = Synthetic.kb_spectrum(energy, spin_fraction=0.0)
        _, hs = Synthetic.kb_spectrum(energy, spin_fraction=1.0)       # 5T2
        # intermediate (~3T, S=1): K-beta' only partially grown AND sitting
        # ~1.5 eV higher than the HS satellite -> a genuinely distinct shape,
        # not just a scaled-down HS spectrum (so SVD resolves two components).
        inter = gaussian(energy, 1.0, FE_KB13_EV - 0.4, 2.9)
        inter += gaussian(energy, 0.09, FE_KBP_EV + 1.5, 3.0)
        ls = area_normalize(energy, ls)
        inter = area_normalize(energy, inter)
        hs = area_normalize(energy, hs)
        d_inter = inter - ls
        d_hs = hs - ls
        # populations: state0 (invisible MLCT) -> state1 (intermediate) -> state2 (HS)
        # taus are the lifetimes of the successive *precursor* states, so the
        # first entry drains the MLCT into the intermediate and the second
        # drains the intermediate into the HS product.
        pops = sequential_populations(delays, 0.0,
                                      [tau_mlct, tau_intermediate], irf_fwhm)
        # states: 0 = MLCT (no Kb' yet, ~LS-like, invisible in difference)
        #         1 = intermediate 3T  -> d_inter
        #         2 = HS 5T2 product   -> d_hs
        M = (np.outer(pops[1], d_inter) + np.outer(pops[2], d_hs))
        rng = np.random.default_rng(seed)
        M = M + rng.normal(0, noise, size=M.shape)
        return energy, np.asarray(delays, float), M


# ----------------------------------------------------------------------
# Plot styling
# ----------------------------------------------------------------------
def use_paper_style():
    """Apply a compact, publication-leaning matplotlib style."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.dpi": 120,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.5,
    })
