"""
Peak-shape (lineshape) functions shared by the static XAS and XES models.

Every spectroscopic peak here is the result of convolving a "natural"
lineshape with an instrument response:

- A core-excited state has a finite lifetime, which by the
  energy-time uncertainty relation gives it a Lorentzian ("cauchy") natural
  linewidth.
- The beamline/spectrometer resolution and any other statistical broadening
  (e.g. Doppler-like spread) is well approximated as Gaussian.
- The convolution of a Lorentzian and a Gaussian is a Voigt profile; when the
  Lorentzian width is negligible, a Gaussian alone is an adequate
  approximation (used for the pre-edge/near-edge features in
  ``xas_model``); when both widths matter (as for the XES emission lines),
  the full Voigt profile is used (``xes_model``).

See ``docs/physics/01_xas_xanes.md`` and ``docs/physics/02_xes_spin_states.md``
for the physical motivation.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.special as special

ArrayLike = npt.ArrayLike


def gaussian(x: ArrayLike, x0: float, sigma: float) -> np.ndarray:
    """
    Unit-height Gaussian centered at ``x0`` with standard deviation
    ``sigma``. A degenerate ``sigma = 0`` (which some stored fit results
    contain, for a Gaussian component the optimizer shrank to zero width)
    evaluates to a zero-width spike, i.e. 0 everywhere.
    """
    x = np.asarray(x, dtype=float)
    if sigma == 0:
        return np.zeros_like(x)
    return np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def cauchy(x: ArrayLike, sigma: float, x0: float, amplitude: float) -> np.ndarray:
    """Lorentzian (Cauchy) profile with half-width-at-half-maximum ``sigma``."""
    x = np.asarray(x, dtype=float)
    return (amplitude / (np.pi * sigma)) * (1 / (1 + ((x - x0) / sigma) ** 2))


def voigt(x: ArrayLike, amplitude: float, x0: float, gamma: float, sigma: float) -> np.ndarray:
    """
    Voigt profile: the convolution of a Gaussian (width ``sigma``) and a
    Lorentzian (HWHM ``gamma``), evaluated via the Faddeeva function.
    """
    x = np.asarray(x, dtype=float)
    z = (x - x0 + 1j * gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(special.wofz(z)) / (sigma * np.sqrt(2 * np.pi))


def bigaussian(x: ArrayLike, *p: float) -> np.ndarray:
    """Sum of two Gaussians plus a constant offset: p = (s1, x1, a1, s2, x2, a2, off)."""
    s1, x1, a1, s2, x2, a2, off = p
    return a1 * gaussian(x, x1, s1) + a2 * gaussian(x, x2, s2) + off


def pseudovoigt(x: ArrayLike, *p: float) -> np.ndarray:
    """Sum of a Gaussian and a Cauchy profile plus offset: p = (*gaussian(3), *cauchy(3), off)."""
    return gaussian(x, *p[1:3]) * p[0] + cauchy(x, *p[3:6]) + p[-1]


def erf_step(x: ArrayLike, x0: float, height: float, sigma: float) -> np.ndarray:
    """
    Gaussian-broadened step (edge jump) of the given ``height`` centered at
    ``x0``.

    This is the *exact* result of convolving an ideal Heaviside step with a
    Gaussian resolution function of standard deviation ``sigma`` -- see
    ``docs/physics/01_xas_xanes.md``. Used for the absorption-edge jump in
    ``xas_model.xas_model``.
    """
    x = np.asarray(x, dtype=float)
    return height * special.erf((x - x0) / (sigma * np.sqrt(2)))


def cdf(x: ArrayLike, sigma: float, x0: float, amplitude: float) -> np.ndarray:
    """Cumulative distribution function of a Gaussian, scaled by ``amplitude``."""
    x = np.asarray(x, dtype=float)
    return 0.5 * amplitude * special.erf((x - x0) / (sigma * np.sqrt(2)))


def ccdf(x: ArrayLike, sigma: float, x0: float, amplitude: float, offset: float) -> np.ndarray:
    """Complementary cumulative distribution function of a Gaussian."""
    x = np.asarray(x, dtype=float)
    return 0.5 * amplitude * special.erfc((x - x0) / (sigma * np.sqrt(2))) + offset


def linear(x: ArrayLike, slope: float, intercept: float) -> np.ndarray:
    """Straight line, ``slope * x + intercept``."""
    x = np.asarray(x, dtype=float)
    return slope * x + intercept


def biexp(x: ArrayLike, tau: float, offset: float, amplitude: float) -> np.ndarray:
    """Single exponential decay with time constant ``tau`` and an offset."""
    x = np.asarray(x, dtype=float)
    return amplitude * np.exp(-x / tau + offset)
