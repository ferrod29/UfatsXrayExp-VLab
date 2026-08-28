"""
Time-resolved (pump-probe) kinetics: convolving an instrument response
function (IRF) with a two-step population cascade, and building transient
XAS/XES maps from it.

The experiment excites a population with a short laser pulse (modeled as a
Gaussian in time, the IRF, of width ``sigma`` centered at ``t0``) which then
relaxes through a cascade of states with lifetime ``tau``. Because the pulse
has finite duration, the population appearing in each state as a function of
pump-probe delay ``t`` is the convolution of the Gaussian pulse with the
state's step-and-decay/step-and-rise population kinetics -- this is the
classic "exponentially modified Gaussian" building block of ultrafast
spectroscopy. See ``docs/physics/03_ultrafast_pump_probe_kinetics.md``.

This module deduplicates the identical convolution machinery that was
previously copy-pasted between the XAS and XES transient notebooks.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.integrate import quad


def _excited_state_kernel(x, delay, t0, i0, sigma, tau):
    return (
        (i0 / (sigma * (2 * np.pi) ** 0.5))
        * np.exp(-(x**2) / (2 * sigma**2))
        * np.heaviside(delay - t0 - x, 0.5)
        * np.exp(-(delay - t0 - x) / tau)
    )


def _product_state_kernel(x, delay, t0, i0, sigma, tau):
    return (
        (i0 / (sigma * (2 * np.pi) ** 0.5))
        * np.exp(-(x**2) / (2 * sigma**2))
        * np.heaviside(delay - t0 - x, 0.5)
        * (1 - np.exp(-(delay - t0 - x) / tau))
    )


#: Half-width, in units of `sigma`, of the integration window used for the
#: convolution integrals below. The integration variable `x` in the kernels
#: above is the pulse-time offset from the Gaussian's own center (0, not
#: `t0`), so beyond a few `sigma` the Gaussian factor makes the integrand
#: negligible regardless of `t0`, `tau`, or which delay is being evaluated --
#: 10 sigma leaves a safety margin many orders of magnitude past that point.
_INTEGRATION_HALF_WIDTH_IN_SIGMA = 10.0


def excited_state_population(
    delays: npt.ArrayLike, t0: float, i0: float, sigma: float, tau: float
) -> np.ndarray:
    """
    Population still present in the (first, directly pumped) excited state
    at each pump-probe ``delays``: rises with the Gaussian pump pulse
    (center ``t0``, width ``sigma``) and decays with lifetime ``tau``.
    """
    delays = np.asarray(delays, dtype=float)
    bound = _INTEGRATION_HALF_WIDTH_IN_SIGMA * sigma
    return np.array(
        [quad(_excited_state_kernel, -bound, bound, args=(t, t0, i0, sigma, tau))[0] for t in delays]
    )


def product_state_population(
    delays: npt.ArrayLike, t0: float, i0: float, sigma: float, tau: float
) -> np.ndarray:
    """
    Population that has transferred out of the excited state into the next
    state down the cascade, at each pump-probe ``delays`` -- the complement
    of :func:`excited_state_population`: it rises towards ``i0`` with the
    same time constant ``tau`` that drains the excited state.
    """
    delays = np.asarray(delays, dtype=float)
    bound = _INTEGRATION_HALF_WIDTH_IN_SIGMA * sigma
    return np.array(
        [quad(_product_state_kernel, -bound, bound, args=(t, t0, i0, sigma, tau))[0] for t in delays]
    )


def build_transient_xas_map(
    ground_state_spectrum: npt.ArrayLike,
    excited_state_spectrum: npt.ArrayLike,
    delays: npt.ArrayLike,
    t0: float,
    i0: float,
    sigma: float,
    tau: float,
) -> np.ndarray:
    """
    Build a (energy x delay) transient (difference) XAS map for a simple
    two-state (ground -> excited) system: at each delay the transient is the
    excited-state population times the (excited - ground) difference
    spectrum.

    This is the pump-induced *change*, the quantity a pump-probe measurement
    reports; the full spectrum at a delay is the ground-state spectrum plus
    the corresponding column of this map.
    """
    population = excited_state_population(delays, t0, i0, sigma, tau)
    difference_spectrum = np.asarray(excited_state_spectrum) - np.asarray(ground_state_spectrum)
    return np.outer(difference_spectrum, population)


def build_transient_xes_map(
    singlet_spectrum: npt.ArrayLike,
    triplet_spectrum: npt.ArrayLike,
    quintet_spectrum: npt.ArrayLike,
    delays: npt.ArrayLike,
    t0: float,
    i0: float,
    sigma: float,
    tau: float,
) -> np.ndarray:
    """
    Build a (energy x delay) transient XES map for the singlet -> triplet ->
    quintet spin-crossover cascade: the triplet grows in and decays with
    lifetime ``tau`` (:func:`excited_state_population`), the quintet grows in
    as the triplet's complement (:func:`product_state_population`), and the
    singlet ground state is depleted by exactly the sum of the two (particle
    number conservation).
    """
    triplet_population = excited_state_population(delays, t0, i0, sigma, tau)
    quintet_population = product_state_population(delays, t0, i0, sigma, tau)
    return (
        np.outer(triplet_spectrum, triplet_population)
        + np.outer(quintet_spectrum, quintet_population)
        - np.outer(singlet_spectrum, triplet_population + quintet_population)
    )


def build_transient_xes_map_with_quintet_decay(
    singlet_spectrum: npt.ArrayLike,
    triplet_spectrum: npt.ArrayLike,
    quintet_spectrum: npt.ArrayLike,
    delays: npt.ArrayLike,
    t0: float,
    i0: float,
    sigma: float,
    tau: float,
    quintet_lifetime: float,
    cascade_window: float = 1000.0,
) -> np.ndarray:
    """
    Extend :func:`build_transient_xes_map` with the eventual decay of the
    long-lived quintet state back to the singlet ground state.

    Within ``cascade_window`` (in the same time unit as ``delays``, default
    picked for femtosecond delays) after ``t0``, the singlet/triplet/quintet
    cascade is evaluated exactly as in :func:`build_transient_xes_map`.
    Beyond that window, the transient spectrum is frozen at its value at the
    window edge and allowed to decay back towards zero (full ground-state
    recovery) with the much longer ``quintet_lifetime``.

    The slow decay is clocked from the cutoff, not from ``t0``, so the map is
    continuous across the window edge for any ratio of ``cascade_window`` to
    ``quintet_lifetime``.
    """
    delays = np.asarray(delays, dtype=float)
    singlet_spectrum = np.asarray(singlet_spectrum)
    triplet_spectrum = np.asarray(triplet_spectrum)
    quintet_spectrum = np.asarray(quintet_spectrum)

    cutoff = t0 + cascade_window
    early = delays < cutoff
    late = ~early

    transient = np.zeros((singlet_spectrum.size, delays.size))

    if early.any():
        triplet_population = excited_state_population(delays[early], t0, i0, sigma, tau)
        quintet_population = product_state_population(delays[early], t0, i0, sigma, tau)
        transient[:, early] = (
            np.outer(triplet_spectrum, triplet_population)
            + np.outer(quintet_spectrum, quintet_population)
            - np.outer(singlet_spectrum, triplet_population + quintet_population)
        )

    if late.any():
        # Freeze the cascade at the cutoff time itself -- rather than at
        # whichever grid point happens to fall just before it, so the result
        # doesn't depend on the delay grid's spacing -- then let it decay.
        triplet_at_cutoff = excited_state_population([cutoff], t0, i0, sigma, tau)[0]
        quintet_at_cutoff = product_state_population([cutoff], t0, i0, sigma, tau)[0]
        spectrum_at_cutoff = (
            triplet_at_cutoff * triplet_spectrum
            + quintet_at_cutoff * quintet_spectrum
            - (triplet_at_cutoff + quintet_at_cutoff) * singlet_spectrum
        )
        transient[:, late] = np.outer(
            spectrum_at_cutoff, np.exp(-(delays[late] - cutoff) / quintet_lifetime)
        )
    return transient
