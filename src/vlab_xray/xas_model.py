"""
Static XAS (X-ray absorption) model: an edge jump plus near-edge resonances.

The model represents a XANES spectrum as:

    model(x) = erf_step(x; edge position, edge height, edge width)
               + sum_i  amplitude_i * gaussian(x; center_i, width_i)
               + constant offset

i.e. one Gaussian-broadened absorption edge (see ``lineshapes.erf_step``)
plus several Gaussian pre-edge/near-edge resonance peaks. Separate fits of
this model to the laser-off (ground state) and laser-on-derived (excited
state) spectra give the two parameter sets stored in
``data/FitResults_GS_ES.pkl``. See ``docs/physics/01_xas_xanes.md``.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import numpy.typing as npt
import scipy.optimize as optimize

from vlab_xray.lineshapes import erf_step, gaussian


def xas_model(x: npt.ArrayLike, params: Sequence[float]) -> np.ndarray:
    """
    Evaluate the step-plus-Gaussians XAS model.

    ``params`` layout: ``[edge_x0, edge_height, edge_sigma,
    (amp_1, x0_1, sigma_1), (amp_2, x0_2, sigma_2), ..., offset]`` -- i.e. 3
    edge parameters, then triplets of Gaussian parameters for as many peaks
    as needed, then a single constant offset.
    """
    x = np.asarray(x, dtype=float)
    y = erf_step(x, params[0], params[1], params[2])
    for i in range(3, len(params) - 1, 3):
        y = y + params[i] * gaussian(x, params[i + 1], params[i + 2])
    y = y + params[-1]
    return y


def _sum_squared_error(params: Sequence[float], x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> float:
    return np.sum((xas_model(x[mask], params) - y[mask]) ** 2)


def fit_static_xas(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    initial_params: Sequence[float],
    bounds: Sequence[tuple[float, float]],
    fit_range: tuple[float, float] | None = None,
):
    """
    Fit the step-plus-Gaussians XAS model to ``(x, y)`` by bounded
    least-squares (``scipy.optimize.minimize``), restricted to
    ``fit_range = (x_min, x_max)`` if given.

    Returns the ``scipy.optimize.OptimizeResult``; the fitted parameters are
    ``result.x`` (same layout as ``xas_model``'s ``params``).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if fit_range is None:
        mask = np.ones_like(x, dtype=bool)
    else:
        mask = (x > fit_range[0]) & (x < fit_range[1])
    return optimize.minimize(
        _sum_squared_error,
        initial_params,
        args=(x, y, mask),
        bounds=bounds,
    )
