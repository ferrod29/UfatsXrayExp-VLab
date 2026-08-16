"""
Generic bounded curve-fitting helper.

A thin, robust wrapper around ``scipy.optimize.curve_fit`` used whenever a
lineshape from ``lineshapes`` needs to be fit to data (e.g. reproducing the
Voigt parameters stored in ``data/Reference_spectra_rec_params.pkl``). Static
and transient XAS/XES fitting in ``xas_model`` uses
``scipy.optimize.minimize`` directly instead, because it needs a masked,
custom sum-of-squares objective (see ``xas_model.fit_static_xas``).
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import scipy.optimize as optimize


def fit(f: Callable, x: Sequence[float], y: Sequence[float], p0: Sequence[float]):
    """
    Fit ``f(x, *params)`` to ``(x, y)`` starting from ``p0``.

    Uses the Trust Region Reflective algorithm with a soft-L1 loss, which is
    robust to the occasional outlier point in experimental spectra without
    needing per-point weights.

    Returns
    -------
    (popt, pcov) : optimal parameters and their covariance matrix, as
        returned by ``scipy.optimize.curve_fit``.
    """
    popt, pcov = optimize.curve_fit(
        f,
        x,
        y,
        p0=p0,
        method="trf",
        tr_solver="lsmr",
        loss="soft_l1",
        jac="3-point",
        max_nfev=1_000_000,
    )
    return popt, pcov
