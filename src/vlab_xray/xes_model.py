"""
Static XES (X-ray emission) reference model: a sum of 4 Voigt profiles.

Each measured Kbeta emission reference spectrum (one per spin state:
singlet, doublet, triplet, quartet, quintet) is modeled as the sum of 4
Voigt profiles (``lineshapes.voigt``), reflecting 4 unresolved
multiplet/spin-orbit components under the Kbeta mainline. See
``docs/physics/02_xes_spin_states.md``.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import numpy.typing as npt

from vlab_xray.lineshapes import voigt


def rec_xes_4voigts(x: npt.ArrayLike, params: Sequence[float]) -> np.ndarray:
    """
    Evaluate the sum-of-4-Voigts XES reference model.

    ``params`` is 16 numbers: 4 consecutive groups of
    ``(amplitude, center, gamma, sigma)``, one group per Voigt component.
    Matches the layout stored in
    ``data/Reference_spectra_rec_params.pkl`` (loaded with
    ``io.load_xes_reference_fit_params``).
    """
    p1, p2, p3, p4 = (params[i : i + 4] for i in range(0, 16, 4))
    return voigt(x, *p1) + voigt(x, *p2) + voigt(x, *p3) + voigt(x, *p4)
