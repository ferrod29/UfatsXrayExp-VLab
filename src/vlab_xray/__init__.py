"""
vlab_xray
=========

Analysis toolkit for the XAFS/XES virtual-lab (vLab) Praktikum.

Submodules
----------
constants        Physical constants and atomic masses.
io               Loaders for every file under ``data/``.
lineshapes       Peak-shape functions (Gaussian, Lorentzian, Voigt, ...).
fitting          Generic bounded curve-fitting helper.
xas_model        Static XAS (step + Gaussians) model and fit routine.
xes_model        Static XES (sum-of-4-Voigts) reference model.
kinetics         IRF-convolved population kinetics and transient map builders.
sample_estimates Sample/beam physics used in the sample-estimates exercise.
vlab_utils       Shared helpers for the Exp21 VLab exercise notebooks
                 (xraydb/CXRO wrappers, Bragg angles, loaders with
                 synthetic fallback).

See ``docs/`` at the repository root for the physical background behind
each module.
"""

from pathlib import Path

#: Repository ``data/`` directory, resolved from this package's own location
#: so it works regardless of the notebook's current working directory.
DATA_DIR = Path(__file__).resolve().parents[2] / "data"

__version__ = "0.1.0"

__all__ = ["DATA_DIR", "__version__"]
