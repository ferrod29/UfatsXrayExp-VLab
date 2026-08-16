# Getting Started

## Install

From the repository root, in a Python >= 3.10 environment:

```bash
pip install -e ".[dev]"
```

This installs `vlab_xray` (from `src/`) in editable mode, plus `pytest` for
the test suite. `-e` means edits to `src/vlab_xray/*.py` take effect
immediately, with no reinstall step.

You don't strictly have to install the package to run the notebooks: every
notebook's first code cell walks up from the current directory until it finds
`src/vlab_xray/`, puts that `src/` at the front of `sys.path` and changes to
the repository root -- so opening a notebook directly from `notebooks/` and
running it top-to-bottom works out of the box, and `data/` and `figures/`
resolve identically however Jupyter was started. Prepending rather than
appending matters: an editable install of another copy of `vlab_xray` (an
older checkout elsewhere on the machine) must not shadow these sources.

## Run the notebooks

```bash
jupyter lab notebooks/
```

Run each notebook top-to-bottom:

1. `01_beam_alignment.ipynb` -- Section 5.1 of the script (Ex. 1-9).
2. `02_preparatory_estimates.ipynb` -- Section 5.2 (Ex. 10-21), the
   sample/beam planning calculations.
3. `03_experiment_acquisition.ipynb` -- Section 5.3 (Ex. 22-29).
4. `04_analysis_interpretation.ipynb` -- Section 5.4 (Ex. 30-35).
5. `05_xes_modeling.ipynb` -- static + transient Kbeta XES modeling.
6. `06_static_and_transient_xas_modeling.ipynb` -- static + transient Fe
   K-edge XAS modeling.

Notebooks 01-04 run on synthetic stand-ins until you drop real VLab exports
into `data/` (see [`../data/README.md`](../data/README.md)); 05-06 work on the
measured data that ships with the repository.

## Run the tests

```bash
pytest
```

Tests in `tests/` exercise the package's lineshape functions, kinetics
convolution integrals, data loaders (against the real files in `data/`), and
sample-estimate formulas.

## Where things live

- `src/vlab_xray/` -- the installable package (all reusable code).
- `notebooks/` -- the six notebooks, importing from `vlab_xray`.
- `data/` -- measured data plus your own VLab exports (see
  [`data_dictionary.md`](data_dictionary.md) and
  [`../data/README.md`](../data/README.md)).
- `docs/` -- this documentation: physics background (`docs/physics/`) and
  software architecture (`docs/software/`).
- `tests/` -- pytest unit tests for `src/vlab_xray/`.
- `figures/` -- PNGs exported by the notebooks.
