# Software Architecture

## Module map

```
src/vlab_xray/
├── constants.py          Physical constants & atomic masses (no internal deps)
├── lineshapes.py          Gaussian/Lorentzian/Voigt/erf-step peak shapes (no internal deps)
├── fitting.py             Generic bounded curve_fit wrapper (no internal deps)
├── io.py                  Loaders for every file in data/ (depends on: constants via DATA_DIR)
├── xas_model.py            Static XAS model + fit routine (depends on: lineshapes)
├── xes_model.py            Static XES reference model (depends on: lineshapes)
├── kinetics.py              IRF-convolved kinetics & transient map builders (no internal deps)
├── sample_estimates.py      Sample/beam planning formulas (depends on: constants)
└── vlab_utils.py            Exp21 exercise helpers (depends on: package DATA_DIR;
                             external: xraydb, molmass)
```

`lineshapes`, `fitting`, `kinetics`, and `constants` are pure math/physics
with no data dependency and no dependency on each other's siblings; `io`,
`xas_model`, `xes_model`, `sample_estimates`, and `vlab_utils` are the
"notebook-facing" layer that each notebook imports from.

The two notebook families reach the package differently. Notebooks 05-06 (the
modeling companions) import the individual model/kinetics modules and work on
the measured files in `data/`. Notebooks 01-04 (the 35 Section-5 exercises)
import `vlab_utils` alone, which bundles the `xraydb`/CXRO wrappers, Bragg
angles, plotting style and the loaders whose synthetic fallback lets an
exercise run before any data has been recorded.

`sample_estimates.py` overlaps with the planning arithmetic in notebook 02,
which does the same physics through `vlab_utils`. It is kept because the
sample-estimates notebooks it was written for live in the older
`Python_Scripts_VLab` checkout and were deliberately not carried over here; it
remains covered by `tests/test_sample_estimates.py`.

## Data flow

```
data/EXAFSterpy_noheader_APS-maybe.txt
        |  io.load_pump_probe_xas
        v
notebook 06  --fit_static_xas-->  data/FitResults_GS_ES.pkl
        |  xas_model.xas_model (ground/excited spectra)
        v
        kinetics.build_transient_xas_map  -->  transient XAS map (energy x delay)

data/Reference_Data.xlsx  +  data/Reference_spectra_rec_params.pkl
        |  io.load_xes_reference_spectra / load_xes_reference_fit_params
        v
notebook 05  --xes_model.rec_xes_4voigts-->  5 spin-state reference spectra
        |
        +--  data/Time_resolved_data.xlsx  (io.load_time_resolved_xes)
        v
        kinetics.build_transient_xes_map[_with_quintet_decay]  -->  transient XES map

data/cxro_*.txt, data/20140203_FeBPY_spectrum_Emilia.txt
        |  io.load_cxro_transmission / load_uvvis_spectrum
        v
        sample_estimates.*  -->  concentration / jet speed / focus size / excited fraction
        (notebook 02 computes the same quantities through vlab_utils, which
         reads the CXRO tabulation from xraydb rather than from data/)

your VLab exports in data/  (see data/README.md)
        |  vlab_utils.load_1d / load_image / load_delay_series
        v                        (synthetic fallback when absent)
notebooks 01-04  -->  beam alignment, planning estimates, acquisition, analysis
```

## Why a `src/` package instead of notebook-local functions

The original notebooks each redefined the same fitting/kinetics functions
with small, undocumented drifts between copies (kept in the `archive/` folder
of the older `Python_Scripts_VLab` checkout -- the `Model_XAS_for_VLAB.ipynb`
vs. `Modified_Model_XAS_for_VLAB.ipynb` step-edge discrepancy in particular).
Factoring the shared logic into an installable package means:

- there is exactly one implementation of e.g. the IRF-convolution kinetics,
  used identically by both the XAS and XES notebooks;
- it can be unit-tested independently of any notebook (`tests/`);
- notebooks stay focused on the narrative (what was measured, what was fit,
  what the result means) rather than utility-function boilerplate.

Every notebook starts with the same bootstrap cell: it walks up from the
current directory to the first ancestor containing `src/vlab_xray/`, prepends
that `src/` to `sys.path` and `chdir`s to the repository root -- see
[`getting_started.md`](../getting_started.md).

## Path resolution

`vlab_xray.DATA_DIR` (in `src/vlab_xray/__init__.py`) resolves the `data/`
directory relative to the package's own location, not the notebook's current
working directory. All `io.py` loaders, and `vlab_utils.DATA_DIR` which is
derived from it, look filenames up under `DATA_DIR`, so notebooks work
whether Jupyter's working directory is `notebooks/`, the repository root, or
anywhere else.

Two failure modes this guards against, both of which had actually bitten this
repository:

- *the wrong copy of the package.* An editable install (`pip install -e .`)
  of an older checkout elsewhere on the machine puts its own `src/` on
  `sys.path` via a `.pth` file. Prepending this repository's `src/` in the
  bootstrap cell, and setting `pythonpath = ["src"]` in `pyproject.toml` for
  pytest, keeps both notebooks and tests on these sources.
- *a `DATA_DIR` pointing at nothing.* The loaders in `vlab_utils` fall back to
  synthetic data when a file is missing, so a `DATA_DIR` resolved one
  directory too high or too low does not raise -- it silently yields synthetic
  results for a notebook the user believes is showing their measurement. The
  loaders return the path they used (`None` when synthetic) precisely so this
  is visible.
