# Ultrafast X-Ray Experiments — VLab (Exp21)

Notebooks and analysis toolkit for the advanced-lab experiment *Ultrafast X-Ray
Experiments* (script v2.0, July 2024): the virtual FXE / European XFEL
experiment measuring time-resolved Fe K$\beta_{1,3}$ X-ray emission of
$[\mathrm{Fe(bipy)_3}]^{2+}$ in aqueous solution after optical MLCT excitation,
plus the static and transient XAS/XES modeling of iron(II) spin-crossover
complexes that backs it up.

## Notebooks

Notebooks 01–04 cover all 35 exercises of **Section 5** of the script, one
section per exercise: a markdown statement (paraphrased from the script), the
computation, and a markdown result/interpretation. Notebooks 05–06 are the
modeling companions, working on measured reference data.

| Notebook | Section | Exercises | Content |
|---|---|---|---|
| [`01_beam_alignment.ipynb`](notebooks/01_beam_alignment.ipynb) | 5.1 | 1–9 | Choice of $E_0$; beam size/position/jitter (BIU2); pink vs. mono spectra & bandwidths (SpA1); IPM jitter; attenuator transmission; CRL focusing |
| [`02_preparatory_estimates.ipynb`](notebooks/02_preparatory_estimates.ipynb) | 5.2 | 10–21 | Molar masses; concentrations; X-ray attenuation budget (Fe vs. water); jet flow speed; pump concentration/focus/fluence & the linear-excitation limit; spectrometer & Bragg angles |
| [`03_experiment_acquisition.ipynb`](notebooks/03_experiment_acquisition.ipynb) | 5.3 | 22–29 | Jet centring; ground-state K$\beta$; pump–probe overlap; the 100 ps transient; figure-of-merit vs. position / intensity / delay |
| [`04_analysis_interpretation.ipynb`](notebooks/04_analysis_interpretation.ipynb) | 5.4 | 30–35 | Reference series & spin-state assignment; time-resolved difference spectra; false-colour map; kinetic traces & HS rise time; SVD + global sequential analysis for the ultrafast intermediate |
| [`05_xes_modeling.ipynb`](notebooks/05_xes_modeling.ipynb) | — | — | Static K$\beta$ XES reference spectra as sums of four Voigt profiles, and the transient singlet → triplet → quintet cascade |
| [`06_static_and_transient_xas_modeling.ipynb`](notebooks/06_static_and_transient_xas_modeling.ipynb) | — | — | Static Fe K-edge XANES (erf edge + near-edge resonances) and its pump-probe transient |

## Quick start

```bash
pip install -e ".[dev]"     # numpy, scipy, pandas, matplotlib, openpyxl, xraydb, molmass
pytest
jupyter lab notebooks/
```

The install is optional: each notebook's first cell puts this repository's
`src/` on `sys.path` and changes to the repository root, so the notebooks also
run from a clean checkout, and `data/` and `figures/` resolve the same way no
matter where Jupyter was started.

X-ray optical constants come from **`xraydb`** — the same tabulation as the CXRO
database referenced throughout the script — so the transmission/attenuation
exercises need **no network access**.

## Working with data

Every notebook runs **out of the box** on physically-motivated *synthetic*
stand-ins, so the analysis can be developed and tested before beamtime. To use
**real data**, export the corresponding trace/image from the VLab (or the online
app) and drop it into `data/` under the filename the loader expects — the
loaders prefer a real file whenever one is present, and return the path they
used so you can tell which is which. See [`data/README.md`](data/README.md) for
the full list of expected filenames per exercise, and
[`docs/data_dictionary.md`](docs/data_dictionary.md) for the measured data that
ships with the repository.

Online app (if accessible):
<https://play.unity.com/en/games/3e3b50d5-b027-473c-adc2-d606c5030baa/xfel-virtual-lab>

## Layout

```
.
├── notebooks/                        # the six notebooks
├── src/vlab_xray/                    # installable package
│   ├── constants.py                  #   physical constants and atomic masses
│   ├── io.py                         #   loaders for the measured data/ files
│   ├── lineshapes.py                 #   Gaussian, Lorentzian, Voigt, erf step
│   ├── fitting.py                    #   bounded curve-fitting helper
│   ├── xas_model.py                  #   static XAS model (step + Gaussians)
│   ├── xes_model.py                  #   static XES model (sum of 4 Voigts)
│   ├── kinetics.py                   #   IRF-convolved kinetics, transient maps
│   ├── sample_estimates.py           #   sample/beam planning physics
│   └── vlab_utils.py                 #   Exp21 helpers: xraydb/CXRO wrappers,
│                                     #   Bragg angles, loaders with synthetic
│                                     #   fallback, plotting style
├── data/                             # measured data + your VLab exports
├── docs/                             # physics background + architecture
├── tests/                            # pytest unit tests
└── figures/                          # PNGs exported by the notebooks
```

See [`docs/index.md`](docs/index.md) for the full documentation, starting with
[`docs/getting_started.md`](docs/getting_started.md).

## Notes

- **Short vs. long version.** The script offers a 3-day short version by skipping
  the grayed-out tasks. Every exercise is implemented here regardless; skip the
  cells you do not need.
- **Values to read from figures.** A few inputs must be taken from the script's
  figures rather than computed — chiefly the molar absorptivity $\epsilon(\lambda)$
  in Ex. 15 (Fig. 3.3). These are exposed as editable variables and flagged in a
  blockquote; replace the placeholder literature values with your own read-off.
- **Exercise 13** is left as `TODO` in the script (v2.0); the notebook records the
  optical-density relation $\mathrm{OD}=\epsilon c d$ it presumably introduces.
- **Synthetic data are teaching stand-ins.** Peak positions and relative
  intensities follow Figs. 4.3 / 4.4 and 3.6 / 3.7, but absolute numbers are
  illustrative. Replace them with VLab exports for a real evaluation.
