# vlab_xray documentation

Analysis toolkit and Praktikum notebooks for the virtual lab (vLab) on
femtosecond pump-probe X-ray absorption (XAS/XANES) and X-ray emission (XES)
spectroscopy of iron(II) spin-crossover complexes.

## Start here

- [`getting_started.md`](getting_started.md) -- install, run the notebooks, run the tests.
- [`data_dictionary.md`](data_dictionary.md) -- what every file in `data/` contains.
- [`software/architecture.md`](software/architecture.md) -- module map and data flow.

## Physical background

1. [`physics/01_xas_xanes.md`](physics/01_xas_xanes.md) -- X-ray absorption, XANES edge + resonance structure, why the edge is modeled as an erf.
2. [`physics/02_xes_spin_states.md`](physics/02_xes_spin_states.md) -- X-ray emission as a spin-state marker, why 4 Voigt profiles per reference spectrum.
3. [`physics/03_ultrafast_pump_probe_kinetics.md`](physics/03_ultrafast_pump_probe_kinetics.md) -- IRF convolution, the singlet->triplet->quintet cascade, the quintet's long-lifetime decay.
4. [`physics/04_sample_and_beam_considerations.md`](physics/04_sample_and_beam_considerations.md) -- concentration, optical density, liquid jet, laser/X-ray focal overlap, excited-state fraction.

## The notebooks

| Notebook | Covers |
| --- | --- |
| `01_beam_alignment.ipynb` | Section 5.1, Ex. 1-9: beam energy, size, position, jitter, bandwidth, focusing |
| `02_preparatory_estimates.ipynb` | Section 5.2, Ex. 10-21: sample/beam planning arithmetic (physics doc 4) |
| `03_experiment_acquisition.ipynb` | Section 5.3, Ex. 22-29: jet centring, pump-probe overlap, figures of merit |
| `04_analysis_interpretation.ipynb` | Section 5.4, Ex. 30-35: difference spectra, kinetic traces, SVD + global analysis |
| `05_xes_modeling.ipynb` | Static + transient Kbeta XES (physics docs 2, 3) |
| `06_static_and_transient_xas_modeling.ipynb` | Static + transient Fe K-edge XAS (physics docs 1, 3) |
