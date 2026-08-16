# X-ray Emission Spectroscopy (XES) as a Spin-State Marker

Code: [`src/vlab_xray/xes_model.py`](../../src/vlab_xray/xes_model.py),
[`src/vlab_xray/lineshapes.py`](../../src/vlab_xray/lineshapes.py) ·
Notebook: [`notebooks/03_xes_modeling.ipynb`](../../notebooks/03_xes_modeling.ipynb)

## 1. What XES measures

X-ray *absorption* promotes a core electron out of the atom, leaving a
core hole. X-ray *emission* is what happens next: a higher-shell electron
(here, a 3p electron) falls into the 1s core hole, emitting a photon whose
energy equals the difference between the two orbital energies. This
particular 3p -> 1s transition is called **K&beta;** emission (the K&beta;
"mainline" studied here sits around 7040-7070 eV for iron, matching the
energy axis of `data/Reference_Data.xlsx` and
`data/Time_resolved_data.xlsx`).

## 2. Why K&beta; is spin-sensitive: the 3p-3d exchange interaction

The emitting 3p electron and the valence 3d electrons of the same atom are
close enough in space to interact via **exchange coupling**. The strength of
this 3p-3d exchange interaction depends on the *number of unpaired 3d
electrons*, i.e. on the total spin `S` of the ion. A higher-spin Fe(II)
configuration (more unpaired 3d electrons) produces a larger exchange
splitting of the final (3p hole) state, which broadens and shifts the K&beta;
emission line-shape -- mainly visible as an increasingly pronounced
low-energy shoulder (historically called K&beta;'). This makes K&beta; XES a
direct, element-specific, **local spin-state marker** -- exactly why this
Praktikum tracks it through a spin-crossover reaction.

## 3. The five spin-state references and the spin-crossover cascade

`data/Reference_Data.xlsx` contains five reference emission spectra, labeled
by their total-spin multiplicity: `singlet`, `doublet`, `triplet`, `quartet`,
`quintet` (2S+1 = 1, 2, 3, 4, 5). For the ferrous ([Fe(bpy)<sub>3</sub>]<sup>2+</sup>-type)
spin-crossover complexes studied here:

- **Singlet** (S=0, low-spin, t2g^6): the thermally stable ground state.
- **Quintet** (S=2, high-spin, t2g^4 eg^2): the long-lived photo-excited
  product state.
- **Triplet** and **quartet/doublet** intermediates: transient states
  visited during the ultrafast intersystem-crossing cascade between singlet
  and quintet (the well-known femtosecond-to-picosecond
  ¹A -> ³T -> ⁵T "spin-crossover cascade" of Fe(II) polypyridyl complexes).

`kinetics.build_transient_xes_map` models the pump-probe signal as a
population-weighted mixture of exactly these reference spectra (see
[`03_ultrafast_pump_probe_kinetics.md`](03_ultrafast_pump_probe_kinetics.md)).

## 4. Why 4 Voigt profiles per reference spectrum

Each reference spectrum is modeled as a sum of 4 Voigt profiles
(`xes_model.rec_xes_4voigts`, using `lineshapes.voigt`). A Voigt profile is
the convolution of:

- a **Lorentzian**, with a natural linewidth set by the core-hole and
  valence-hole lifetimes (energy-time uncertainty: a state with lifetime
  `tau` has an intrinsic linewidth `~ hbar/tau`), and
- a **Gaussian**, from instrumental/spectrometer broadening.

Four components are used because the K&beta; mainline of a 3d transition
metal is not a single transition: multiplet structure and spin-orbit/exchange
splitting of the final 3p⁵3d^n state produce several overlapping,
unresolved sub-peaks, which four Voigt profiles are enough to reproduce
empirically (the `I0, x0, gamma, sigma` parameters for all 4 components x 5
spin states are pre-fit and stored in
`data/Reference_spectra_rec_params.pkl`, loaded by
`io.load_xes_reference_fit_params`).
