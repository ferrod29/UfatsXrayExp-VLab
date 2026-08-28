# X-ray Absorption (XAS/XANES) of the Fe K-edge

Code: [`src/vlab_xray/xas_model.py`](../../src/vlab_xray/xas_model.py),
[`src/vlab_xray/lineshapes.py`](../../src/vlab_xray/lineshapes.py) ·
Notebook: [`notebooks/06_static_and_transient_xas_modeling.ipynb`](../../notebooks/06_static_and_transient_xas_modeling.ipynb)

## 1. Why X-rays absorb: the photoelectric effect and Beer-Lambert law

When an X-ray photon has enough energy to eject a core electron from an atom
(photoelectric absorption), the transmitted intensity through a sample of
thickness `l` follows the Beer-Lambert law:

```
I(E) = I0 * exp(-mu(E) * l)
```

where `mu(E)` is the energy-dependent linear absorption coefficient. `mu(E)`
rises sharply whenever the photon energy crosses the binding energy of a core
level -- an **absorption edge**. The data in
`data/EXAFSterpy_noheader_APS-maybe.txt` scans across the iron **K-edge**
(1s binding energy, ~7112 eV for metallic Fe, shifted somewhat by chemical
environment) of the spin-crossover complex [Fe(terpy)<sub>2</sub>]<sup>2+</sup>.

## 2. Anatomy of a XANES spectrum: edge jump + near-edge resonances

X-ray Absorption Near-Edge Structure (XANES) spans roughly the 50 eV around
an edge. Two physically distinct contributions dominate it:

1. **The edge jump itself**: above the binding energy, the core electron is
   promoted to the continuum (ionized). To first approximation this is a
   step function in energy -- absorption is "off" below the edge and "on"
   above it.
2. **Pre-edge and near-edge resonances**: dipole-allowed transitions from the
   1s core level into *bound* valence and low-lying continuum states (for a
   3d transition metal, mainly 1s -> 3d (formally weak, quadrupole/dipole
   mixed via 3d-4p hybridization) and 1s -> 4p transitions) appear as
   discrete peaks superimposed on/near the edge, below and just above the
   ionization threshold.

`xas_model.xas_model` captures exactly this structure:

```python
model(x) = erf_step(x; edge_position, edge_height, edge_width)
           + sum_i amplitude_i * gaussian(x; center_i, width_i)
           + offset
```

one broadened step for the edge jump, and a handful of Gaussians for the
resonances (their number and starting positions come from prior knowledge of
where the pre-edge/near-edge features of this complex lie).

## 3. Why the edge is an error function (`erf`), not a raw step

An ideal absorption edge is a Heaviside step. Any real measurement, however,
convolves the "true" spectrum with the instrument's resolution function
(monochromator bandwidth, core-hole lifetime broadening, etc.), which to good
approximation is Gaussian. The convolution of a step function with a Gaussian
of standard deviation `sigma` is *exactly* an error function:

```
step(x) convolved with Gaussian(sigma)  =  0.5 * (1 + erf((x - x0) / (sigma*sqrt(2))))
```

`lineshapes.erf_step` implements the `erf` part of this as
`height * erf((x - x0) / (sigma*sqrt(2)))`: the additive constant is absorbed
into the model's `offset` term, and the factor of one half into `height`.
Two consequences worth keeping in mind when reading fitted parameters: the
curve runs from `-height` to `+height`, so the **total edge jump is
`2 * edge_height`**, and the fitted `edge_height` in
`data/FitResults_GS_ES.pkl` is therefore half the jump.

This erf form is also why `Modified_Model_XAS_for_VLAB.ipynb` (in the
`archive/` folder of the older `Python_Scripts_VLab` checkout, not carried
over here) switched away from a logistic sigmoid -- a convenient but ad hoc
approximation to a broadened step: the erf form is the analytically correct
broadening, not an approximation, and is the version this project's
`xas_model.py` uses.

## 4. Ground state, excited state, and the pump-probe difference signal

The measured file has three signal columns:

- `OFF`: X-rays only, no laser pump -- the **ground-state** spectrum.
- `ON`: X-rays with the laser pump at some fixed pump-probe delay -- a
  *mixture* of ground- and excited-state molecules (only a fraction of the
  sample absorbs a pump photon).
- `Diff = ON - OFF`: the pump-induced difference signal, small because only
  a fraction `excited_state_fraction` of the sample is excited.

To recover the *pure* excited-state spectrum (as if 100% of the sample were
excited), `io.load_pump_probe_xas` extrapolates:

```
ES = OFF + Diff / excited_state_fraction
```

This assumes the measured signal is a linear combination
`ON = (1 - f)*OFF + f*ES` of the two pure spectra with excited fraction `f`,
which rearranges to the formula above. `xas_model.fit_static_xas` is then
fit independently to `OFF` and `ES` to obtain the two reference spectra
`params_opt` (ground state) and `params_opt_es` (excited state) stored in
`data/FitResults_GS_ES.pkl` -- these are the two "pure" spectra that the
transient kinetics model (see
[`03_ultrafast_pump_probe_kinetics.md`](03_ultrafast_pump_probe_kinetics.md))
interpolates between at each pump-probe delay.
