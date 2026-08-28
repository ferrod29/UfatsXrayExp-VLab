# Sample and Beam Considerations for a Pump-Probe X-ray Experiment

Code: [`src/vlab_xray/sample_estimates.py`](../../src/vlab_xray/sample_estimates.py),
[`src/vlab_xray/constants.py`](../../src/vlab_xray/constants.py) ·
Notebook: [`notebooks/02_preparatory_estimates.ipynb`](../../notebooks/02_preparatory_estimates.ipynb)
(Ex. 10-19, which works the same physics through `vlab_utils`)

Before running a femtosecond pump-probe X-ray experiment, a set of practical
questions need quantitative answers: how much sample is needed, how
concentrated does it need to be, how fast must it flow, and how should the
laser and X-ray pulses be matched. This is the physics behind those
estimates, worked out for the [Fe(bpy)<sub>3</sub>]<sup>2+</sup> spin-crossover
complex.

## 1. Molar mass and concentration

Straightforward stoichiometry (`sample_estimates.molar_mass`,
`sample_estimates.concentration`): sum atomic masses over the molecular
formula, then `concentration = (mass / molar_mass) / volume`. The sample
salt, [Fe(bpy)<sub>3</sub>]Cl<sub>2</sub>&middot;3H<sub>2</sub>O, crystallizes
with three waters of hydration, which must be included in the molar mass
used to convert a weighed-out mass to a solution concentration.

## 2. Beer-Lambert law, optical density, and X-ray transmission

The same Beer-Lambert law as in [`01_xas_xanes.md`](01_xas_xanes.md) governs
how much of the beam a sample of a given concentration and thickness
absorbs:

```
OD = -log10(I/I0) = epsilon * c * l      (optical density form)
I/I0 = exp(-mu * l)                      (X-ray transmission form)
```

`data/cxro_water_25u.txt` and `data/cxro_iron_00595u.txt` are pre-computed
X-ray transmission curves (from the CXRO/Henke optical-constants database,
https://henke.lbl.gov/optical_constants/filter2.html) for a 25&mu;m water jet
and a fictive iron foil, in the same 7000-7200 eV window as the Fe K-edge XAS
scan. **Photoelectric absorption dominates both**: at 7-7.2 keV it accounts
for ~97% of water's total attenuation (the rest being ~2% coherent/Rayleigh
and ~1% incoherent/Compton scattering -- Compton only takes over for light
elements well above ~30 keV). What distinguishes the two is not the
mechanism but the *energy dependence*: water is far above the O K-edge, so
its attenuation is a smooth, featureless `E^-3` background, whereas iron's
cross-section jumps at 7112 eV and carries the near-edge structure being
measured. Comparing the curves is what shows the sample must be thin enough,
and concentrated enough, that the dilute iron's edge step is not swamped by
that flat water background.

For the optical (UV-vis) excitation side, `sample_estimates.beer_lambert_concentration`
inverts the optical-density form to answer "what concentration gives OD = 1
at this laser wavelength, through this sample thickness", using the sample's
measured molar absorptivity spectrum (`data/20140203_FeBPY_spectrum_Emilia.txt`).

## 3. Liquid jet sample delivery

Pump-probe X-ray experiments typically flow the sample as a thin free-falling
liquid jet, so that every laser/X-ray pulse pair hits **fresh, unphotolyzed
sample**. The minimum jet flow speed is set by requiring the jet to advance
at least one beam-spot-worth of thickness between consecutive pulses:

```
v_min = jet_thickness / pulse_period = jet_thickness * repetition_rate
```

(`sample_estimates.jet_flow_speed_min`) -- if the jet were slower, the next
pulse would probe sample that the previous pulse (and its heat/photoproduct)
had already affected.

## 4. Laser focus size and laser/X-ray spatial overlap

The diffraction-limited laser focus diameter follows the Rayleigh criterion,

```
d = 2.44 * lambda * f / D
```

(`sample_estimates.rayleigh_focus_diameter`, with `f` the focusing lens'
focal length and `D` its aperture diameter). The factor is 2.44 because this
is the *diameter* of the Airy disc out to its first zero; the more familiar
`1.22 * lambda * f / D` is its radius. For a clean pump-probe
measurement, the laser spot should be significantly *larger* than the X-ray
spot: the X-ray probe then only ever samples a sub-region of the pumped
volume where the excitation density is close to uniform, avoiding a mixture
of differently-excited (or entirely unexcited) volumes in a single
measurement.

## 5. Excited-state fraction and photon-vs-molecule counting

Two related numbers matter for planning the pump pulse energy: how many
sample molecules sit in the X-ray-probed, laser-illuminated volume
(`sample_estimates.molecules_in_volume`, via Avogadro's number), and what
pulse energy would deliver exactly one photon per molecule
(`sample_estimates.pulse_energy_for_molecule_count`, via the photon energy
`E = hc/lambda`, `sample_estimates.photon_energy`). In practice only a
fraction of molecules that see a photon actually absorb and photoreact (see
`excited_state_fraction` in
[`01_xas_xanes.md`](01_xas_xanes.md), fixed at 0.76 for the measured
[Fe(terpy)<sub>2</sub>]<sup>2+</sup> data set), so the "one photon per
molecule" pulse energy is a lower bound, not a guarantee of full conversion.
