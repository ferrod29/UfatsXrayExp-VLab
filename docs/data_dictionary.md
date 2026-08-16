# Data Dictionary

All paths below are relative to `data/`. Every file has a matching loader in
[`src/vlab_xray/io.py`](../src/vlab_xray/io.py) that returns it as a tidy,
labeled `pandas` object -- use the loader rather than re-parsing the raw
file.

## `EXAFSterpy_noheader_APS-maybe.txt`

Static/pump-probe Fe K-edge XAS scan of [Fe(terpy)<sub>2</sub>]<sup>2+</sup>,
whitespace-separated, no header, 4 columns: `Energy` (keV), `ON` (laser-on),
`OFF` (laser-off / ground state), `Diff` (`ON - OFF`). Measured at an
APS synchrotron beamline (per filename; exact beamline/run details are not
otherwise recorded). Loader: `io.load_pump_probe_xas` (also adds an `ES`
extrapolated-excited-state column -- see
[`physics/01_xas_xanes.md`](physics/01_xas_xanes.md)).

## `cxro_water_25u.txt`, `cxro_iron_00595u.txt`

X-ray transmission curves for a 25&micro;m water layer and a 0.0595&micro;m
iron foil respectively, over 7000-7200 eV (2-line header, then `Energy`
[eV], `Transmission` [0-1]). Computed from the CXRO/Henke atomic scattering
factor tables (https://henke.lbl.gov/optical_constants/filter2.html), used
to compare the iron K-edge absorption signal against the water background --
see [`physics/04_sample_and_beam_considerations.md`](physics/04_sample_and_beam_considerations.md).
Loader: `io.load_cxro_transmission`.

## `20140203_FeBPY_spectrum_Emilia.txt`

UV-vis absorption spectrum of the [Fe(bpy)<sub>3</sub>]<sup>2+</sup> sample
(2-line header, then `Wavelength` [nm], raw instrument `Abs.` reading).
Loader: `io.load_uvvis_spectrum` multiplies the raw reading by a fixed
calibration factor (5e4) to convert it to molar absorptivity `epsilon`
[L mol<sup>-1</sup> cm<sup>-1</sup>].

## `Reference_Data.xlsx`

Measured static Kbeta XES reference spectra for the five spin states of the
spin-crossover cascade. Columns: `emission energy` (eV, becomes the index),
`singlet`, `doublet`, `triplet`, `quartet`, `quintet` (intensity, arb.
units). See [`physics/02_xes_spin_states.md`](physics/02_xes_spin_states.md).
Loader: `io.load_xes_reference_spectra`.

## `Reference_spectra_rec_params.pkl`

Pickled `pandas.DataFrame`: pre-fit 4-Voigt parameters reproducing
`Reference_Data.xlsx`. 16 rows (`I0_i, x0_i, gamma_i, sigma_i` for
`i in 1..4`), one column per spin state. Consumed by
`xes_model.rec_xes_4voigts`. Loader: `io.load_xes_reference_fit_params`.
These parameters were obtained from a fit performed in a separate script,
not included in this project.

## `Time_resolved_data.xlsx`

Time-resolved (pump-probe) XES difference map. First column: emission energy
(eV). Remaining columns: pump-probe delay, given in the raw file in
picoseconds; `io.load_time_resolved_xes` converts these to femtoseconds and
relabels the columns (e.g. `"-567 fs"`) accordingly. Cell values are the
(laser-on minus laser-off) difference signal at that energy/delay.

## `FitResults_GS_ES.pkl`

Pickled `pandas.DataFrame` with two rows, `Ground State` and
`Excited State`: the fitted `xas_model.xas_model` parameters
(`[edge_x0, edge_height, edge_sigma, (amp, x0, sigma) * N, offset]`)
reproducing the `OFF` and `ES` columns of
`EXAFSterpy_noheader_APS-maybe.txt`. The `Ground State` row has fewer
resonance Gaussians than `Excited State` (an extra pre-edge feature appears
upon excitation), so it is padded with trailing `NaN`s to match column
count -- always `.dropna()` after indexing a row before passing it to
`xas_model.xas_model`. Loader: `io.load_fit_params`.
