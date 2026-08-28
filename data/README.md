# `data/`

Two kinds of file live here.

**Measured data shipped with the repository** — the Fe K-edge XAS scan, the
Kβ XES reference spectra, the time-resolved XES map and the CXRO transmission
curves used by notebooks 05–06. Each has a dedicated loader in
[`src/vlab_xray/io.py`](../src/vlab_xray/io.py); their columns, units and
provenance are documented in
[`docs/data_dictionary.md`](../docs/data_dictionary.md). Use the loader rather
than re-parsing the raw file.

**VLab exports you record yourself** — the exercise notebooks 01–04 read these
through `vlab_utils.load_1d` / `load_image` / `load_delay_series`, which look
for the file first and fall back to a physically-motivated synthetic generator
when it is absent. That is why every notebook runs end-to-end before beamtime:
drop a real export in with the expected name and the same cell switches to your
data, no code change needed. The loader returns the path it used (`None` for
synthetic), so you can always tell which one you are looking at.

## Expected filenames

| Notebook | Exercise | File | Loader | Content |
|---|---|---|---|---|
| 01 | 2 | `biu2_stack.npy` | `load_image` | Stack of BIU2 beam images (n × ny × nx) for position/jitter statistics |
| 01 | 3 | `biu2_pink.npy`, `biu2_mono.npy` | `load_image` | Single BIU2 beam image, pink and monochromatic beam |
| 01 | 4 | `spa1_pink_shot0.txt`, `…shot1`, `…shot2` | `load_1d` | Single-shot SASE spectra from SpA1 |
| 01 | 4 | `spa1_pink_avg.txt` | `load_1d` | Averaged pink-beam spectrum |
| 01 | 5 | `spa1_mono.txt` | `load_1d` | Monochromatised beam spectrum |
| 01 | 7 | `ipm.txt` | `np.loadtxt` | Four-diode IPM signals, one row per pulse: `I0 I1 I2 I3` |
| 01 | 9 | `xray_eye.npy` | `load_image` | X-ray eye image at the sample position |
| 03 | 22 | `jet_scan.txt` | `load_1d` | Signal vs. jet position, for centring the liquid jet |
| 03 | 23 | `kb_ground.txt` | `load_1d` | Ground-state Kβ emission spectrum |
| 03 | 24 | `xray_on_laser_off.npy`, `xray_off_laser_on.npy` | `load_image` | Spatial-overlap images for pump–probe alignment |
| 03 | 26 | `fom_xscan.txt`, `fom_yscan.txt` | `np.loadtxt` | Figure of merit vs. pump displacement (position [µm], FOM) |
| 03 | 27 | `fom_intensity.txt` | `np.loadtxt` | Figure of merit vs. pump pulse energy (energy [µJ], FOM) |
| 04 | 30 | `references.txt` | `np.loadtxt` | Kβ reference series: energy [eV], then singlet…quintet |
| 04 | 30–35 | `transients/` | `load_delay_series` | One 2-column file per delay, named `<delay_fs>.txt` (e.g. `-150.txt`, `100.txt`) |

The few files read with `np.loadtxt` rather than a `vlab_utils` loader are
read relative to the repository root, which every notebook's first cell
`chdir`s to. Exercise 28 has no file of its own: its FOM-vs-delay curve is
integrated from the Ex. 29/32 delay series.

## Formats

- **`load_1d`** — two whitespace- or comma-separated columns, *x* then *y*.
  Extensions tried in order: none, `.txt`, `.csv`, `.dat`. Leading text
  headers (up to four lines) are skipped automatically.
- **`load_image`** — `.npy` array, or a text file holding the pixel matrix.
  Extensions tried: none, `.npy`, `.txt`, `.csv`.
- **`load_delay_series`** — a subdirectory holding one `load_1d`-style file per
  pump–probe delay, the filename being the delay in femtoseconds. Returns
  `(energy, delays_fs, matrix, path)` sorted by delay.

`load_1d` and `load_image` likewise return the path they read as their last
element (`None` when the synthetic generator was used).

Energies are in eV and delays in fs throughout, matching the conventions in
[`src/vlab_xray/vlab_utils.py`](../src/vlab_xray/vlab_utils.py).
