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
through `vlab_utils.load_1d` / `load_image` / `load_table` /
`load_delay_series`, which look for the file first and fall back to a
physically-motivated synthetic generator when it is absent. That is why every notebook runs end-to-end before beamtime:
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

## `raw/` and `processed/` — your own VLab sessions

Both directories are **git-ignored**: they hold measurement data, not source.
A session runs to a few hundred megabytes, so it stays on your machine and is
rebuilt from the raw exports rather than versioned.

Each session lives in its own dated folder:

```
data/raw/20260901/...          exports straight out of the VLab
data/processed/20260901/*.npy  the arrays built from them
```

Build them with

```bash
python scripts/build_processed_data.py            # the newest session
python scripts/build_processed_data.py 20260901   # a specific one
python scripts/build_processed_data.py --all      # every session
```

The loaders search the newest built session, then `data/` itself, and accept
`.npy` as well as text, so the notebooks pick the arrays up with no edits.
To analyse an older session instead:

```python
V.use_session("20260226")     # or set VLAB_SESSION=20260226 in the environment
```

Session folders are ordered by date, not by name — one session is named
`20263108`, which is 31 August written `YYYYDDMM` and would otherwise sort
after September.

### Array names

| Array | Shape | Exercise | Built from |
|---|---|---|---|
| `biu2_stack.npy` | (n, 420, 420) | 2 | one imager frame per pulse |
| `biu2_pink.npy`, `biu2_mono.npy` | (420, 420) | 3 | a full-beam and a monochromatised frame |
| `spa1_pink_shot{0,1,2}.npy`, `spa1_pink_avg.npy` | (420, 2) | 4 | the SpA1 pink shots |
| `spa1_mono.npy` | (420, 2) | 5 | the SpA1 monochromatic shots |
| `spa1_*_unattenuated.npy` | (420, 2) | — | any ~100× brighter shot in a set |
| `ipm.npy`, `ipm_mono.npy`, `ipm_scan_position.npy` | (n, 4), (n,) | 7 | the IPM export, else the imager frame headers |
| `xray_eye.npy` | (420, 420) | 9 | the focused spot at the sample |
| `xray_eye_unfocused_stack.npy` | (n, 420, 420) | 9 | the same detector before focusing |
| `biu3_focus_{tight,wide}_stack.npy` | (n, 420, 420) | 9 | BIU3 at two CRL settings |
| `tad_{pink,mono}.npy` | (420, 2) | — | the arrival-time monitor's step edge |
| `jet_scan.npy` | (420, 2) | 22 | the LPD scan across the jet |
| `kb_ground.npy`, `kb_laser_off_null.npy` | (420, 2) | 23 | laser-off von Hamos: the one with signal, and the null |
| `transients_{energy,delays,matrix}.npy` | (420,), (n,), (n, 420) | 29, 32–35 | the delay scan, repeats averaged |
| `transients_{sem,n_repeats}.npy` | (n, 420), (n,) | 29, 32–35 | spread across those repeats |

Not every session holds every measurement; the build prints what it skipped,
and the notebooks fall back to synthetic stand-ins for the rest.

### What the build works out for itself

**Which folder is the pink beam.** Never taken from the folder name. A
monochromator can neither broaden its input nor add flux, so the pink beam is
the set that is spectrally wider, varies shot to shot, and carries more
photons; `classify_beam` checks all three and prints what it found. This
matters: in the **20260226** session the two folders are *swapped* — the set
in `pink_X-ray/` is narrow (ΔE/E = 3.2e-4), identical shot to shot
(correlation 1.000) and dim, while `monochromated_x-ray/` is broad
(2.4–4.5e-3), spiky (0.52) and 39× brighter. The **20260901** session's labels
are correct, and the build says so.

**A changed spectrometer range.** If the delay scan switched range part-way
through, the delays sit on different energy grids; they are interpolated onto
the overlap and the build reports it.

Two things the exports do not record, so the script assumes them — change them
at the top of `scripts/build_processed_data.py` if your session differed:

- **Delay units.** Delay folder names carry no unit and are read as
  femtoseconds (`DELAY_UNIT_FS`).
- **Attenuator outliers.** An acquisition integrating ~100× above the rest of
  its set is kept aside rather than averaged in (`OUTLIER_FACTOR`).
