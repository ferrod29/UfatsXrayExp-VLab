"""
Loaders for every data file under ``data/``.

Each function returns a tidy ``pandas`` object with descriptive column
names/index, so notebooks never need to re-derive "which column is which"
or re-parse a file format. See ``docs/data_dictionary.md`` for the
provenance and units of each file.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from vlab_xray import DATA_DIR


def _resolve(path: str | Path) -> Path:
    """Resolve ``path`` against DATA_DIR unless it is already absolute/exists as given."""
    path = Path(path)
    if path.exists() or path.is_absolute():
        return path
    return DATA_DIR / path


def load_pump_probe_xas(
    path: str | Path = "EXAFSterpy_noheader_APS-maybe.txt",
    excited_state_fraction: float = 0.76,
) -> pd.DataFrame:
    """
    Load the static-laser-on/laser-off XAS pump-probe scan of
    [Fe(terpy)2]2+ (whitespace-separated, no header).

    Columns: ``Energy`` [keV], ``ON`` (laser-on transmission), ``OFF``
    (laser-off / ground-state transmission), ``Diff`` (ON - OFF). An ``ES``
    (excited-state) column is added, extrapolating the difference signal to
    a 100%-excited sample: ``ES = OFF + Diff / excited_state_fraction`` (see
    ``docs/physics/01_xas_xanes.md``).
    """
    data = pd.read_csv(_resolve(path), sep=r"\s+", header=None)
    data.columns = ["Energy", "ON", "OFF", "Diff"]
    data["ES"] = data["OFF"] + data["Diff"] / excited_state_fraction
    return data


def load_cxro_transmission(path: str | Path) -> pd.DataFrame:
    """
    Load a CXRO/Henke-tables filter transmission curve (``cxro_*.txt``).

    Columns: ``Energy`` [eV], ``Transmission`` (0-1).
    """
    return pd.read_csv(
        _resolve(path),
        sep=r"\s+",
        names=["Energy", "Transmission"],
        skiprows=2,
    )


def load_uvvis_spectrum(
    path: str | Path = "20140203_FeBPY_spectrum_Emilia.txt",
    calibration_factor: float = 5e4,
) -> pd.DataFrame:
    """
    Load a UV-vis absorption spectrum and calibrate it to molar absorptivity.

    Columns: ``Wavelength`` [nm], ``Absorption`` (calibrated molar
    absorptivity epsilon, [L mol^-1 cm^-1], i.e. the raw instrument reading
    multiplied by ``calibration_factor``).
    """
    data = pd.read_csv(
        _resolve(path),
        sep=r"\s+",
        names=["Wavelength", "Absorption"],
        skiprows=2,
    )
    data["Absorption"] = calibration_factor * data["Absorption"]
    return data


def load_xes_reference_spectra(path: str | Path = "Reference_Data.xlsx") -> pd.DataFrame:
    """
    Load the measured static XES reference spectra for the five spin states.

    Index: emission energy [eV]. Columns: ``singlet``, ``doublet``,
    ``triplet``, ``quartet``, ``quintet`` (intensity, arb. units).
    """
    data = pd.read_excel(_resolve(path))
    return data.set_index("emission energy")


def load_xes_reference_fit_params(
    path: str | Path = "Reference_spectra_rec_params.pkl",
) -> pd.DataFrame:
    """
    Load the pre-fitted 4-Voigt parameters for each spin-state reference
    spectrum (see ``xes_model.rec_xes_4voigts``).

    Index: ``I0_i, x0_i, gamma_i, sigma_i`` for ``i in 1..4``. Columns:
    ``singlet``, ``doublet``, ``triplet``, ``quartet``, ``quintet``.
    """
    return pd.read_pickle(_resolve(path))


def load_time_resolved_xes(path: str | Path = "Time_resolved_data.xlsx") -> pd.DataFrame:
    """
    Load the time-resolved (pump-probe) XES difference map.

    Index: emission energy [eV]. Columns: pump-probe delay, labeled in
    femtoseconds (the source file stores delays in picoseconds; this loader
    converts and relabels them).
    """
    data = pd.read_excel(_resolve(path), skiprows=1)
    delays_ps = data.columns[1:].to_numpy(dtype=float)
    delay_labels = [f"{round(1000 * d)} fs" for d in delays_ps]
    data.columns = ["emission_energy", *delay_labels]
    return data.set_index("emission_energy")


def save_fit_params(params: dict[str, "pd.Series | list[float]"], path: str | Path) -> None:
    """Save a ``{label: params}`` mapping (e.g. ``{"Ground State": [...], ...}``) as a pickled DataFrame."""
    pd.DataFrame(params).T.to_pickle(_resolve(path))


def load_fit_params(path: str | Path = "FitResults_GS_ES.pkl") -> pd.DataFrame:
    """Load fitted static-XAS model parameters, indexed by ``Ground State`` / ``Excited State``."""
    return pd.read_pickle(_resolve(path))
