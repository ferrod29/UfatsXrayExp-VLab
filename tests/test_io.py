import numpy as np

from vlab_xray import io as vio
from vlab_xray.xas_model import xas_model


def test_load_pump_probe_xas_has_expected_columns():
    data = vio.load_pump_probe_xas()
    assert list(data.columns) == ["Energy", "ON", "OFF", "Diff", "ES"]
    assert len(data) > 0
    assert 7.0 < data["Energy"].min() < data["Energy"].max() < 8.0


def test_load_cxro_transmission_water_and_iron():
    water = vio.load_cxro_transmission("cxro_water_25u.txt")
    iron = vio.load_cxro_transmission("cxro_iron_00595u.txt")
    for df in (water, iron):
        assert list(df.columns) == ["Energy", "Transmission"]
        assert (df["Transmission"] >= 0).all()
        assert (df["Transmission"] <= 1).all()


def test_load_uvvis_spectrum_is_calibrated():
    data = vio.load_uvvis_spectrum()
    assert list(data.columns) == ["Wavelength", "Absorption"]
    assert data["Wavelength"].min() >= 0


def test_load_xes_reference_spectra_has_five_spin_states():
    data = vio.load_xes_reference_spectra()
    assert list(data.columns) == ["singlet", "doublet", "triplet", "quartet", "quintet"]
    assert data.index.name == "emission energy"


def test_load_xes_reference_fit_params_shape():
    params = vio.load_xes_reference_fit_params()
    assert params.shape == (16, 5)


def test_load_time_resolved_xes_columns_are_fs_labels():
    data = vio.load_time_resolved_xes()
    assert all(col.endswith(" fs") for col in data.columns)


def test_load_fit_params_has_ground_and_excited_state():
    params = vio.load_fit_params()
    assert "Ground State" in params.index
    assert "Excited State" in params.index


def test_stored_fit_params_actually_reproduce_the_measured_spectra():
    """The stored fit must match the *current* ``xas_model`` edge convention.

    ``FitResults_GS_ES.pkl`` is a derived artefact: parameters fit to the
    ``OFF``/``ES`` columns. Nothing else notices if it drifts out of step with
    the model that consumes it -- an earlier revision was fit against a
    logistic-sigmoid edge and, read back through the erf edge, produced a
    model off by a third of the edge jump while still loading cleanly.
    """
    data = vio.load_pump_probe_xas()
    x = data["Energy"].to_numpy()
    params = vio.load_fit_params()
    window = (x > 7.1) & (x < 7.62)
    for row, column in (("Ground State", "OFF"), ("Excited State", "ES")):
        y = data[column].to_numpy()
        residual = xas_model(x[window], params.loc[row].dropna().values) - y[window]
        rms = np.sqrt(np.mean(residual**2))
        assert rms < 0.01, f"{row}: RMS residual {rms:.4f} (edge jump is ~0.2)"


def test_save_fit_params_round_trips_parameter_vectors_of_different_lengths(tmp_path):
    # A ground-state fit uses one resonance Gaussian fewer than the excited
    # state, so the saved rows are ragged and the short one is NaN-padded.
    path = tmp_path / "fit.pkl"
    vio.save_fit_params({"Short": [1.0, 2.0, 3.0], "Long": [4.0, 5.0, 6.0, 7.0]}, path)
    frame = vio.load_fit_params(path)
    assert frame.shape == (2, 4)
    assert list(frame.loc["Short"].dropna()) == [1.0, 2.0, 3.0]
    assert list(frame.loc["Long"].dropna()) == [4.0, 5.0, 6.0, 7.0]
