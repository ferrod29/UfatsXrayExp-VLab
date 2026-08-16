from vlab_xray import io as vio


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
