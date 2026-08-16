from vlab_xray.constants import ATOMIC_MASS
from vlab_xray.sample_estimates import (
    beer_lambert_concentration,
    jet_flow_speed_min,
    molar_mass,
    photon_energy,
    rayleigh_focus_diameter,
)


def test_molar_mass_of_febpy3_matches_hand_calculation():
    # BPY = C10H8N2; hand-calculated in the original vlab_sample_estimates.ipynb
    m_bpy = molar_mass({"C": 10, "H": 8, "N": 2}, ATOMIC_MASS)
    m_febpy3 = ATOMIC_MASS["Fe"] + 3 * m_bpy
    assert abs(m_bpy - 156.188) < 1e-3
    assert abs(m_febpy3 - 524.404) < 1e-3


def test_beer_lambert_concentration_inverts_optical_density():
    path_length_cm, epsilon = 25e-4, 5000.0
    c = beer_lambert_concentration(1.0, path_length_cm, epsilon)
    assert abs(epsilon * c * path_length_cm - 1.0) < 1e-9


def test_jet_flow_speed_matches_thickness_over_period():
    thickness_m, rep_rate_hz = 10e-6, 1.1e6
    v_min = jet_flow_speed_min(thickness_m, rep_rate_hz)
    assert abs(v_min - thickness_m * rep_rate_hz) < 1e-12


def test_photon_energy_400nm_is_about_3ev():
    e_joule = photon_energy(400e-9)
    e_ev = e_joule / 1.602176634e-19
    assert 3.0 < e_ev < 3.2


def test_rayleigh_focus_diameter_scales_linearly_with_wavelength():
    d1 = rayleigh_focus_diameter(400e-9, 0.25, 0.01)
    d2 = rayleigh_focus_diameter(800e-9, 0.25, 0.01)
    assert abs(d2 - 2 * d1) < 1e-12
