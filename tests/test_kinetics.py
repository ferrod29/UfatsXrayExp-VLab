import numpy as np

from vlab_xray.kinetics import (
    build_transient_xas_map,
    build_transient_xes_map_with_quintet_decay,
    excited_state_population,
    product_state_population,
)


def test_population_is_negligible_before_the_pulse_arrives():
    delays = np.linspace(-500, 500, 501)
    excited = excited_state_population(delays, t0=0, i0=1, sigma=5, tau=20)
    assert excited[0] < 1e-6  # far (100 sigma) before t0


def test_excited_state_decays_while_product_state_saturates():
    # Long after the pump pulse and several lifetimes later, essentially all
    # population should have moved from the "excited" to the "product" state.
    delays = np.linspace(-200, 400, 601)
    i0, tau = 2.0, 10
    excited = excited_state_population(delays, t0=0, i0=i0, sigma=2, tau=tau)
    product = product_state_population(delays, t0=0, i0=i0, sigma=2, tau=tau)
    assert excited[-1] < 1e-3
    assert abs(product[-1] - i0) < 1e-3


def test_transient_xas_map_vanishes_if_ground_equals_excited_spectrum():
    energy = np.linspace(7.1, 7.2, 50)
    spectrum = np.sin(energy * 100)  # arbitrary but identical for GS/ES
    delays = np.linspace(-100, 500, 20)
    transient = build_transient_xas_map(spectrum, spectrum, delays, t0=0, i0=1, sigma=20, tau=100)
    assert np.allclose(transient, 0.0)


def test_quintet_decay_extension_matches_exponential_after_cutoff():
    energy = np.linspace(7030, 7080, 30)
    singlet = np.ones_like(energy)
    triplet = 2 * np.ones_like(energy)
    quintet = 3 * np.ones_like(energy)
    t0, cascade_window, lifetime = 0.0, 1000.0, 50_000.0

    delays = np.array([cascade_window + 10_000, cascade_window + 20_000])
    transient = build_transient_xes_map_with_quintet_decay(
        singlet, triplet, quintet, delays,
        t0=t0, i0=1, sigma=100, tau=100, quintet_lifetime=lifetime,
        cascade_window=cascade_window,
    )
    ratio = transient[:, 1] / transient[:, 0]
    expected_ratio = np.exp(-((delays[1] - t0) - (delays[0] - t0)) / lifetime)
    assert np.allclose(ratio, expected_ratio, rtol=1e-6)
