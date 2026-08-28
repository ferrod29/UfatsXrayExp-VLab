import warnings

import numpy as np
import pytest

from vlab_xray.lineshapes import cauchy, erf_step, gaussian, voigt


def test_gaussian_peak_is_unit_height_at_center():
    assert gaussian(np.array([0.0]), 0.0, 1.0)[0] == 1.0


def test_gaussian_zero_sigma_is_zero_everywhere():
    # A degenerate sigma=0 shows up in real fit results (see FitResults_GS_ES.pkl);
    # it must not raise or emit a divide-by-zero warning.
    result = gaussian(np.array([0.0, 1.0, -1.0]), 0.5, 0.0)
    assert np.all(result == 0.0)


def test_erf_step_zero_sigma_is_an_unbroadened_step_without_warnings():
    # fit_static_xas bounds the edge width from below by zero, so the optimizer
    # does probe sigma=0; that must give the exact limit rather than a nan
    # objective (which silently derails the fit).
    x = np.array([-1.0, 0.0, 1.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        y = erf_step(x, 0.0, 0.5, 0.0)
    assert np.all(np.isfinite(y))
    assert list(y) == [-0.5, 0.0, 0.5]


def test_erf_step_is_antisymmetric_around_edge():
    x0, height, sigma = 7.1, 0.5, 0.01
    x = np.array([x0 - 1, x0, x0 + 1])
    y = erf_step(x, x0, height, sigma)
    assert y[1] == 0.0
    assert y[0] == -y[2]
    assert abs(y[2] - height) < 1e-9  # far above the edge, saturates to `height`


def test_voigt_reduces_to_gaussian_when_lorentzian_width_vanishes():
    x = np.linspace(-5, 5, 201)
    amplitude, x0, sigma = 2.0, 0.3, 1.5
    v = voigt(x, amplitude, x0, gamma=1e-8, sigma=sigma)
    g = amplitude * gaussian(x, x0, sigma) / (sigma * np.sqrt(2 * np.pi))
    assert np.allclose(v, g, atol=1e-4)


def test_voigt_reduces_to_lorentzian_when_gaussian_width_vanishes():
    x = np.linspace(-20, 20, 4001)
    amplitude, x0, gamma = 3.0, -0.2, 0.7
    v = voigt(x, amplitude, x0, gamma=gamma, sigma=1e-4)
    c = cauchy(x, gamma, x0, amplitude)
    assert np.allclose(v, c, atol=1e-2)
