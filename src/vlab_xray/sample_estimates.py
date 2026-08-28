"""
Physics helpers for planning a pump-probe X-ray experiment: sample
concentration, optical density, jet speed, focal spot size, and
excited-state fraction.

These are the named building blocks behind the planning arithmetic of
``notebooks/02_preparatory_estimates.ipynb`` (Ex. 10-19), which works the
same physics through ``vlab_utils``. See
``docs/physics/04_sample_and_beam_considerations.md`` for the physical
background of each formula.
"""

from __future__ import annotations

from vlab_xray.constants import AVOGADRO_NUMBER, PLANCK_CONSTANT, SPEED_OF_LIGHT


def molar_mass(atom_counts: dict[str, float], atomic_mass: dict[str, float]) -> float:
    """
    Molar mass [g/mol] of a compound given as ``{element: count}``, e.g.
    ``molar_mass({"C": 10, "H": 8, "N": 2}, ATOMIC_MASS)`` for one bipyridine
    (bpy) ligand.
    """
    return sum(count * atomic_mass[element] for element, count in atom_counts.items())


def concentration(mass_g: float, molar_mass_g_per_mol: float, volume_l: float) -> float:
    """Molar concentration [mol/L] of ``mass_g`` grams dissolved in ``volume_l`` liters."""
    return (mass_g / molar_mass_g_per_mol) / volume_l


def beer_lambert_concentration(
    optical_density: float, path_length_cm: float, molar_absorptivity: float
) -> float:
    """
    Concentration [mol/L] required to reach a given optical density
    ``OD = epsilon * c * l`` (Beer-Lambert law), for a sample of thickness
    ``path_length_cm`` and molar absorptivity ``molar_absorptivity``
    [L mol^-1 cm^-1].
    """
    return optical_density / (path_length_cm * molar_absorptivity)


def jet_flow_speed_min(jet_thickness_m: float, repetition_rate_hz: float) -> float:
    """
    Minimum liquid-jet flow speed [m/s] so that a fresh sample volume is
    delivered before the next X-ray pulse arrives, i.e. the jet must move at
    least one ``jet_thickness_m`` between pulses spaced ``1/repetition_rate_hz``
    apart.
    """
    pulse_period_s = 1 / repetition_rate_hz
    return jet_thickness_m / pulse_period_s


def photon_energy(wavelength_m: float) -> float:
    """Energy [J] of a single photon of the given ``wavelength_m``, E = h c / lambda."""
    return PLANCK_CONSTANT * SPEED_OF_LIGHT / wavelength_m


def rayleigh_focus_diameter(wavelength_m: float, focal_length_m: float, aperture_diameter_m: float) -> float:
    """
    Diffraction-limited focus *diameter* [m] -- the full width of the Airy
    disc out to its first zero, ``d = 2.44 * lambda * f / D``.

    The familiar ``1.22 * lambda * f / D`` is the Airy *radius*; the diameter
    is twice that. This matches the focal spot computed in Ex. 17 of
    ``notebooks/02_preparatory_estimates.ipynb``.
    """
    return 2.44 * wavelength_m * focal_length_m / aperture_diameter_m


def molecules_in_volume(concentration_mol_per_l: float, volume_l: float) -> float:
    """Number of molecules present in ``volume_l`` liters at the given molar concentration."""
    return AVOGADRO_NUMBER * volume_l * concentration_mol_per_l


def pulse_energy_for_molecule_count(wavelength_m: float, n_molecules: float) -> float:
    """
    Pulse energy [J] such that the number of photons in the pulse equals
    ``n_molecules`` -- the pulse energy at which, in principle, every
    molecule in the excitation volume could absorb one photon.
    """
    return photon_energy(wavelength_m) * n_molecules
