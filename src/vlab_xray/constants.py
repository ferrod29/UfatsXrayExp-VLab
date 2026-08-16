"""
Physical constants and atomic masses used across the sample-estimates and
spectroscopy calculations.

Values are standard atomic weights (IUPAC) and CODATA physical constants,
given to the precision actually needed for the Praktikum calculations. See
``docs/physics/04_sample_and_beam_considerations.md`` for how these feed
into the sample-estimates exercise.
"""

#: Planck constant [J s]
PLANCK_CONSTANT = 6.62607015e-34

#: Speed of light in vacuum [m/s]
SPEED_OF_LIGHT = 2.99792458e8

#: Avogadro constant [1/mol]
AVOGADRO_NUMBER = 6.02214076e23

#: Standard atomic weights [g/mol], for the elements appearing in
#: [Fe(bpy)3]Cl2 . 3H2O and [Fe(terpy)2]2+.
ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "Fe": 55.84,
    "Cl": 35.45,
}

#: Molar mass of water [g/mol], from ATOMIC_MASS.
MOLAR_MASS_WATER = 2 * ATOMIC_MASS["H"] + ATOMIC_MASS["O"]

#: Density of metallic iron [g/cm^3]
DENSITY_IRON = 7.874
