import pytest
import numpy as np

from planet import Planet


@pytest.fixture
def planet_data():
    """Fixture for creating a sample planet with default values."""
    radius_data = {"models": ["radius_original"], "radius": np.nan}
    worb_data = {"star_mass": 1.0, "worb": np.nan}
    luminosity_data = {"models": ["Burrows"], "luminosity": np.nan, "star_age": 4.5}
    return {
        "name": "TestPlanet",
        "mass": 1.0,
        "radius": radius_data,
        "distance": 5.0,
        "worb": worb_data,
        "luminosity": luminosity_data,
        "detection_method": "transit",
        "wrot": 1.0,
    }


@pytest.fixture
def planet_instance(planet_data):
    """Create a default Planet instance for testing."""
    return Planet(**planet_data)


# Test for the Constructor and Property Getters
def test_constructor(planet_instance):
    assert planet_instance.name == "TestPlanet"
    assert planet_instance.mass == 1.0
    assert planet_instance.stardist == 5.0
    assert planet_instance.detection_method == "transit"
    assert isinstance(planet_instance.radius, float)
    assert isinstance(planet_instance.orbitperiod, float)
    assert isinstance(planet_instance.luminosity, float)


# Test for Unnormalized Mass Calculation
def test_unnormalize_mass(planet_instance):
    expected_mass = 1.0 * 1.8986e27  # mass in kg
    assert planet_instance.unnormalize_mass() == pytest.approx(expected_mass, rel=1e-5)


# Test for Unnormalized Radius Calculation
def test_unnormalize_radius(planet_instance):
    jupiter_radius = 71492e3  # in meters
    expected_radius = planet_instance.radius * jupiter_radius
    assert planet_instance.unnormalize_radius() == pytest.approx(expected_radius, rel=1e-5)


# Test for Radius Expansion Factor Calculation
def test_radius_expansion(planet_instance):
    Teq = planet_instance.radius_expansion(luminosity=1.0, eccentricity=0.1)
    assert isinstance(Teq, float)
    assert Teq > 0
    assert planet_instance.radius > 0


# Test for Radius Expansion (Old Method)
def test_radius_expansion_old(planet_instance):
    Teq = planet_instance.radius_expansion_old(luminosity=1.0, eccentricity=0.1)
    assert isinstance(Teq, float)
    assert Teq > 0
    assert planet_instance.radius > 0


# Test for Orbital Period Calculation
def test_unnormalize_orbital_period(planet_instance):
    jupiter_orbital_period = 1.77e-4  # in s^-1
    expected_orbit = planet_instance.orbitperiod * jupiter_orbital_period
    assert planet_instance.unnormalize_orbital_period() == pytest.approx(expected_orbit, rel=1e-5)


# Test for Rotation Rate Calculation
def test_unnormalize_rotrate(planet_instance):
    jupiter_rotation_rate = 1.77e-4  # in s^-1
    expected_rotrate = planet_instance.rotrate * jupiter_rotation_rate
    assert planet_instance.unnormalize_rotrate() == pytest.approx(expected_rotrate, rel=1e-5)


# Test for Tidal Locking
def test_tidal_locking(planet_instance):
    star_mass = 1.0
    age = 4.5
    initial_rotation_rate = planet_instance.rotrate
    planet_instance.tidal_locking(age=age, star_mass=star_mass)
    if planet_instance.stardist <= planet_instance._calculate_synchro_dist(
        planet_instance.rotrate, Qp=3.16e5, tsync=age, planet_mass=planet_instance.mass,
        planet_radius=planet_instance.radius, star_mass=star_mass):
        assert planet_instance.rotrate == planet_instance.orbitperiod
    else:
        assert planet_instance.rotrate == initial_rotation_rate


# Test _calculate_radius Static Method
def test_calculate_radius():
    radius = Planet._calculate_radius(models=["radius_original"], mass=1.0)
    assert isinstance(radius, float)
    assert radius > 0


# Test _calculate_luminosity Static Method
def test_calculate_luminosity():
    luminosity = Planet._calculate_luminosity(models=["Burrows"], planet_mass=1.0, star_age=4.5)
    assert isinstance(luminosity, float)
    assert luminosity > 0


# Test _calculate_synchro_dist Static Method
def test_calculate_synchro_dist():
    distance = Planet._calculate_synchro_dist(
        wrot=1.0, Qp=3.16e5, tsync=4.5, planet_mass=1.0, planet_radius=1.0, star_mass=1.0
    )
    assert isinstance(distance, float)
    assert distance > 0
