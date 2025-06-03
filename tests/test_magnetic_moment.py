# -*- coding: utf-8 -*-
from palantir.prediction_tools import (
    DynamoRegion,
    Planet,
    Star,
    MagneticMoment,
    dynamo_region
)
import numpy as np
import pytest

# --------------------------------------------------------- #
# ------------------- Physical constants ------------------ #

MS = 1.989e30  # kg
RS = 6.96342e8  # m
AS = 4.6  # yo
BSsw = 1  # T
LS = 3.826e26  # W


MJ = 1.8986e27  # kg
RJ = 69911e3  # m
wJ = 1.77e-4  # s-1

@pytest.fixture
def planet_data():
    """Fixture for creating a sample planet with default values."""
    radius_data = {"models": ["radius_original"], "radius": np.nan}
    worb_data = {"star_mass": 1.0, "worb": np.nan}
    luminosity_data = {"models": ["Baraffe_noirrad"], "luminosity": np.nan, "star_age": 4.5}
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


# --------------------------------------------------------- #
# ------------------- TESTS ------------------ #

def test_magnetic_moment():
    Jupiter = Planet(
        name="Jupiter",
        mass=1.0,
        radius={"models": ["radius_original"], "radius": 1.0},
        distance=5.2,
        worb={"star_mass": MS, "worb": 1.0},
        luminosity={
            "models": ["Baraffe_noirrad"],
            "luminosity": np.nan,
            "star_age": 4.6,
        },
        wrot=1.0,
    )
    planet = Planet(
        name="TestPlanet",
        mass=1.0,
        radius={"models": ["radius_original"], "radius": 1.0},
        distance=5.2,
        worb={"star_mass": MS, "worb": 1.0},
        luminosity={
            "models": ["Baraffe_noirrad"],
            "luminosity": np.nan,
            "star_age": 4.6,
        },
        wrot=1.0,
    )
    star = Star(
        name="Sun",
        mass=1.0,
        radius={"models": ['Tout'], "radius": 1.0},
        age=AS,
        obs_dist=1.0,
        sp_type ='GV',
    )
    dynamo = DynamoRegion.from_planet(planet, rhocrit=700)

    assert isinstance(planet, Planet)
    assert isinstance(star, Star)
    assert isinstance(dynamo, DynamoRegion)

    planet.tidal_locking(age=star.age, star_mass=star.mass)

    m_mag = MagneticMoment(["busse","sano"], Mm=1.56e27, Rs=1.0)

    m_mag.magnetic_moment(dynamo=dynamo, planet=planet, jup=planet, dynamo_jup=dynamo_region)

    with pytest.raises(ValueError):
        m_mag.magnetic_moment(
            dynamo=dynamo, planet=planet, jup=planet, dynamo_jup=dynamo_region, Mmean=True, Mmax=True
        )