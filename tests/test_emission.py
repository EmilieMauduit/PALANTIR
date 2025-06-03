# -*- coding: utf-8 -*-

from palantir.prediction_tools import (
    DynamoRegion,
    Planet,
    Star,
    MagneticMoment,
    StellarWind,
    Emission,
    dynamo_region
)

import pytest
import numpy as np

@pytest.fixture
def star_data():
    """Fixture for creating a sample star with default values."""
    radius_data = {"models": ["Tout"], "radius": np.nan}
    return {
        "name": "TestStar",
        "mass": 1.0,
        "radius": radius_data,
        "age": 4.6,
        "obs_dist" : 1.0,
        "sp_type": "GV",
    }


@pytest.fixture
def star_instance(star_data):
    """Create a default Star instance for testing."""
    return Star(**star_data)

@pytest.fixture
def planet_data():
    """Fixture for creating a sample planet with default values."""
    radius_data = {"models": ["radius_original"], "radius": np.nan}
    worb_data = {"star_mass": 1.0, "worb": np.nan}
    luminosity_data = {"models": ["Burrows","Baraffe_irrad","Baraffe_noirrad"], "luminosity": np.nan, "star_age": 4.5}
    return {
        "name": "TestPlanet",
        "mass": 1.0,
        "radius": radius_data,
        "distance": 5.0,
        "eccentricity" : 0,
        "worb": worb_data,
        "luminosity": luminosity_data,
        "detection_method": "transit",
        "wrot": 1.0,
    }


@pytest.fixture
def planet_instance(planet_data):
    """Create a default Planet instance for testing."""
    return Planet(**planet_data)

def test_emission(planet_instance, star_instance):
    

    dynamo = DynamoRegion.from_planet(planet_instance, rhocrit=700)

    assert isinstance(planet_instance, Planet)
    assert isinstance(star_instance, Star)
    assert isinstance(dynamo, DynamoRegion)

    planet_instance.tidal_locking(age=star_instance.age, star_mass=star_instance.age)

    m_mag = MagneticMoment(["busse", "sano"], 1.0, 1.0)

    assert isinstance(m_mag, MagneticMoment)
    m_mag.magnetic_moment(dynamo=dynamo, planet=planet_instance, jup=planet_instance, dynamo_jup=dynamo_region)

    s_wind = StellarWind.from_system(star=star_instance, planet=planet_instance)
    assert isinstance(s_wind, StellarWind)

    my_emission = Emission(
        name=planet_instance.name,
        mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
        pow_emission={
            "planet": planet_instance,
            "star": star_instance,
            "magnetic_moment": m_mag,
            "stellar_wind": s_wind,
        },
        pow_received={"star": star_instance},
        fmax_star=1.0,
    )

    assert isinstance(my_emission, Emission)

    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "star": star_instance,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "magnetic_moment": m_mag,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet_instance,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_emission = Emission(
            name=planet_instance.name,
            mag_field={"planet": planet_instance, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet_instance,
                "star": star_instance,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
