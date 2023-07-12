# -*- coding: utf-8 -*-
from dynamo_region import DynamoRegion
from planet import Planet
from star import Star
from magnetic_moment import MagneticMoment
import pytest


def test_magnetic_moment():
    planet = Planet(
        name="Jupiter",
        mass=1.0,
        radius={"models": ["original_radius"], "radius": 1.0},
        distance=1.0,
        worb={"star_mass": 1.989e30, "worb": 1.0},
        wrot=1.0,
    )
    star = Star(
        name="Soleil",
        mass=1.0,
        radius={"models": ["Tout"], "radius": 1.0},
        age=4.6,
        obs_dist=1.0,
        magfield=1.0,
    )
    dynamo = DynamoRegion.from_planet(planet, rhocrit=700)

    assert isinstance(planet, Planet)
    assert isinstance(star, Star)
    assert isinstance(dynamo, DynamoRegion)

    planet.tidal_locking(age=star.age, star_mass=star.age)

    m_mag = MagneticMoment(["busse", "sano"], 1.0, 1.0)

    m_mag.magnetic_moment(dynamo=dynamo, planet=planet, star=star)

    with pytest.raises(ValueError):
        m_mag.magnetic_moment(
            dynamo=dynamo, planet=planet, star=star, Mmean=True, Mmax=True
        )

    with pytest.warns()
