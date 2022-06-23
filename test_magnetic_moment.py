# -*- coding: utf-8 -*-
from dynamo_region import DynamoRegion
from planet import Planet
from star import Star
from magnetic_moment import MagneticMoment
import pytest


def test_magnetic_moment():
    planet = Planet("Jupiter", 1.8986e27, 71492e3, 5.2 * 149597870700, 1.0, 1.77e-4)
    star = Star("Sun", 1.989e30, 6.96342e8, 4.6e9, 1.0, 1.0, 3.826e26)
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
