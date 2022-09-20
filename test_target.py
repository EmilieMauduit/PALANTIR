# -*- coding: utf-8 -*-


from dynamo_region import DynamoRegion
from planet import Planet
from star import Star
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind
from target import Target
import pytest


def test_target():
    planet = Planet("Jupiter", 1.8986e27, 71492e3, 5.2 * 149597870700, 1.0, 1.77e-4)
    star = Star("Sun", 1.989e30, 6.96342e8, 4.6, 1.0, 1.0, 3.826e26)
    dynamo = DynamoRegion.from_planet(planet, rhocrit=700)

    assert isinstance(planet, Planet)
    assert isinstance(star, Star)
    assert isinstance(dynamo, DynamoRegion)

    planet.tidal_locking(age=star.age, star_mass=star.age)

    m_mag = MagneticMoment(["busse", "sano"], 1.0, 1.0)

    assert isinstance(m_mag, MagneticMoment)
    m_mag.magnetic_moment(dynamo=dynamo, planet=planet, star=star)

    s_wind = StellarWind.from_system(star=star, planet=planet)
    assert isinstance(s_wind, StellarWind)

    my_target = Target(
        name=planet.name,
        mag_field={"planet": planet, "magnetic_moment": m_mag},
        flux={
            "planet": planet,
            "star": star,
            "magnetic_moment": m_mag,
            "stellar_wind": s_wind,
        },
        pow_received={
            "planet": planet,
            "magnetic_moment": m_mag,
            "stellar_wind": s_wind,
        },
        Pem=1.0,
        fmax_star=1.0,
    )
    assert isinstance(my_target, Target)

    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "star": star,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            flux={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            Pem=1.0,
            fmax_star=1.0,
        )
