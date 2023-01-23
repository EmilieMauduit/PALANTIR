# -*- coding: utf-8 -*-


from dynamo_region import DynamoRegion
from planet import Planet
from star import Star
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind
from target import Target
import pytest


def test_target():
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

    assert isinstance(m_mag, MagneticMoment)
    m_mag.magnetic_moment(dynamo=dynamo, planet=planet, star=star)

    s_wind = StellarWind.from_system(star=star, planet=planet)
    assert isinstance(s_wind, StellarWind)

    my_target = Target(
        name=planet.name,
        mag_field={"planet": planet, "magnetic_moment": m_mag},
        pow_emission={
            "planet": planet,
            "star": star,
            "magnetic_moment": m_mag,
            "stellar_wind": s_wind,
        },
        pow_received={"star": star},
        fmax_star=1.0,
    )

    assert isinstance(my_target, Target)

    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet},
            pow_emission={
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
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"magnetic_moment": m_mag},
            pow_emission={
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
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet,
                "star": star,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "star": star,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "magnetic_moment": m_mag,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "planet": planet,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
    with pytest.raises(KeyError):
        my_target = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": m_mag},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            pow_received={
                "magnetic_moment": m_mag,
                "stellar_wind": s_wind,
            },
            fmax_star=1.0,
        )
