#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ i mports ------------------------ #

import pandas as pd
from star import Star
from planet import Planet
from dynamo_region import DynamoRegion
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind
from target import Target

# --------------------------------------------------------- #
# ------------------- Physical constants ------------------ #

MS = 1.989e30  # kg
RS = 6.96342e8  # m
AS = 4.6e9  # yo
BSsw = 1  # T
LS = 3.826e26  # W


MJ = 1.8986e27  # kg
RJ = 69911e3  # m
wJ = 1.77e-4  # s-1

ME = 5.97237e24  # kg
RE = 6371.0e3  # m
wE = 7.27e-5  # s-1

rc = 0.85 * RJ
rhoc = 1800  # kg/m3


# --------------------------------------------------------- #
# -------------- Configuration settings input ------------- #

config = pd.read_csv(
    r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/parametres.csv",
    delimiter=";",
)

# Setting the models to use

dynamo_density_models = []
magnetic_moment_models = []
rho_crit = config.value[8]

for i in range(6):
    if config.value[i] == 1:
        magnetic_moment_models.append(config.setting[i])
for i in range(6, 8):
    if config.value[i] == 1:
        dynamo_density_models.append(config.setting[i])

if config.value[9] == 1:
    talk = True
else:
    talk = False

# Criterions for target selection

selection_config = pd.read_csv(
    r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/selection.csv",
    delimiter=";",
)


# --------------------------------------------------------- #
# ---------------------- Data input ----------------------- #

# data =pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/exoplanet.eu_catalog.csv')
test_targets = pd.DataFrame(
    {
        "star_radius": [1.0, 1.0, 1.18, 1.12, 1.48, 1.48, 1.48],
        "star_mass": [1.0, 1.0, 1.06, 1.1, 1.42, 1.42, 1.42],
        "star_age": [4.6, 4.6, 5.5, 3.0, 1.0, 1.0, 1.0],
        "obs_dist" : [1.0, 1.0, 47.3, 1500., 15.6, 15.6, 15.6],
        "planet_radius": [1.0, 0.84, 1.42, 1.25, 1.2, 1.58, 1.48],
        "planet_mass": [1.0, 0.3, 0.69, 1.18, 4.4, 7.0, 10.0],
        "planet_orbital_frequency": [1.0, 0.93, 0.12, 0.34, 0.12, 0.12, 0.12],
        "planet_distance": [5.2, 9.5, 0.045, 0.0225, 0.0489, 0.0489, 0.0489],
    },
    index=[
        "Jupiter",
        "Saturn",
        "HD 209458b",
        "OGLE-TR-56b",
        "T-bootes (light)",
        "T-bootes (medium)",
        "T-bootes (heavy)",
    ],
)


# --------------------------------------------------------- #
# ------------------------ Main --------------------------- #
# --------------------------------------------------------- #

sun = Star(name="Soleil", M=1.0, R=1.0, t=AS, s=1.0, B=BSsw, L=1.0)
# sol.talk(talk=talk)
jup = Planet(
    name="Jupiter", mass=1.0, radius=1.0, distance=1.0, axis=1.0, worb=1.0, wrot=1.0
)
jup.talk(talk=talk)
jup.tidal_locking(age=4.6e9, star_mass=1.0)
dyn_region_jup = DynamoRegion.from_planet(planet=jup, rhocrit=rho_crit)
mag_moment_jup = MagneticMoment(models=magnetic_moment_models, Mm=1.56e27, Rs=1.0)
# mag_moment_jup.magnetic_moment(
# dynamo=dyn_region_jup, planet=jup, star=sun, normalize=False
# )
mag_moment_jup.talk(talk=True)
selected_targets = []

i = 0

for target in test_targets.itertuples():
    planet = Planet(
        name=target.Index,
        mass=target.planet_mass,
        radius=target.planet_radius,
        distance=target.planet_distance,
        axis=1.0,
        worb=target.planet_orbital_frequency,
    )
    star = Star(
        name="star",
        M=target.star_mass,
        R=target.star_radius,
        t=target.star_age,
        s=target.obs_dist,
        B=1.0,
        L=1.0,
    )
    planet.tidal_locking(age=star.age, star_mass=star.mass)
    planet.talk(talk=talk)
    star.talk(talk=talk)
    dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=rho_crit)
    dyn_region.normalize(other=dyn_region_jup)
    dyn_region.talk(talk=talk)
    stellar_wind = StellarWind.from_system(star=star, planet=planet)
    stellar_wind.talk(talk=talk)
    magnetic_moment = MagneticMoment(models=magnetic_moment_models, Mm=1.0, Rs=1.0)
    magnetic_moment.magnetic_moment(dynamo=dyn_region, planet=planet, star=star)
    magnetic_moment.magnetosphere_radius(mag_moment_jup, stellar_wind=stellar_wind)
    magnetic_moment.talk(talk=talk)
    target = Target(
        name=planet.name,
        mag_field={"planet": planet, "magnetic_moment": magnetic_moment},
        flux={
            "planet": planet,
            "star": star,
            "magnetic_moment": magnetic_moment,
            "stellar_wind": stellar_wind,
        },
        pow_received={
            "planet": planet,
            "magnetic_moment": magnetic_moment,
            "stellar_wind": stellar_wind,
        },
        Pem=1.0,
        fmax_star=1.0,
    )
    target.talk(talk=talk)
    if i == 0:
        rhoc_jup = dyn_region.density
        df_target = pd.DataFrame(
            {
                "name": target.name,
                "planet_mass": planet.mass,
                "planet_radius": planet.radius,
                "star_planet_distance": planet.stardist,
                "semi_major_axis": planet.sm_axis,
                "planet_rotation_rate": planet.rotrate,
                "planet_orbital_period": planet.orbitperiod,
                "star_mass": star.mass,
                "star_radius": star.radius,
                "star_age": star.age,
                "earth_distance": star.obs_dist,
                "luminosity": star.luminosity,
                "dynamo_density": dyn_region.density / rhoc_jup,
                "dynamo_radius": dyn_region.radius / planet.radius,
                "magnetic_moment": magnetic_moment.mag_moment,
                "standoff_distance": magnetic_moment.normalize_standoff_dist(
                    planet=planet
                ),
                "sw_density": stellar_wind.density,
                "sw_velocity": stellar_wind.effective_velocity,
                "coronal_temperature": stellar_wind.corona_temperature,
                "sw_magfield": stellar_wind.mag_field,
                "magnetic_field": target._mag_field,
                "freq_max": target._freq_max,
                "flux": target._flux,
                "pow_received": target._pow_received,
                "pow_emission": target.pow_emission,
                "freq_max_star": target.freq_max_star,
            },
            index=[i],
        )
    else:
        df2 = pd.DataFrame(
            {
                "name": target.name,
                "planet_mass": planet.mass,
                "planet_radius": planet.radius,
                "star_planet_distance": planet.stardist,
                "semi_major_axis": planet.sm_axis,
                "planet_rotation_rate": planet.rotrate,
                "planet_orbital_period": planet.orbitperiod,
                "star_mass": star.mass,
                "star_radius": star.radius,
                "star_age": star.age,
                "earth_distance": star.obs_dist,
                "luminosity": star.luminosity,
                "dynamo_density": dyn_region.density / rhoc_jup,
                "dynamo_radius": dyn_region.radius / planet.radius,
                "magnetic_moment": magnetic_moment.mag_moment,
                "standoff_distance": magnetic_moment.normalize_standoff_dist(
                    planet=planet
                ),
                "sw_density": stellar_wind.density,
                "sw_velocity": stellar_wind.effective_velocity,
                "coronal_temperature": stellar_wind.corona_temperature,
                "sw_magfield": stellar_wind.mag_field,
                "magnetic_field": target._mag_field,
                "freq_max": target._freq_max / 1e9,
                "flux": target._flux / 1e-26,
                "pow_received": target._pow_received,
                "pow_emission": target.pow_emission,
                "freq_max_star": target.freq_max_star,
            },
            index=[i],
        )
        df_target = pd.concat([df_target, df2], ignore_index=True)
    i += 1
    if target.select_target():
        selected_targets.append(target)

df_target.to_csv("test_main.csv", sep=";", index=False)
