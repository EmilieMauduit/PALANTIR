#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ imports ------------------------ #

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

# sol = Star(name="Soleil",M=MS,R=RS,t=AS,s=1.0,B=BSsw,L=LS)
# sol.talk(talk=talk)
# jup = Planet(name="Jupiter", mass=MJ,radius=RJ,distance=1.0,axis=1.0,worb=wJ,wrot=wJ)
# jup.talk(talk=talk)

selected_targets = []

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
        s=1.0,
        B=1.0,
        L=1.0,
    )
    planet.tidal_locking(age=star.age, star_mass=star.mass)
    planet.talk(talk=talk)
    star.talk(talk=talk)
    dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=rho_crit)
    dyn_region.talk(talk=talk)
    stellar_wind = StellarWind.from_system(star=star, planet=planet)
    stellar_wind.talk(talk=talk)
    magnetic_moment = MagneticMoment(models=magnetic_moment_models, Mm=1.0, Rs=1.0)
    magnetic_moment.magnetic_moment(dynamo=dyn_region, planet=planet, star=star)
    magnetic_moment.magnetosphere_radius(stellar_wind=stellar_wind)
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
    if target.select_target():
        selected_targets.append(target)
