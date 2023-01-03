#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ i mports ------------------------ #

import pandas as pd
import numpy as np
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
planet_radius_models = []
star_radius_models = []
rho_crit = config.value[8]

for i in range(6):
    if config.value[i] == 1:
        magnetic_moment_models.append(config.setting[i])
for i in range(6, 8):
    if config.value[i] == 1:
        dynamo_density_models.append(config.setting[i])

for i in range(9,10) :
    if config.value[i] == 1 :
        planet_radius_models.append(config.setting[i])

for i in range(10,11) :
    if config.value[i] == 1 :
        star_radius_models.append(config.setting[i])

if config.value[11] == 1:
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

#data =pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/exoplanet.eu_catalog.csv', index_col=0)
data=pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/planet_test.csv',delimiter= ';', index_col=0)
# test_targets = pd.DataFrame(
#      {
#          "star_radius": [1.0, 1.0, 1.18, 1.12, 1.48, 1.48, 1.48,0.683,0.805],
#          "star_mass": [1.0, 1.0, 1.06, 1.1, 1.42, 1.42, 1.42,0.809,0.8],
#          "star_age": [4.6, 4.6, 5.5, 3.0, 1.0, 1.0, 1.0,6.5,0.6],
#          "star_distance" : [1.0, 1.0, 47.3, 1500., 15.6, 15.6, 15.6,37.89,19.3],
#          "star_magnetic_field" : [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0],
#          "radius": [1.0, 0.84, 1.42, 1.25, 1.2, 1.58, 1.48,0.389,1.138],
#          "mass": [1.0, 0.3, 0.69, 0.69, 4.4, 7.0, 10.0,0.0736,1.138],
#          "mass_sini" : [1.0, 0.3, 0.69, 1.18, 4.4, 7.0, 10.0,0.0736,1.138],
#          "orbital_period": [1.0, 0.93, 0.12, 0.34, 0.12, 0.12, 0.12,1.371,0.6224165],
#          "semi_major_axis": [5.2, 9.5, 0.045, 0.0225, 0.0489, 0.0489, 0.0489,0.053,0.031],
#          "detection_type" : [ None, None, None, None, None, None, None, None, None],
#          "eccentricity" : [0.0, 0.0, 0.0082 ,0.0, 0.0787, 0.0787, 0.0787, 0.265 ,0.0 ]
#      },
#      index=[
#          "Jupiter",
#          "Saturn",
#          "HD 209458b",
#          "OGLE-TR-56b",
#          "T-bootes (light)",
#          "T-bootes (medium)",
#          "T-bootes (heavy)",
#          "HAT-P-11",
#          "HD_189733",
#      ],
#  )


# --------------------------------------------------------- #
# ------------------------ Main --------------------------- #
# --------------------------------------------------------- #

sun = Star(name="Soleil", M=1.0, R={ "models" : star_radius_models, "radius" : 1.0}, t=AS, s=1.0, B=BSsw, L=1.0)
# sol.talk(talk=talk)
jup = Planet(
    name="Jupiter", mass=1.0, radius={ "models" : planet_radius_models, "radius" : 1.0}, distance=1.0, worb={"star_mass" : MS, "worb" : 1.0}, wrot=1.0
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

for target in data.itertuples():
    print(target.Index)
    if not np.isnan(target.semi_major_axis) and (not np.isnan(target.mass) or not np.isnan(target.mass_sini)) and not np.isnan(target.star_mass):    

        if np.isnan(target.mass) :
            planet_mass = target.mass_sini * np.sqrt(4/3.)
        else : 
            planet_mass = target.mass
        if not np.isnan(target.eccentricity) :
            planet_distance = target.semi_major_axis * (1-target.eccentricity)
        else :
            planet_distance = target.semi_major_axis

        if np.isnan(target.star_age) :
            star_age = 5.2
        else : 
            if target.star_age < 0.5 :
                star_age = 0.5
            else :
                star_age= target.star_age
            
        planet = Planet(
            name=target.Index,
            mass=planet_mass,
            radius={"models" : planet_radius_models, "radius" : target.radius},
            distance=planet_distance,
            worb={"star_mass" : target.star_mass, "worb" : target.orbital_period},
            detection_method=target.detection_type,
            wrot=target.rot_rate
        )
        star = Star(
            name="star",
            M=target.star_mass,
            R={"models" : star_radius_models, "radius" : target.star_radius},
            t=star_age,
            s=target.star_distance
        )
#B=target.star_magnetic_field,
        planet.tidal_locking(age=star.age, star_mass=star.mass, Qpp=target.tidal_Q)
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
        mytarget = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": magnetic_moment},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": magnetic_moment,
                "stellar_wind": stellar_wind,
            },
            pow_received={"star": star},
            fmax_star=1.0,
        )
        mytarget.talk(talk=talk)
        if i == 0:
            rhoc_jup = dyn_region.density
            df_target = pd.DataFrame(
                {
                    "name": mytarget.name,
                    "planet_mass": planet.mass,
                    "planet_radius": planet.radius,
                    "star_planet_distance": planet.stardist,
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
                    "magnetic_field": mytarget._mag_field,
                    "freq_max": mytarget._freq_max /1e6,
                    "pow_emission": mytarget._pow_emission / 1e14 ,
                    "flux au": mytarget.flux/ 1e-26 / 1e10,
                    "pow_received": mytarget._pow_received * 1e3 /1e-26,
                    "freq_max_star": mytarget.freq_max_star,
                },
                index=[i],
            )
        else:
            df2 = pd.DataFrame(
                {
                    "name": mytarget.name,
                    "planet_mass": planet.mass,
                    "planet_radius": planet.radius,
                    "star_planet_distance": planet.stardist,
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
                    "magnetic_field": mytarget._mag_field,
                    "freq_max": mytarget._freq_max / 1e6,
                    "pow_emission": mytarget.pow_emission / 1e14,
                    "flux au": mytarget.flux / 1e-26 / 1e10,
                    "pow_received": mytarget._pow_received * 1e3 / 1e-26,
                    "freq_max_star": mytarget.freq_max_star,
                },
                index=[i],
            )
            df_target = pd.concat([df_target, df2], ignore_index=True)
        i += 1
    #if mytarget.select_target():
            #selected_targets.append(mytarget)

df_target.to_csv("test_main.csv", sep=";", index=False)
