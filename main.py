#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

import pandas as pd
import numpy as np
from star import Star
from planet import Planet
from dynamo_region import DynamoRegion
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind,parker
from target import Target

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

ME = 5.97237e24  # kg
RE = 6371.0e3  # m
wE = 7.27e-5  # s-1

#rc = 0.85 * RJ
#rhoc = 1800  # kg/m3


# --------------------------------------------------------- #
# -------------- Configuration settings input ------------- #

config = pd.read_csv(
    r"/Users/emauduit/Documents/These/Sélection des cibles/Programmes/parametres.csv",
    delimiter=";",
)

# Setting the models to use

dynamo_density_models = []
magnetic_moment_models = []
planet_radius_models = []
star_radius_models = []
rho_crit = config.value[8]

magnetic_moment_models = [config.setting[i] for i in range(6) if config.value[i] == 1]
dynamo_density_models = [config.setting[i] for i in range(6,8) if config.value[i] == 1]
planet_radius_models = [config.setting[i] for i in range(9,10) if config.value[i] == 1]
star_radius_models = [config.setting[i] for i in range(10,11) if config.value[i] == 1]
planet_luminosity_models = [config.setting[i] for i in range(11,14) if config.value[i] == 1]

if config.value[14] == 1:
    talk = True
else:
    talk = False

# Criterions for target selection

selection_config = pd.read_csv(
    r"/Users/emauduit/Documents/These/Sélection des cibles/Programmes/selection.csv",
    delimiter=";",
)


# --------------------------------------------------------- #
# ---------------------- Data input ----------------------- #

#data =pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/planet_test_bis.csv', index_col=0)
#data =pd.read_csv(r'/Users/emauduit/Documents/These/Sélection des cibles/Proceeding PRE/M_VS_Mmag.csv', index_col=0)
data = pd.read_csv(
    r"/Users/emauduit/Documents/These/Sélection des cibles/Programmes/planet_test_unique.csv",
    delimiter=";",
    index_col=0,
    )
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

sun = Star(
    name="Soleil",
    mass=1.0,
    radius={"models": star_radius_models, "radius": 1.0},
    age=AS,
    obs_dist=1.0
)
#sun.talk(talk=talk)
jup = Planet(
    name="Jupiter",
    mass=1.0,
    radius={"models": planet_radius_models, "radius": 1.0},
    distance=5.2,
    worb={"star_mass": MS, "worb": 1.0},
    luminosity = {"models" : planet_luminosity_models, "luminosity" : np.nan, 'star_age' : 4.6},
    wrot=1.0,
)
jup.talk(talk=talk)
jup.tidal_locking(age=4.6e9, star_mass=1.0)
dyn_region_jup = DynamoRegion.from_planet(planet=jup, rhocrit=rho_crit)
mag_moment_jup = MagneticMoment(models=magnetic_moment_models, Mm=1.56e27, Rs=1.0)
#mag_moment_jup.talk(talk=True)
vjup,vejup,nejup,Tjup=parker(star = sun, planet= jup, T=0.81e6)
sw_jup = StellarWind(ne=nejup,ve=vejup,Tcor=Tjup,Bsw = {"planet" : jup, "star" : sun, "vsw" : vjup})
sw_jup.talk(talk=talk)
selected_targets = []

i = 0

for target in data.itertuples():
    print(target.Index)
    #if 'PSR' in target.Index :
    #    print('Warning : {} has been skipped'.format(target.Index))
    #    continue
    
    if (
        not np.isnan(target.semi_major_axis)
        and (not np.isnan(target.mass) or not np.isnan(target.mass_sini))
        and not np.isnan(target.star_mass)
    ):

        if np.isnan(target.mass):
            planet_mass = target.mass_sini * np.sqrt(4 / 3.0)
        else:
            planet_mass = target.mass
        if not np.isnan(target.eccentricity):
            planet_distance = target.semi_major_axis * (1 - target.eccentricity)
        else:
            planet_distance = target.semi_major_axis

        if np.isnan(target.star_age):
            star_age = 5.2
        else:
            if target.star_age <= 1e-7 :
                raise ValueError ("Stellar age is too small (< 100 yr) !!!")
            elif target.star_age < 0.5:
                star_age = 0.5
            else:
                star_age = target.star_age

        planet = Planet(
            name=target.Index,
            mass=planet_mass,
            radius={"models": planet_radius_models, "radius": target.radius},
            distance=planet_distance,
            worb={"star_mass": target.star_mass, "worb": np.nan},
            luminosity = {"models" : planet_luminosity_models, "luminosity" : np.nan, 'star_age' : star_age},
            detection_method=target.detection_type,
            wrot=1.0,
        )
        star = Star(
            name="star",
            mass=target.star_mass,
            radius={"models": star_radius_models, "radius": target.star_radius},
            age=star_age,
            obs_dist=target.star_distance,
        )
    
        if np.isnan(target.radius):
            planet.radius_expansion(
                luminosity=star.luminosity, eccentricity=target.eccentricity
            )
        try :
            planet.tidal_locking(age=star.age, star_mass=star.mass)
        except OverflowError :
            continue
        planet.talk(talk=True)
        star.talk(talk=talk)
        dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=rho_crit)
        dyn_region.normalize(other=dyn_region_jup)
        dyn_region.talk(talk=talk)
        try :
            stellar_wind = StellarWind.from_system(star=star, planet=planet)
            stellar_wind.talk(talk=talk)
        except ValueError :
            continue
        magnetic_moment = MagneticMoment(models=magnetic_moment_models, Mm=1.0, Rs=1.0)
        magnetic_moment.magnetic_moment(dynamo=dyn_region, planet=planet, star=star,jup=jup)
        magnetic_moment.magnetosphere_radius(mag_moment_jup, stellar_wind=stellar_wind)
        magnetic_moment.talk(talk=talk)
        if magnetic_moment.normalize_standoff_dist(planet) < 1 :
            magnetic_moment.standoff_dist = planet.unnormalize_radius()
            print('Warning ! Magnetosphere radius lower than 1.')
        mytarget = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": magnetic_moment},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": magnetic_moment,
                "stellar_wind": stellar_wind,
                "sw_jupiter" : sw_jup,
            },
            pow_received={"star": star},
            fmax_star={"star" : star},
        )
        mytarget.talk(talk=talk)
        if i == 0:
            
            df_target = pd.DataFrame(
                {
                    "name": mytarget.name,
                    "planet_mass": planet.mass,
                    "planet_radius": planet.radius,
                    "planet_luminosity": planet.luminosity,
                    "star_planet_distance": planet.stardist,
                    "planet_rotation_rate": planet.rotrate,
                    "planet_orbital_period": planet.orbitperiod,
                    "star_mass": star.mass,
                    "star_radius": star.radius,
                    "star_age": star.age,
                    "earth_distance": star.obs_dist,
                    "star_magfield" : star.magfield,
                    "star_rotperiod" : star.rotperiod,
                    "star_luminosity": star.luminosity,
                    "dynamo_density": dyn_region.density,
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
                    "pow_emission_kinetic": mytarget._pow_emission_kinetic / 1e14,
                    "pow_emission_magnetic": mytarget._pow_emission_magnetic / 1e14,
                    "flux kinetic au": mytarget.flux_kinetic_au / 1e-26 / 1e10,
                    "flux magnetic au": mytarget.flux_magnetic_au / 1e-26 / 1e10,
                    "pow_received_kinetic": mytarget._pow_received_kinetic * 1e3 / 1e-26,
                    "pow_received_magnetic": mytarget._pow_received_magnetic * 1e3 / 1e-26,
                    "star_magnetic_field" : star._magfield,
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
                    "planet_luminosity": planet.luminosity,
                    "star_planet_distance": planet.stardist,
                    "planet_rotation_rate": planet.rotrate,
                    "planet_orbital_period": planet.orbitperiod,
                    "star_mass": star.mass,
                    "star_radius": star.radius,
                    "star_age": star.age,
                    "earth_distance": star.obs_dist,
                    "star_magfield" : star.magfield,
                    "star_rotperiod" : star.rotperiod,
                    "star_luminosity": star.luminosity,
                    "dynamo_density": dyn_region.density ,
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
                    "pow_emission_kinetic": mytarget._pow_emission_kinetic / 1e14,
                    "pow_emission_magnetic": mytarget._pow_emission_magnetic / 1e14,
                    "flux kinetic au": mytarget.flux_kinetic_au / 1e-26 / 1e10,
                    "flux magnetic au": mytarget.flux_magnetic_au / 1e-26 / 1e10,
                    "pow_received_kinetic": mytarget._pow_received_kinetic * 1e3 / 1e-26,
                    "pow_received_magnetic": mytarget._pow_received_magnetic * 1e3 / 1e-26,
                    "star_magnetic_field" : star._magfield,
                    "freq_max_star": mytarget.freq_max_star,
                },
                index=[i],
            )
            df_target = pd.concat([df_target, df2], ignore_index=True)
        i += 1
    # if mytarget.select_target():
    # selected_targets.append(mytarget)

#df_target.to_csv("/Users/emauduit/Documents/These/Sélection des cibles/Programmes/M_VS_Mmag_result_RC.csv", sep=";", index=False)

## Sorting values by power of emission and by frequency

#df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic","freq_max"], ascending=[False,False])
#df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic"], ascending=[False])
#df_target_sorted.to_csv("/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/main_MSB_test_sorted.csv", sep=";", index=False)
#df_target_small=df_target_sorted[df_target_sorted['freq_max'] >= 4.99].iloc[0:50]
#df_target_small.to_csv("/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/main_MSB_top50.csv", sep=";", index=False)
