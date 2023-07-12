#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

from warnings import WarningMessage
import pandas as pd
import numpy as np
from star import Star
from planet import Planet
from dynamo_region import DynamoRegion
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind, parker
from target import Target
from target_selection import Config

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

# rc = 0.85 * RJ
# rhoc = 1800  # kg/m3


# --------------------------------------------------------- #
# -------------- Configuration settings input ------------- #

config_param = Config()
print(config_param.rc_dyn)

# Criterions for target selection

selection_config = pd.read_csv(
    r"/Users/emauduit/Documents/These/Sélection des cibles/Programmes/selection.csv",
    delimiter=";",
)


# --------------------------------------------------------- #
# ---------------------- Data input ----------------------- #

#data =pd.read_csv(r'/Users/emauduit/Documents/These/Sélection des cibles/Programmes/exoplanet.eu_catalog.csv')
#, index_col=0)
#data =pd.read_csv(r'/Users/emauduit/Documents/These/Articles/Proceeding PRE/M_VS_Mmag.csv')#, index_col=0)
data = pd.read_csv(
    r"/Users/emauduit/Documents/These/Sélection des cibles/Programmes/planet_test_unique.csv",
    delimiter=";")

data = config_param.param_names(data=data)

# --------------------------------------------------------- #
# ------------------------ Main --------------------------- #
# --------------------------------------------------------- #

sun = Star(
    name="Soleil",
    mass=1.0,
    radius={"models": config_param.star_radius_models, "radius": 1.0},
    age=AS,
    obs_dist=1.0,
)
# sun.talk(talk=talk)
jup = Planet(
    name="Jupiter",
    mass=1.0,
    radius={"models": config_param.planet_radius_models, "radius": 1.0},
    distance=5.2,
    worb={"star_mass": MS, "worb": 1.0},
    luminosity={
        "models": config_param.planet_luminosity_models,
        "luminosity": np.nan,
        "star_age": 4.6,
    },
    wrot=1.0,
)
jup.talk(talk=config_param.talk)
jup.tidal_locking(age=4.6e9, star_mass=1.0)
dyn_region_jup = DynamoRegion.from_planet(planet=jup, rhocrit=config_param.rho_crit)
dyn_region_jup.magnetic_field(planet=jup,rc_dyn=config_param.rc_dyn, jup=True)
mag_moment_jup = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.56e27, Rs=1.0)
# mag_moment_jup.talk(talk=True)
vjup, vejup, nejup, Tjup = parker(star=sun, planet=jup, T=0.81e6)
sw_jup = StellarWind(
    ne=nejup, ve=vejup, Tcor=Tjup, Bsw={"planet": jup, "star": sun, "vsw": vjup}
)
sw_jup.talk(talk=config_param.talk)
selected_targets = []

i = 0

skipped_targets = []

for target in data.itertuples():
    print(target.pl_name)
    #if 'PSR' in target.pl_name :
    #    print('Warning : {} has been skipped'.format(target.pl_name))
    #    skipped_targets.append(target.pl_name)
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
            if target.star_age <= 1e-7:
                print("Stellar age is too small (< 100 yr) !!!")
                skipped_targets.append(target.pl_name)
                continue
            elif target.star_age < 0.5:
                star_age = 0.5
            else:
                star_age = target.star_age

        planet = Planet(
            name=target.pl_name,
            mass=planet_mass,
            radius={"models": config_param.planet_radius_models, "radius": target.radius},
            distance=planet_distance,
            worb={"star_mass": target.star_mass, "worb": np.nan},
            luminosity={
                "models": config_param.planet_luminosity_models,
                "luminosity": np.nan,
                "star_age": star_age,
            },
            detection_method=target.detection_type,
            wrot=1.0,
        )
        star = Star(
            name="star",
            mass=target.star_mass,
            radius={"models": config_param.star_radius_models, "radius": target.star_radius},
            age=star_age,
            obs_dist=target.star_distance,
        )

        if np.isnan(target.radius):
            planet.radius_expansion(
                luminosity=star.luminosity, eccentricity=target.eccentricity
            )
        try:
            planet.tidal_locking(age=star.age, star_mass=star.mass)
        except OverflowError:
            skipped_targets.append(target.pl_name)
            continue
        planet.talk(talk=config_param.talk)
        star.talk(talk=config_param.talk)
        dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=config_param.rho_crit)
        dyn_region.normalize(other=dyn_region_jup)
        dyn_region.magnetic_field(planet=planet, rc_dyn=config_param.rc_dyn)
        dyn_region_jup.talk(talk=True)
        dyn_region.talk(talk=config_param.talk)
        try:
            stellar_wind = StellarWind.from_system(star=star, planet=planet)
            stellar_wind.talk(talk=config_param.talk)
        except ValueError:
            skipped_targets.append(target.pl_name)
            continue
        magnetic_moment = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.0, Rs=1.0)
        magnetic_moment.magnetic_moment(
            dynamo=dyn_region, planet=planet, jup=jup, dynamo_jup = dyn_region_jup
        )
        magnetic_moment.magnetosphere_radius(mag_moment_jup, stellar_wind=stellar_wind)
        magnetic_moment.talk(talk=config_param.talk)
        if magnetic_moment.normalize_standoff_dist(planet) < 1:
            magnetic_moment.standoff_dist = planet.unnormalize_radius()
            print("Magnetosphere radius lower than 1.")
            
        mytarget = Target(
            name=planet.name,
            mag_field={"planet": planet, "magnetic_moment": magnetic_moment},
            pow_emission={
                "planet": planet,
                "star": star,
                "magnetic_moment": magnetic_moment,
                "stellar_wind": stellar_wind,
                "sw_jupiter": sw_jup,
            },
            pow_received={"star": star},
            fmax_star={"star": star},
        )
        mytarget.talk(talk=config_param.talk)
        if i == 0:

            df_target = pd.DataFrame(
                {
                    "name": mytarget.name,
                    #"ra" : target.ra,
                    #"dec" : target.dec,
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
                    "star_magfield": star.magfield,
                    "star_rotperiod": star.rotperiod,
                    "star_luminosity": star.luminosity,
                    "dynamo_density": dyn_region.density,
                    "dynamo_radius": dyn_region.radius, #/ planet.unnormalize_radius(),
                    "B_dyn" : dyn_region.mag_field_dynamo,
                    "B_eq" : dyn_region.mag_field_equatorial,
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
                    "pow_received_kinetic": mytarget._pow_received_kinetic
                    * 1e3
                    / 1e-26,
                    "pow_received_magnetic": mytarget._pow_received_magnetic
                    * 1e3
                    / 1e-26,
                    "star_magnetic_field": star._magfield,
                    "freq_max_star": mytarget.freq_max_star,
                },
                index=[i],
            )
        else:
            df2 = pd.DataFrame(
                {
                    "name": mytarget.name,
                    #"ra" : target.ra,
                    #"dec" : target.dec,
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
                    "star_magfield": star.magfield,
                    "star_rotperiod": star.rotperiod,
                    "star_luminosity": star.luminosity,
                    "dynamo_density": dyn_region.density,
                    "dynamo_radius": dyn_region.radius, #/ planet.radius,
                    "B_dyn" : dyn_region.mag_field_dynamo,
                    "B_eq" : dyn_region.mag_field_equatorial,
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
                    "pow_received_kinetic": mytarget._pow_received_kinetic
                    * 1e3
                    / 1e-26,
                    "pow_received_magnetic": mytarget._pow_received_magnetic
                    * 1e3
                    / 1e-26,
                    "star_magnetic_field": star._magfield,
                    "freq_max_star": mytarget.freq_max_star,
                },
                index=[i],
            )
            df_target = pd.concat([df_target, df2], ignore_index=True)
        i += 1
    
    else :
        skipped_targets.append(target.pl_name)
    # if mytarget.select_target():
    # selected_targets.append(mytarget)

df_target.to_csv("/Users/emauduit/Documents/These/Sélection des cibles/Programmes/main_RC_test.csv", sep=";", index=False)

# Sorting values by power of emission and by frequency

#df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic","freq_max"], ascending=[False,False])
df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic"], ascending=[False])
df_target_sorted.to_csv("/Users/emauduit/Documents/These/Sélection des cibles/Programmes/main_RC_test_sorted.csv", sep=";", index=False)
# df_target_small=df_target_sorted[df_target_sorted['freq_max'] >= 4.99].iloc[0:50]
# df_target_small.to_csv("/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/main_MSB_top50.csv", sep=";", index=False)

print(skipped_targets)