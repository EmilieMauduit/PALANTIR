#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

import os
import pandas as pd
import numpy as np
import datetime
import palantir
from importlib_resources import files

from palantir.prediction_tools.star import Star
from palantir.prediction_tools.planet import Planet
from palantir.prediction_tools.dynamo_region import DynamoRegion
from palantir.prediction_tools.magnetic_moment import MagneticMoment
from palantir.prediction_tools.stellar_wind import StellarWind
from palantir.prediction_tools.emission import Emission
from palantir.prediction_tools.target_selection import Config

import logging
log = logging.getLogger('palantir.scripts.main')

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

dateoftheday = datetime.datetime.today().isoformat().split(':')
dateofrun = dateoftheday[0] + 'h' + dateoftheday[1]
maps_dir = files("palantir.scripts.input_files")

os.system('mkdir '+config_param.output_path +'/'+dateofrun)

palantir.setup_logging(log_filepath=config_param.output_path + '/' + dateofrun + '/',verbose=config_param.talk)
log.info('This run was made with version {} of PALANTIR.'.format(palantir.__version__))
config_param.log_current_run_parameters()

# Criterions for target selection

selection_config = pd.read_csv(
    r"/Users/emauduit/Documents/These/target_selection/Programmes/selection.csv",
    delimiter=";",
)


# --------------------------------------------------------- #
# ---------------------- Data input ----------------------- #

dict_data = { 'nasa_data' : 'exoplanet_catalog_NASA.csv',
            'exoplanet_data' : 'exoplanet.eu_catalog.csv',
            'custom_data' : 'custom_dataset.csv'}


data =pd.read_csv(maps_dir / dict_data[config_param.database])
data = config_param.param_names(data=data)

# --------------------------------------------------------- #
# ------------------------ Main --------------------------- #
# --------------------------------------------------------- #

sun = Star(
    name="Sun",
    mass=1.0,
    radius={"models": config_param.star_radius_models, "radius": 1.0},
    age=AS,
    obs_dist=1.0,
    sp_type ='GV',
    )
jup = Planet(
    name="Jupiter",
    mass=1.0,
    radius={"models": config_param.planet_radius_models, "radius": 1.0},
    distance=5.2,
    eccentricity=0.0487,
    worb={"star_mass": MS, "worb": 1.0},
    luminosity={
        "models": config_param.planet_luminosity_models,
        "luminosity": np.nan,
        "star_age": 4.6,
    },
    wrot=1.0,
)
jup.tidal_locking(age=4.6e9, star_mass=1.0)
sun.compute_effective_temperature(np.nan)
sun.compute_magnetic_field(value= {'model': config_param.star_magfield_models, 'mag_field' : 1.435})
dyn_region_jup = DynamoRegion.from_planet(planet=jup, rhocrit=config_param.rho_crit)
dyn_region_jup.magnetic_field(planet=jup,rc_dyn=config_param.rc_dyn, jup=True)
mag_moment_jup = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.56e27, Rs=1.0)
vjup, vejup, nejup, Tjup = StellarWind._Parker(star=sun, planet=jup, T=0.81e6)
sw_jup = StellarWind(
    ne=nejup, ve=vejup, Tcor=Tjup, Bsw={"planet": jup, "star": sun, "vsw": vjup}
)

selected_targets = []

i = 1

skipped_targets = open('skipped_targets.txt', "w")
#df_target=pd.DataFrame()
#df_target['0'] = pd.Series(config_param.output_params_units, index = config_param.output_params)

new_rows = [config_param.output_params_units]

for target in data.itertuples():
    log.info("Planet : {}".format(target.pl_name))
    if ('PSR' in target.pl_name) or ('Pulsar' in target.pl_name) :
        print('Warning : {} has been skipped'.format(target.pl_name))
        skipped_targets.write(target.pl_name + ': pulsar \n')
        continue

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
                log.info("Stellar age is too small (< 100 yr) !!!")
                skipped_targets.write(target.pl_name + ': stellar age is too small\n')
                continue
            elif target.star_age < 0.7:
                star_age = 0.7
            else:
                star_age = target.star_age
        
        if np.isnan(target.orbital_period) :
            log.info("Planetary orbital period unknown")
            skipped_targets.write(target.pl_name + ': planetary orbital period unknown\n')
            continue

        planet = Planet(
            name=target.pl_name,
            mass=planet_mass,
            radius={"models": config_param.planet_radius_models, "radius": target.radius},
            distance=planet_distance,
            eccentricity=target.eccentricity,
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
            name=target.star_name,
            mass=target.star_mass,
            radius={"models": config_param.star_radius_models, "radius": target.star_radius},
            age=star_age,
            obs_dist=target.star_distance,
            sp_type = config_param.retrieve_spectral_type(star_name = target.star_name,sp_type = str(target.star_sp_type)),
            )

        if np.isnan(target.radius) and config_param.radius_expansion :
            planet.radius_expansion(luminosity=star.luminosity)
        try:
            planet.tidal_locking(age=star.age, star_mass=star.mass)
        except OverflowError:
            log.info("Divergence in tidal locking")
            skipped_targets.write(target.pl_name + ': divergence in tidal locking\n')
            continue
        
        if star.sp_type_code > config_param.sp_type_code :
            log.info("Star spectral type is not consistent with configuration parameters.")
            skipped_targets.write(target.pl_name + ': star spectral type is not consistent with configuration parameters.\n')
            continue

        #### Computing stellar magnetic field

        if config_param.star_magfield_catalog_only :
            catalog_Bstar = pd.read_csv(maps_dir / 'Bstar_catalog.csv', delimiter=';')
            crossmatch = catalog_Bstar[catalog_Bstar['Planet_Name']==target.pl_name]
            if (crossmatch.size < 1) or (config_param.star_magfield_catalog_only and not np.asarray(crossmatch['True_pred'])[0]):
                log.info("KNN prediction for B* could not be done since too many parameters were missing.")
                skipped_targets.write(target.pl_name + ': KNN prediction for B* could not be done since too many parameters were missing.\n')
                continue
            mag_field = np.asarray(crossmatch['B_G'])[0]

        else :
            catalog_Bstar = pd.read_csv(maps_dir / 'crossmatch_mag_exo.csv', delimiter=',')
            crossmatch = catalog_Bstar[catalog_Bstar['Simbad_ID']==star.main_id]
            mag_field = np.asarray(crossmatch['Bestim_G'])[0] if (crossmatch.size>0) else np.nan

        Teff = target.star_teff
        star.compute_effective_temperature(Teff)

        try :
            star.compute_magnetic_field(value= {'model': config_param.star_magfield_models, 'mag_field' : mag_field})
        except OverflowError :
            log.info("Divergence in stellar magnetic field estimate")
            skipped_targets.write(target.pl_name + ': divergence in stellar magnetic field estimate\n')
            continue

        dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=config_param.rho_crit)
        dyn_region.normalize(other=dyn_region_jup)
        dyn_region.magnetic_field(planet=planet, rc_dyn=config_param.rc_dyn)

        try:
            stellar_wind = StellarWind.from_system(star=star, planet=planet)
        except ValueError:
            log.info("Divergence in stellar wind calculation")
            skipped_targets.write(target.pl_name + ' : divergence in stellar wind calculation\n')
            continue

        magnetic_moment = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.0, Rs=1.0)
        magnetic_moment.magnetic_moment(
            dynamo=dyn_region, planet=planet, jup=jup, dynamo_jup = dyn_region_jup
        )
        magnetic_moment.magnetosphere_radius(mag_moment_jup, stellar_wind=stellar_wind)
        if magnetic_moment.normalize_standoff_dist(planet) < 1:
            magnetic_moment.standoff_dist = planet.unnormalize_radius()
            log.info("Magnetosphere radius lower than 1.")
            
        target_emission = Emission(
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
        
        if config_param.talk :
            print(planet)
            print(star)
            print(stellar_wind)
            print(dyn_region)
            print(magnetic_moment)
            print(target_emission)
        
        new_rows.append([
            target_emission.name,
            target.ra,
            target.dec,
            planet.mass,
            planet.radius,
            planet.luminosity,
            planet.stardist,
            planet.rotrate,
            planet.orbitperiod,
            star.main_id,
            star.mass,
            star.radius,
            star.age,
            star.obs_dist,
            star.magfield,
            star.rotperiod,
            star.luminosity,
            star.sp_type,
            star.sp_type_code,
            star.effective_temperature,
            dyn_region.density,
            dyn_region.radius, #/ planet.unnormalize_radius(),
            dyn_region.mag_field_dynamo,
            dyn_region.mag_field_equatorial,
            magnetic_moment.mag_moment,
            magnetic_moment.normalize_standoff_dist(planet=planet),
            stellar_wind.density,
            stellar_wind.effective_velocity,
            stellar_wind.corona_temperature,
            stellar_wind.mag_field,
            stellar_wind.alfven_velocity,
            target_emission._mag_field_planet,
            target_emission._freq_max_planet / 1e6,
            target_emission._pow_emission_kinetic / 1e14,
            target_emission._pow_emission_magnetic / 1e14,
            target_emission._pow_emission_spi / 1e14,
            target_emission.flux_kinetic_au / 1e-26 / 1e10,
            target_emission.flux_magnetic_au / 1e-26 / 1e10,
            target_emission.flux_spi_au / 1e-26 / 1e10,
            target_emission._pow_received_kinetic* 1e3/ 1e-26,
            target_emission._pow_received_magnetic* 1e3/ 1e-26,
            target_emission._pow_received_spi* 1e3/ 1e-26,
            target_emission.freq_max_star/ 1e6,
        ]#, index = config_param.output_params)
        )
        i+=1
    
    else :
        skipped_targets.write(target.pl_name + ' : semi-major axis, planetary mass or star mass unknown.\n')
    # if mytarget.select_target():
    # selected_targets.append(mytarget)
print(sun)
print(sw_jup)

skipped_targets.close()
#df_target = pd.concat([df_target]+new_columns, axis=1)
#df_target = df_target.transpose()
df_target = pd.DataFrame(new_rows,columns=config_param.output_params)
# --------------------------------------------------------- #
# -------- Saving input and output in one folder  --------- #

os.system('cp skipped_targets.txt '+config_param.output_path +'/'+dateofrun+'/skipped_targets.txt')
os.system('rm skipped_targets.txt')

data.to_csv(config_param.output_path +"/"+dateofrun+"/catalog_input.csv", sep=";", index=False)
df_target.to_csv(config_param.output_path +"/"+dateofrun+"/main_output.csv", sep=";", index=False)

# Sorting values by power of emission and by frequency

#df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic","freq_max"], ascending=[False,False])
#df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic"], ascending=[False])
#df_target_sorted.to_csv("/Users/emauduit/Documents/These/target_selection/Programmes/MvsMmag_rcdyn_sorted.csv", sep=";", index=False)
# df_target_small=df_target_sorted[df_target_sorted['freq_max'] >= 4.99].iloc[0:50]
# df_target_small.to_csv("/Users/emauduit/Documents/These/target_selection/Programmes/main_MSB_top50.csv", sep=";", index=False)



