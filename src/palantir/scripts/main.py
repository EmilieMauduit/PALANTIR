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
import html
import datetime
import palantir
from importlib_resources import files
from astropy.coordinates import SkyCoord
import astropy.units as u
from pyinstrument import Profiler

from palantir.prediction_tools.star import Star
from palantir.prediction_tools.planet import Planet
from palantir.prediction_tools.dynamo_region import DynamoRegion
from palantir.prediction_tools.magnetic_moment import MagneticMoment
from palantir.prediction_tools.stellar_wind import ParkerGrid, StellarWind
from palantir.prediction_tools.emission import Emission
from palantir.prediction_tools.target_selection import Config,FlagTracker, XRayFluxCalculatorNEXXUS

import logging
log = logging.getLogger('palantir.scripts.main')

# --------------------------------------------------------- #
# ------------------- Physical constants ------------------ #

MS = 1.989e30  # kg
RS = 6.96342e8  # m
AS = 4.6  # Gyo
BSsw = 1  # T
LS = 3.826e26  # W


MJ = 1.8986e27  # kg
RJ = 69911e3  # m
wrotJ = 1.77e-4  # rad.s-1
worbJ = 1.68e-8  # rad.s-1

ME = 5.97237e24  # kg
RE = 6371.0e3  # m
wE = 7.27e-5  # s-1

dua = 1.49597870700e11  # m
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

Xray_calculator = XRayFluxCalculatorNEXXUS()
parker_grid = ParkerGrid()
if config_param.star_magfield_models[0] == 'Bstar_catalog' :
    catalog_Bstar = pd.read_csv(maps_dir / 'Bstar_catalog_latest.csv', delimiter=';')

# --------------------------------------------------------- #
# ---------------------- Data input ----------------------- #
# --------------------------------------------------------- #

dict_data = { 'nasa_data' : 'exoplanet_catalog_NASA.csv',
            'exoplanet_data' : 'exoplanet.eu_catalog.csv',
            'custom_data' : 'custom_dataset.csv'}


data = pd.read_csv(maps_dir / dict_data[config_param.database])
data = config_param.param_names(data=data)
os.system('cp ' + str(maps_dir / dict_data[config_param.database]) + ' ' + config_param.output_path +"/"+dateofrun+"/catalog_input.csv")

# --------------------------------------------------------- #
# ------------------------ Main --------------------------- #
# --------------------------------------------------------- #


# Instanciating stellar and planetary references as the Sun and Jupiter

sun = Star(
    name="Sun",
    main_id="Sun",
    mass={"mass" : 1.0, "radius" : 1.0, "sptype" : "GV"},
    coordinates= SkyCoord(ra=286.123*u.deg, dec=63.87*u.deg, frame='icrs'),
    radius={"models": config_param.star_radius_models, "radius": 1.0},
    age=AS,
    Teff = np.nan,
    rot_period={"period": np.nan, "vsini" : np.nan, "flux_B" : np.nan, "flux_V" : np.nan},
    obs_dist=1.0,
    sp_type ='GV',
    Xray_calculator=Xray_calculator,
    magnetic_field={'model': config_param.star_magfield_models, 'mag_field' : 1.435, 'Bstar_catalog' : config_param.Bstar_database}
    )
jupiter  = Planet(
    name="Jupiter",
    mass={"mass" : 1.0, "mass_sini" : np.nan, "radius" : 1.0},
    radius={"models": config_param.planet_radius_models, "radius": 1.0},
    semi_major_axis=5.4546,
    distance=5.2,
    eccentricity=0.0487,
    Torb={"star_mass": 1., "Torb": 4332.01},
    luminosity={
        "models": config_param.planet_luminosity_models,
        "luminosity": np.nan,
        "star_age": 4.6,
    },
    Trot=0.41351,
)
jupiter.tidal_locking(age=4.6e9, star_mass=1.0)
dyn_region_jup = DynamoRegion.from_planet(planet=jupiter, rhocrit=config_param.rho_crit)
dyn_region_jup.magnetic_field(planet=jupiter,rc_dyn=config_param.rc_dyn, jup=True)
mag_moment_jup = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.56e27, Rm=1.0)
sw_jup = StellarWind(Tcor={"star" : sun, "Tcor" : 1.77e6}, d_alfven_point={ "star": sun, "eccentricity" : jupiter.eccentricity, "parker_grid" : parker_grid})
sw_jup.compute_Parker_solution(planet = jupiter, star=sun, parker_grid=parker_grid)
sw_jup.compute_B_imf_components(planet = jupiter, star=sun)

selected_targets = []

i = 1

skipped_targets = open('skipped_targets.txt', "w")
new_rows = [config_param.output_params_units]

profiling_pyinstrument = False
if profiling_pyinstrument:
    profiler = Profiler()
    profiler.start()

try :

    for target in data.itertuples():
        flag = 0
        log.info("Planet : {}".format(target.pl_name))
        if ('PSR' in target.pl_name) or ('Pulsar' in target.pl_name) :
            print('Warning : {} has been skipped'.format(target.pl_name))
            skipped_targets.write(target.pl_name + ': pulsar \n')
            continue

        if (not np.isnan(target.semi_major_axis)):
            flag_tracker = FlagTracker(nmax_hypothesis=7)
            if np.isnan(target.mass) and np.isnan(target.mass_sini) and np.isnan(target.radius):
                log.info("Planetary mass and radius were not available, those parameters were set equal to Jupiter's.")
                flag_tracker.activate(hypothesis_number=0)

            planet_distance = target.semi_major_axis if np.isnan(target.eccentricity) else target.semi_major_axis * (1 - target.eccentricity)

            if np.isnan(target.star_age):
                star_age = 5.2
                flag_tracker.activate(hypothesis_number=1)
            else:
                if target.star_age <= 0.0001 :
                    log.info("Stellar age is too small")
                    skipped_targets.write(target.pl_name + ': stellar age is too small\n')
                    continue
                else:
                    star_age = target.star_age
            
            # ==================================================== #
            # ---------------------- Planet ---------------------- #
            # ==================================================== #

            planet = Planet(
                name=target.pl_name,
                mass={"mass" : target.mass, "mass_sini" : target.mass_sini, "radius" : target.radius},
                radius={"models": config_param.planet_radius_models, "radius": target.radius},
                semi_major_axis=target.semi_major_axis,
                distance=planet_distance,
                eccentricity=target.eccentricity,
                Torb={"star_mass": target.star_mass, "Torb": target.orbital_period},
                luminosity={
                    "models": config_param.planet_luminosity_models,
                    "luminosity": np.nan,
                    "star_age": star_age,
                },
                detection_method=target.detection_type,
                Trot=jupiter.rotperiod,
            )

            # ================================================== #
            # ---------------------- Star ---------------------- #
            # ================================================== #

            simbad_query_result = config_param.query_simbad_star_param(star_name=html.unescape(target.star_name), sp_type = str(target.star_sp_type), star_alternate_names = target.star_alternate_names)

            if config_param.star_magfield_models[0] == 'Bstar_catalog' :
                crossmatch = catalog_Bstar[catalog_Bstar['Planet_Name']==target.pl_name]
                if (crossmatch.size < 1):
                    log.info("KNN prediction for B* could not be done since too many parameters were missing.")
                    skipped_targets.write(target.pl_name + ': KNN prediction for B* could not be done since too many parameters were missing.\n')
                    continue
                mag_field = np.asarray(crossmatch['B_G'])[0]
            else : 
                mag_field = np.nan

            try : 
                star = Star(
                    name=html.unescape(target.star_name),
                    main_id= simbad_query_result['main_id'],
                    coordinates = SkyCoord(ra=target.ra*u.deg, dec=target.dec*u.deg, frame='icrs'),
                    mass={"mass" : target.star_mass, "radius" : target.star_radius, "sptype" : simbad_query_result['sp_type']},
                    radius={"models": config_param.star_radius_models, "radius": target.star_radius},
                    age=star_age,
                    Teff = target.star_teff,
                    rot_period = simbad_query_result,
                    obs_dist=target.star_distance,
                    sp_type = simbad_query_result['sp_type'],
                    Xray_calculator=Xray_calculator,
                    magnetic_field={'model': config_param.star_magfield_models, 'mag_field' : mag_field, 'Bstar_catalog' : config_param.Bstar_database}
                    )
            except (OverflowError, ValueError) :
                log.info("Divergence in stellar magnetic field estimate")
                skipped_targets.write(target.pl_name + ': divergence in stellar magnetic field estimate\n')
                continue
            
            if np.isnan(star.mass) :
                log.info("Unable to estimate stellar mass.")
                skipped_targets.write(target.pl_name + ': stellar mass could not be estimated.\n')
                continue
            
            if np.isnan(target.star_mass) and np.isnan(target.star_radius) :
                log.info("No value for stellar mass or radius in the catalog. Stellar spectral type was used to estimate the mass.")
                flag_tracker.activate(hypothesis_number = 2)

            if star.rotperiod < 0.042 : #rotation period of less than 1h is unlikely.
                log.info("Stellar rotation period is less than 1h.")
                flag_tracker.activate(hypothesis_number = 3)

            if star.sp_type_code > config_param.sp_type_code :
                log.info("Star spectral type is not consistent with configuration parameters.")
                skipped_targets.write(target.pl_name + ': star spectral type is not consistent with configuration parameters {}.\n'.format(target.star_sp_type))
                continue

            if np.isnan(target.radius) and config_param.radius_expansion :
                planet.radius_expansion(luminosity=star.luminosity)

            if (np.isnan(target.eccentricity)) or (target.eccentricity == 0.) :
                try:
                    planet.tidal_locking(age=star.age, star_mass=star.mass)
                except OverflowError:
                    log.info("Divergence in tidal locking")
                    skipped_targets.write(target.pl_name + ': divergence in tidal locking\n')
                    continue
            else :
                planet.tidally_locked = False
            
            if planet.rotrate == jupiter.rotrate :
                log.info("No tidal locking for this planet, rotation period assumed to be Jupiter's.")
                flag_tracker.activate(hypothesis_number = 4)

            # =========================================================== #
            # ---------------------- Dynamo Region ---------------------- #
            # =========================================================== #

            try : 
                dyn_region = DynamoRegion.from_planet(planet=planet, rhocrit=config_param.rho_crit)
            except IndexError :
                log.info("No solution found for LaneEmden equation.")
                dyn_region = DynamoRegion(rhocrit = config_param.rho_crit, rhoc = np.nan, rc = np.nan)
                flag_tracker.activate(hypothesis_number=5)

            dyn_region.magnetic_field(planet=planet, rc_dyn=config_param.rc_dyn)

            # ========================================================== #
            # ---------------------- Stellar Wind ---------------------- #
            # ========================================================== #

            try:
                stellar_wind = StellarWind(Tcor={"star" : star, "Tcor" : np.nan}, d_alfven_point={ "star": star, "eccentricity" : planet.eccentricity, "parker_grid" : parker_grid})
                stellar_wind.compute_Parker_solution(planet = planet, star=star, parker_grid=parker_grid)
                stellar_wind.compute_B_imf_components(planet = planet, star=star)
            except (ValueError, RuntimeError):
                log.info("Divergence in stellar wind calculation")
                skipped_targets.write(target.pl_name + ': divergence in stellar wind calculation\n')
                continue

            # ============================================================= #
            # ---------------------- Magnetic Moment ---------------------- #
            # ============================================================= #

            magnetic_moment = MagneticMoment(models=config_param.magnetic_moment_models, Mm=1.0, Rm=1.0)
            magnetic_moment.magnetic_moment(
                dynamo=dyn_region, planet=planet, jup=jupiter, dynamo_jup = dyn_region_jup
            )
            magnetic_moment.calc_magnetosphere_radius(mag_moment_jup, stellar_wind=stellar_wind)
            if magnetic_moment.normalize_standoff_dist(planet) < 1:
                magnetic_moment.magnetosphere_radius = planet.unnormalize_radius()
                log.info("Predicted magnetosphere radius lower than 1 Rp. Forced to 1. Rp")
                flag_tracker.activate(hypothesis_number=6)

            if not config_param.rc_dyn :
                dyn_region.mag_field_equatorial = magnetic_moment.mag_moment * mag_moment_jup.mag_moment * 1e-7 /np.power(planet.unnormalize_radius(),3)
                dyn_region.mag_field_dynamo = dyn_region.mag_field_equatorial * 2 * np.sqrt(2)

            # ====================================================== #
            # ---------------------- Emission ---------------------- #
            # ====================================================== #

            target_emission = Emission(
                name=planet.name,
                planetary_magnetic_field={"planet": planet, "magnetic_moment": magnetic_moment},
                pow_emission={
                    "planet": planet,
                    "star": star,
                    "magnetic_moment": magnetic_moment,
                    "stellar_wind": stellar_wind,
                    "sw_jupiter": sw_jup,
                },
                flux_received={"star": star, "stellar_wind" : stellar_wind},
                free_free_spi={"star": star, "stellar_wind" : stellar_wind},
                free_free_ms={"star": star, "stellar_wind" : stellar_wind, "planet":planet},
                fcmax_star={"star": star, "planet": planet, "stellar_wind" : stellar_wind, "T_corona": stellar_wind.corona_temperature, "parker_grid":parker_grid},
                fp_planet={"ne" : stellar_wind.density_planet},
                fp_star={"ne" : stellar_wind.density_star}
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
                planet.semi_major_axis,
                planet.rotperiod * 24.0,
                planet.orbitperiod,
                planet.tidally_locked,
                star.main_id,
                star.mass,
                star.radius,
                star.age,
                star.obs_dist,
                star.magfield,
                star.rotperiod,
                star.luminosity,
                star.Xray_flux,
                star.sp_type,
                star.sp_type_code,
                star.effective_temperature,
                dyn_region.density,
                dyn_region.radius / planet.unnormalize_radius(),
                dyn_region.mag_field_dynamo,
                dyn_region.mag_field_equatorial,
                magnetic_moment.mag_moment,
                magnetic_moment.normalize_standoff_dist(planet=planet),
                stellar_wind.density_planet,
                stellar_wind.effective_velocity,
                stellar_wind.velocity_sw,
                stellar_wind.mass_loss_rate,
                stellar_wind.corona_temperature,
                stellar_wind.radial_mag_field,
                stellar_wind.azimuthal_mag_field,
                stellar_wind.total_mag_field,
                stellar_wind.perp_mag_field,
                stellar_wind.distance_alfven_point,
                stellar_wind.alfven_velocity,
                target_emission._mag_field_planet,
                target_emission._freq_c_max_planet / 1e6,
                target_emission._freq_p_planet / 1e6,
                target_emission._pow_emission_kinetic,
                target_emission._pow_emission_magnetic,
                target_emission._pow_emission_spi,
                target_emission.flux_kinetic_au / 1e-26,
                target_emission.flux_magnetic_au / 1e-26,
                target_emission.flux_spi_au / 1e-26,
                target_emission._flux_received_kinetic* 1e3/ 1e-26,
                target_emission._flux_received_magnetic* 1e3/ 1e-26,
                target_emission._flux_received_spi* 1e3/ 1e-26,
                target_emission.freq_c_max_star/ 1e6,
                target_emission.freq_p_star/ 1e6,
                target_emission.test_escaping_spi["radius"],
                target_emission.test_escaping_spi['density'],
                target_emission.test_escaping_spi["fc_star"]/1e6,
                target_emission.test_escaping_spi["fp_star"]/1e6,
                target_emission.tau_free_free_ms,
                target_emission.tau_free_free_spi,
                flag_tracker.flags
            ]#, index = config_param.output_params)
            )
            i+=1
        
        else :
            skipped_targets.write(target.pl_name + ': semi-major axis unknown.\n')
            continue

except KeyboardInterrupt : 
    log.info("Interrupted by user.")

finally :
    skipped_targets.close()
    #df_target = pd.concat([df_target]+new_columns, axis=1)
    #df_target = df_target.transpose()
    df_target = pd.DataFrame(new_rows,columns=config_param.output_params)
    # --------------------------------------------------------- #
    # -------- Saving input and output in one folder  --------- #

    os.system('cp skipped_targets.txt '+config_param.output_path +'/'+dateofrun+'/skipped_targets.txt')
    os.system('rm skipped_targets.txt')


    #data.to_csv(config_param.output_path +"/"+dateofrun+"/catalog_input.csv", sep=",", index=False)
    df_target.to_csv(config_param.output_path +"/"+dateofrun+"/main_output.csv", sep=";", index=False)

    if profiling_pyinstrument:
            profiler.stop()
            profiler.print()

    # Sorting values by power of emission and by frequency

    #df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic","freq_max"], ascending=[False,False])
    #df_target_sorted=df_target.sort_values(by = ["pow_received_magnetic"], ascending=[False])
    #df_target_sorted.to_csv("/Users/emauduit/Documents/These/target_selection/Programmes/MvsMmag_rcdyn_sorted.csv", sep=";", index=False)
    # df_target_small=df_target_sorted[df_target_sorted['freq_max'] >= 4.99].iloc[0:50]
    # df_target_small.to_csv("/Users/emauduit/Documents/These/target_selection/Programmes/main_MSB_top50.csv", sep=";", index=False)



