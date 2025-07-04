#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

import pandas as pd
from astroquery.simbad import Simbad
from importlib_resources import files

import logging
logger = logging.getLogger('palantir.prediction_tools.target_selection')

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

# ============================================================= #
# --------------------------- Config -------------------------- #
# ============================================================= #
class Config:
    def __init__(self, config_file : str = None):
        if config_file is None :
            maps_dir = files("palantir.scripts.input_files")
            config = pd.read_csv(maps_dir / "run_parameters.csv", delimiter=';')
        else : 
            config = pd.read_csv(config_file, delimiter=';')

        database = [config.setting[i] for i in range(3) if int(config.value[i]) == 1]
        if len(database) > 1 :
            logger.error('ValueError : Only one database can be used, two were given in the configuration file')
            raise ValueError('Only one database can be used, two were given in the configuration file')
        else :
            self.database = database[0]
        self.magnetic_moment_models = [config.setting[i] for i in range(3,10) if int(config.value[i]) == 1]
        self.dynamo_density_models = [config.setting[i] for i in range(10, 12) if int(config.value[i]) == 1]
        self.planet_radius_models = [config.setting[i] for i in range(13, 15) if int(config.value[i]) == 1]
        self.star_radius_models = [config.setting[i] for i in range(16, 17) if int(config.value[i]) == 1]
        self.planet_luminosity_models = [
            config.setting[i] for i in range(17, 20) if int(config.value[i]) == 1
        ]
        self.star_magfield_models = [config.setting[i] for i in range(21, 23) if int(config.value[i]) == 1]
        self.star_magfield_catalog_only = True if (int(config.value[23]) == 1) else False
        self.rho_crit = int(config.value[12])
        self.rc_dyn = True if int(config.value[6]) == 1 else False
        self.radius_expansion = True if int(config.value[15]) == 1 else False
        self.sp_type_code = int(config.value[20])
        self.talk = True if int(config.value[24]) == 1 else False
        self.output_path = config.value[25]
        self.output_params = ["name","ra","dec","planet_mass", "planet_radius", "planet_luminosity", "star_planet_distance",
            "planet_rotation_rate", "planet_orbital_period","star_simbad_id", "star_mass","star_radius","star_age","earth_distance",
            "star_magfield","star_rotperiod","star_luminosity","spectral_type", "spectral_type_code","star_effective_temp","dynamo_density",
            "dynamo_radius","B_dyn" ,"B_eq","magnetic_moment","standoff_distance","sw_density","sw_velocity",
            "coronal_temperature","sw_magfield","alfven_velocity","magnetic_field_planet","freq_max_planet",
            "pow_emission_kinetic","pow_emission_magnetic","pow_emission_spi", "flux_kinetic_au", "flux_magnetic_au",
            "flux_spi_au", "pow_received_kinetic","pow_received_magnetic", "pow_received_spi",
            "freq_max_star"]
        self.output_params_units = ["","hh:mm:ss", "dd:mm:ss", "MJ", "RJ", "LS", "AU","wJ","w_orb_J","", "MS", "RS", "yr", "pc", "G", "days",
            "LS", "", "","K", "rho_dyn_J", "r_dyn_J", "T", "T", "MmagJ", "Rp", "m-3", "m.s-1", "K", "T", "m.s-1", "T", "MHz", "10^14W", 
            "10^14W", "10^14W","10^10Jy", "10^10Jy", "10^10Jy", "mJy","mJy", "mJy", "MHz"]
        logger.info('Configuration parameters succesfully initialized.')

        

    ############ Methods ############

    def param_names(self,data : pd.DataFrame ) :
        if self.database == 'nasa_data' :
            new_names = {'pl_massj':'mass','pl_massjerr2':'mass_error_min','pl_massjerr1':'mass_error_max',
                'pl_msinij':'mass_sini','pl_msinijerr2':'mass_sini_error_min','pl_msinijerr1':'mass_sini_error_max',
                'pl_radj':'radius','pl_radjerr2':'radius_error_min','pl_radjerr1':'radius_error_max','pl_orbper':'orbital_period',
                'pl_orbpererr2':'orbital_period_error_min','pl_orbpererr1':'orbital_period_error_max','pl_orbsmax':'semi_major_axis',
                'pl_orbsmaxerr2':'semi_major_axis_error_min','pl_orbsmaxerr1':'semi_major_axis_error_max','pl_orbeccen':'eccentricity',
                'pl_orbeccenerr2':'eccentricity_error_min','pl_orbeccenerr1':'eccentricity_error_max','pl_orbincl':'inclination',
                'pl_orbinclerr2':'inclination_error_min','pl_orbinclerr1':'inclination_error_max','disc_year':'discovered',
                'rowupdate':'updated','pl_imppar':'impact_parameter','pl_impparerr2':'impact_parameter_error_min',
                'pl_impparerr1':'impact_parameter_error_max','st_logg':'log_g','discoverymethod':'detection_type',
                'pl_bmassprov':'mass_detection_type','hostname':'star_name','rastr':'ra','decstr':'dec','sy_vmag':'mag_v','sy_kmag':'mag_k',
                'sy_dist':'star_distance','sy_disterr2':'star_distance_error_min','sy_disterr1':'star_distance_error_max',
                'st_met':'star_metallicity','st_meterr2':'star_metallicity_error_min','st_meterr1':'star_metallicity_error_max',
                'st_mass':'star_mass','st_masserr2':'star_mass_error_min','st_masserr1':'star_mass_error_max','st_rad':'star_radius',
                'st_raderr2':'star_radius_error_min','st_raderr1':'star_radius_error_max','st_spectype':'star_sp_type','st_age':'star_age',
                'st_ageerr2':'star_age_error_min','st_ageerr1':'star_age_error_max','st_teff':'star_teff','st_tefferr2':'star_teff_error_min',
                'st_tefferr1':'star_teff_error_max'}
            data = data.rename(columns=new_names)
        elif self.database == 'exoplanet_data':
            new_names={'name' : 'pl_name'}
            data = data.rename(columns=new_names)
        return data

    def retrieve_spectral_type(self,star_name,sp_type):
    # Configure Simbad to display the spectral type
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('sptype')

        # Query Simbad for the star
        result_table = custom_simbad.query_object(star_name)

        # Check if the result is not empty and contains the spectral type
        if result_table is not None and 'SP_TYPE' in result_table.colnames:
            spectral_type = result_table['SP_TYPE'][0]
            if spectral_type:
                return spectral_type
            else:
                return sp_type
        else:
            return sp_type

    def log_current_run_parameters(self):
        logger.info("Database used for this run : {}".format(self.database))
        logger.info("Models used for planetary magnetic moment predictions : {}".format(self.magnetic_moment_models))
        logger.info("Models used for planetary dynamo density : {}".format(self.dynamo_density_models))
        logger.info("Critical density value used : {} g.cm-3".format(self.rho_crit))
        logger.info("Model used for planetary radius prediction : {}".format(self.planet_radius_models))
        logger.info("Was planetary radius expansion used ? {}".format(self.radius_expansion))
        logger.info("Model used for stellar radius prediction : {}".format(self.star_radius_models))
        logger.info("Tables used for planetary apparent luminosity : {}".format(self.planet_luminosity_models))
        logger.info("Star spectral type criterion : {}".format(self.sp_type_code))
        logger.info("Models used for stellar magnetic field prediction : {}".format(self.star_magfield_models))
        logger.info("Were stellar magnetic_field based only on catalog ? {}".format(self.star_magfield_catalog_only))
        logger.info("Path given to store the outputs : {}".format(self.output_path))



# ============================================================= #
# ------------------------- Prediction ------------------------ #
# ============================================================= #

class Prediction:

    def __init__(self):
        pass
