#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

from operator import methodcaller
import pandas as pd
import numpy as np
from star import Star
from planet import Planet
from dynamo_region import DynamoRegion
from magnetic_moment import MagneticMoment
from stellar_wind import StellarWind, parker
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

# rc = 0.85 * RJ
# rhoc = 1800  # kg/m3

# --------------------------------------------------------- #
# -------------- Configuration settings input ------------- #

class Config:
    def __init__(self, config_file : str = "parametres.csv"):
        config = pd.read_csv(config_file, delimiter=';')
        self.database = [config.setting[i] for i in range(2) if config.value[i] == 1]
        self.magnetic_moment_models = [config.setting[i] for i in range(2,9) if config.value[i] == 1]
        self.dynamo_density_models = [config.setting[i] for i in range(9, 11) if config.value[i] == 1]
        self.planet_radius_models = [config.setting[i] for i in range(12, 13) if config.value[i] == 1]
        self.star_radius_models = [config.setting[i] for i in range(13, 14) if config.value[i] == 1]
        self.planet_luminosity_models = [
            config.setting[i] for i in range(14, 17) if config.value[i] == 1
        ]
        self.rho_crit = config.value[11]
        if config.value[5] == 1 :
            self.rc_dyn = True
        else :
            self.rc_dyn = False
        if config.value[17] == 1:
            self.talk = True
        else:
            self.talk = False

    ############ Methods ############

    def param_names(self,data : pd.DataFrame ) :
        if len(self.database) > 1 :
            raise ValueError('Only one database can be used, two were given in the configuration file')
        for database in self.database :
            if database == 'nasa_data' :
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
            elif database == 'exoplanet_data':
                new_names={'# name' : 'pl_name'}
                data = data.rename(columns=new_names)
        return data


# --------------------------------------------------------- #
# ----------------------- Prediction ---------------------- #

class Prediction:

    def __init__(self):
        pass
