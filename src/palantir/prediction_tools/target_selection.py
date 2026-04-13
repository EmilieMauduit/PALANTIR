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
import numpy as np
import astropy.units as u
import pyvo
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import time


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
        self.star_radius_models = [config.setting[i] for i in range(16, 18) if int(config.value[i]) == 1]
        self.planet_luminosity_models = [
            config.setting[i] for i in range(18, 21) if int(config.value[i]) == 1
        ]
        self.star_magfield_models = [config.setting[i] for i in range(22, 25) if int(config.value[i]) == 1]
        self.star_magfield_catalog_only = True if (int(config.value[25]) == 1) else False
        self.rho_crit = int(config.value[12])
        self.rc_dyn = True if int(config.value[6]) == 1 else False
        self.radius_expansion = True if int(config.value[15]) == 1 else False
        self.sp_type_code = int(config.value[21])
        self.talk = True if int(config.value[26]) == 1 else False
        self.output_path = config.value[27]
        self.output_params = ["name","ra","dec","planet_mass", "planet_radius", "planet_luminosity", "star_planet_distance", "semi_major_axis",
            "planet_rotation_period", "planet_orbital_period","tidally_locked","star_simbad_id", "star_mass","star_radius","star_age","earth_distance",
            "star_magfield","star_rotperiod","star_luminosity","star_Xray_flux","spectral_type", "spectral_type_code","star_effective_temp","dynamo_density",
            "dynamo_radius","B_dyn" ,"B_eq","magnetic_moment","magnetosphere_radius","sw_density","sw_effective_velocity","sw_velocity",
            "coronal_temperature","sw_radial_magfield_planet","sw_azimuthal_magfield_planet","sw_total_magfield_planet","sw_perp_magfield_planet",
            "distance_alfven_point","alfven_velocity","magnetic_field_planet","fc_max_planet","fp_planet",
            "pow_emission_kinetic","pow_emission_magnetic","pow_emission_spi", "flux_kinetic_au", "flux_magnetic_au",
            "flux_spi_au", "flux_received_kinetic","flux_received_magnetic", "flux_received_spi",
            "fc_max_star", "fp_star"]
        self.output_params_units = ["","deg", "deg", "MJ", "RJ", "LS", "AU", "AU","hr","days","","", "MS", "RS", "yr", "pc", "T", "days",
            "LS", "erg.cm-2.s-1","", "","K", "rho_dyn_J", "r_dyn_J", "T", "T", "MmagJ", "Rp", "m-3", "m.s-1", "m.s-1", "K", "T","T","T","T","AU", "m.s-1", "T", "MHz", "MHz", "W", 
            "W", "W","Jy", "Jy", "Jy", "mJy","mJy", "mJy", "MHz", "MHz"]
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
# ------------------- ROSAT Catalog queries ------------------- #
# ============================================================= #

class XRayFluxCalculator:
    def __init__(self, catalog="IX/30", radius=10*u.arcsec):
        """
        Instanciate a calculator to avoid multiple/repetitive queries to the Vizier database.

        :param catalog:
            Name of the catalog to load from Vizier. Default is the ROSAT catalog for X-ray measurements: "IX/30"
        :type catalog:
            str
        :param radius:
            The radius of acceptance for a match. Default is 10 arcsec
        :type radius:
            :class::`~a`stropy.Quantity
        """
        self.catalog_name = catalog
        self.catalogue_ivoid = f"ivo://CDS.VizieR/{self.catalog_name}"
        self.voresource = pyvo.registry.search(ivoid=self.catalogue_ivoid)[0]
        self.url = self.voresource.reference_url
        self.tables = self.voresource.get_tables()
        self.radius = radius
        self.vizier_table = None
        self._load_catalog()
    
    def _load_catalog(self):
        """Load the ROSAT catalog once."""

        query_columns = f"""
                SELECT RAJ2000, DEJ2000, Crate, HR1
                FROM "{list(self.tables.keys())[0]}"
                """
        
        tap_service = self.voresource.get_service("tap")
        query_result = tap_service.search(query_columns)
        if query_result:
            self.vizier_table = query_result
            logger.info(f"ROSAT catalog loaded with {len(self.vizier_table)} stars.")
        else:
            logger.info("No data found in the ROSAT catalog.")
            self.vizier_table = None
    
    def _find_closest_star(self, target_coord, max_radius=None):
        """
        Find the closest star in the catalog from given coordinates.

        :param target_coord:
            The target coordinates.
        :type target_coord:
            :class:`~astropy.coordinates.SkyCoord`
        :param max_radius:
            The maximum radius (in arcsec), to look around the target.
        :type max_radius:
            :class:`~astropy.Quantity`, optional

        :returns:
            A tuple with the index of the match found in the catalog and the distance from the target.
        :rtype:
            Tuple[float,float]
        """
        if self.vizier_table is None or len(self.vizier_table) == 0:
            return None, None

        catalog_coords = SkyCoord(ra=self.vizier_table.getcolumn(name='RAJ2000').data, 
                                dec=self.vizier_table.getcolumn(name='DEJ2000').data, 
                                unit=(u.deg, u.deg), 
                                frame='icrs')

        separations = target_coord.separation(catalog_coords)

        search_radius = max_radius if max_radius is not None else self.radius
        mask = separations <= search_radius
        
        if not np.any(mask):
            return None, None

        closest_idx = np.argmin(separations)
        min_distance = separations[closest_idx]
        
        return closest_idx, min_distance
    
    def compute_xray_flux_from_coords(self, target_coord, star_name=None):
        """
        Computes the XRay flux of a star by using a previously loaded ROSAT catalog from Vizier, but using coordinates and not make a query to the server.

        :param target_coord:
            The coordinates of the target.
        :type target_coord:
            :class:`~astropy.coordinates.SkyCoord`
        :param star_name:
            The name of the star, only used for logging purposes. Default is None.
        :type star_name:
            str

        :returns:
            XRay flux in erg.cm-2.s-1, or np.nan if no match is found.
        :rtype:
            float
        """
        closest_idx, distance = self._find_closest_star(target_coord)
        
        if closest_idx is None:
            name_info = f" ({star_name})" if star_name else ""
            logger.debug(f"No star found for {target_coord.ra.value:.6f}, {target_coord.dec.value:.6f}{name_info}")
            return np.nan

        HR1 = self.vizier_table.getcolumn(name='HR1').data[closest_idx]
        Counts = self.vizier_table.getcolumn(name='Crate').data[closest_idx]

        ECF_rosat = (5.30 * HR1 + 8.7) * 1e-12  # erg.cm-2.counts-1
        Xray_flux = ECF_rosat * Counts  # erg.cm-2.s-1

        name_info = f" ({star_name})" if star_name else ""
        logger.debug(f"Flux X-ray found for {target_coord.ra.value:.6f}, {target_coord.dec.value:.6f}{name_info}: distance={distance.arcsec:.2f}arcsec, flux={Xray_flux:.2e} erg.cm-2.s-1")
        
        return Xray_flux
    
    def compute_xray_flux_from_name(self, star_name):
        """
        Computes the XRay flux of a star by making a query with its name to Vizier on the ROSAT catalog.

        :param star_name:
            The name of the star, preferably its main ID on SIMBAD.
        :type star_name:
            str

        :returns:
            XRay flux in erg.cm-2.s-1, or np.nan if no match is found.
        :rtype:
            float
        """
        vizier = Vizier(columns=["HR1", "Crate", "_r", "RAJ2000", "DEJ2000"])
        query_result = vizier.query_object(star_name, 
                                        radius=self.radius,
                                        catalog=self.catalog)
        
        time.sleep(1) 
        
        if not query_result:
            logger.debug(f"No match found in ROSAT catalog for star: {star_name}")
            return np.nan
        
        table = query_result
        if ('HR1' and 'Crate') in table.keys():
            imin = np.argmin(table["_r"])
            HR1 = table['HR1'][imin]
            ECF_rosat = (5.30 * HR1 + 8.7) * 1e-12
            Counts = table['Crate'][imin]
            Xray_flux = ECF_rosat * Counts
            return Xray_flux
        else:
            logger.debug(f"Columns HR1/Crate not found for {star_name}")
            return np.nan

# ============================================================= #
# ------------------------- Prediction ------------------------ #
# ============================================================= #

class Prediction:

    def __init__(self):
        pass
