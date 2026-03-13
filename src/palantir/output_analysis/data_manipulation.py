#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  15 10:09:06 2024

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

from cgi import test
import pandas as pd
import numpy as np
from typing import Tuple
from nenupy.instru.interferometer import ObservingMode
from nenupy.instru import NenuFAR
import astropy.units as u
from typing import List
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import logging
log = logging.getLogger('palantir.output_analysis.data_manipulation')


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

# ============================================================= #
# ---------------------- DataManipulation --------------------- #
# ============================================================= #

class DataManipulation:

    def __init__(self, data_base : pd.DataFrame, instrument_name : str) -> None:
        """
        Class to deal with output of Palantir.
        :param data_base:
        The data base to use, designed for the specific output of Palantir, or a database with the same structure.
        :type data_base:
        pd.DataFrame
        
        :param instrument_coordinates:
        Coordinates of the instrument selected [latitude, longitude, altitude] en deg, deg, meters
        :type instrument_coordinates:
        List
        
        :param instrument_sensitivity:
        Sensitivity (mJy) with respect to frequency (MHz) of the selected instrument
        :type intrument_sensitivity:
        np.ndarray
        """
        self.data_base = data_base
        self._dict_axis_label = None
        self.instrument_name = instrument_name
        self._instrument_coordinates = None
        self._instrument_sensitivity = None

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    #### Setter ####
    @property
    def dict_axis_label(self):
        if self._dict_axis_label is None :
            dict_axis_label = {
                'name' : ['name','[]'],
                'ra' : ['Right Ascension', '[]'], 
                'dec' : ['Declination', '[]'], 
                'planet_mass' : ['$M_p$' , '[$M_J$]'], 
                'planet_radius' : ['$R_p$', '[$R_J$]'],
                'planet_luminosity' : ['$L_p$', '[$L_S$]'],
                'star_planet_distance' : ['$d_{*-p}$', '[AU]'],
                'semi_major_axis' : ['a', '[AU]'],
                'planet_rotation_period' : ['$T_{rot,p}$', '[hr]'],
                'planet_orbital_period' : ['$T_{orb,p}$', '[days]'], 
                'star_simbad_id' : ['Simbad ID', '[]'],
                'star_mass' : ['$M_*$', '[$M_S$]'], 
                'star_radius' : ['$R_*$', '[$R_S$]'],
                'star_age' : ['t_*', '[yr]'], 
                'earth_distance' : ['s', '[pc]'],
                'star_magfield' : ['$B_*$', '[T]'],
                'star_rotperiod' : ['$\omega_{rot,*}$', '[days]'],
                'star_luminosity' : ['$L_*$', '[$L_S$]'],
                'spectral_type' : ['Spectral type', '[]'],
                'spectral_type_code' : ['Spectral type code', '[]'],
                'star_effective_temp' : ['$T_{eff,*}$', '[K]'], 
                'dynamo_density' : [r"$\rho_{dyn}$", r"[$\rho_{dyn,J}$]"],
                'dynamo_radius' : ['$r_{c,dyn}$', '[$r_{c,dyn,J}$]'],
                'B_dyn' : ['$B_{p,dyn}$', '[T]'],
                'B_eq' : ['$B_{p,eq}$', '[T]'], 
                'magnetic_moment' : ['$M_{mag,p}$', '[$M_{mag,J}$]'],
                'magnetosphere_radius' : ['$R_m$', '[$R_p$]'],
                'sw_density' : ['$n_{e,SW}$', '[$m^{-3}$]'],
                'sw_effective_velocity' : ['$v_{eff,SW}$', '[$m.s^{-1}$]'],
                'sw_velocity' : ['$v_{SW}$', '[$m.s^{-1}$]'],
                'coronal_temperature' : ['$T_{corona}$', '[K]'],
                'sw_radial_magfield_planet' : ['$B_{radial}(d_{*-p})$', '[T]'],
                'sw_azimuthal_magfield_planet' : ['$B_{\phi}(d_{*-p})$', '[T]'],
                'sw_total_magfield_planet' : ['$B_{total}(d_{*-p})$', '[T]'],
                'sw_perp_magfield_planet' : ['$B_{SW}$', '[T]'],
                'distance_alfven_point' : ['$d_{A}$', '[AU]'], 
                'alfven_velocity' : ['$v_A$', '[$m.s^{-1}$]'],
                'magnetic_field_planet' : ['$B_p$', '[T]'],
                'fc_max_planet' : ['$f_{c,p}^{max}$', '[MHz]'],
                'fp_planet' : ['$f_{p,p}$', '[MHz]'], 
                'pow_emission_kinetic' : ['$P_{kin}^{em}$', '[W]'],
                'pow_emission_magnetic' : ['$P_{mag}^{em}$', '[W]'], 
                'pow_emission_spi' : ['$P_{SPI}^{em}$', '[W]'],
                'flux_kinetic_au' : ['$\Phi_{kin}^{AU}$', '[Jy]'],
                'flux_magnetic_au' : ['$\Phi_{mag}^{AU}$', '[Jy]'],
                'flux_spi_au' : ['$\Phi_{SPI}^{AU}$', '[Jy]'], 
                'flux_received_kinetic' : ['$\Phi_{kin}$', '[mJy]'],
                'flux_received_magnetic' : ['$\Phi_{mag}$', '[mJy]'], 
                'flux_received_spi' : ['$\Phi_{SPI}$', '[mJy]'], 
                'fc_max_star' : ['$f_{c,*}^{max}$', '[MHz]'],
                'fp_star' : ['$f_{p,*}$', '[MHz]']
            }
            self._dict_axis_label = dict_axis_label
        return self._dict_axis_label

    @property
    def instrument_coordinates(self):
        if self._instrument_coordinates is None : 
            coordinates_dict = {'NenuFAR' : [47.380564,2.193236,140],
                    'LOFAR' : [53.095278,-7.915,50],
                    'SKA1 low' : [np.nan,np.nan,np.nan],
                    'SKA2 low' : [np.nan,np.nan,np.nan],
                    'GMRT' : [np.nan,np.nan,np.nan],
                    'VLA' : [np.nan,np.nan,np.nan],
                    'UTR-2' : [49.6357875,36.9395702,1460]}
            if self.instrument_name not in coordinates_dict.keys():
                raise ValueError("Invalid instrument name. Available instruments are : 'NenuFAR', 'LOFAR','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'.")
            self._instrument_coordinates = coordinates_dict[self.instrument_name]

        return self._instrument_coordinates

    @property
    def instrument_sensitivity(self):
        if self._instrument_sensitivity is None :
            NenuFAR_sensitivity = []
            frequencies = np.linspace(10,80,100)
            nenufar = NenuFAR()
            NenuFAR_sensitivity = nenufar.sensitivity( 
                frequency=frequencies*u.MHz,
                mode=ObservingMode.BEAMFORMING,
                dt=10*60*u.s,
                df=2*u.MHz,
                elevation=90*u.deg,
                efficiency=1.,
                decoherence=1.,
                lna_filter=0)

            sensitivity_dict={'NenuFAR' : np.array([frequencies,NenuFAR_sensitivity*1e3]),
                            'LOFAR low' : np.array([[14,19,30,60,75],[2.5e2,1e2,4e1,2e1,5e1]]),
                            'LOFAR high' : np.array([[130,160,190,200,250],[1,9e-1,1,2,6]]),
                            'SKA1 low' : np.array([[50,90,200,350],[5e-1,3e-1,1.5e-1,6e-2]]),
                            'SKA2 low' : np.array([[50,90,200,350],[5e-2,3e-2,1.5e-2,6e-3]]),
                            'GMRT' : np.array([[150],[3]]),
                            'VLA' : np.array([[70,80],[9e1,9e1]]),
                            'UTR-2' : np.array([[10,30],[1e1,1e1]])}
            
            if self.instrument_name not in sensitivity_dict.keys():
                raise ValueError("Invalid instrument name. Available instruments are : 'NenuFAR', 'LOFAR','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'.")
        
            self._instrument_sensitivity = sensitivity_dict[self.instrument_name]
        return self._instrument_sensitivity

    ##### Class methods ####

    @classmethod
    def from_file(cls, filename : str, instrument_name :  str):
        data_base = pd.read_csv(filename, delimiter=';')
        return cls(data_base,instrument_name)
        
    ##### PLOTS ####

    def plot_frequency_flux(self,
            interaction : str = 'MS',
            instruments : List[str] = None,
            test_alfven_velocity : bool = False,
            test_escaping : bool = False,
            **kwargs) :
        """ 
        This method allows to produce plots of predicted flux vs maximum cyclotron frequency. 
        
        :param interaction:
            As different two types of interactions are considered, it is possible to plot the results 
            concerning Star-Planet Interactions ('SPI') or Magnetosphere-Stellar wind interactions ('MS').
            Default is 'MS'.
        :type interaction:
            str

        :param instruments:
            It is possible to overplot the sensitivity of one or more radiotelescopes. In this version, 
            available ones are : 'NenuFAR', 'LOFAR low', 'LOFAR high','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'. One ore 
            more can be chosen. By default none will be overplotted.
        :type instruments:
            List[str]

        :param test_alfven_velocity:
            A parameter to distinguish, or not, points where the velocity of the emission is superior or inferior to the alfven velocity.
        :type test_alfven_velocity:
            bool
        """

        if (instruments is not None) :
            test_names = np.where(np.isin(instruments,['NenuFAR', 'LOFAR low', 'LOFAR high', 'SKA1 low','SKA2 low','GMRT','VLA','UTR-2']),True,False)
            if np.all(test_names):
                NenuFAR_sensitivity = []
                frequencies = np.linspace(10,80,100)
                nenufar = NenuFAR()
                NenuFAR_sensitivity = nenufar.sensitivity( 
                    frequency=frequencies*u.MHz,
                    mode=ObservingMode.BEAMFORMING,
                    dt=10*60*u.s,
                    df=2*u.MHz,
                    elevation=90*u.deg,
                    efficiency=1.,
                    decoherence=1.,
                    lna_filter=0)

                sensitivity_dict={'NenuFAR' : np.array([frequencies,NenuFAR_sensitivity*1e3]),
                                'LOFAR low' : np.array([[14,19,30,60,75],[2.5e2,1e2,4e1,2e1,5e1]]),
                                'LOFAR high' : np.array([[130,160,190,200,250],[1,9e-1,1,2,6]]),
                                'SKA1 low' : np.array([[50,90,200,350],[5e-1,3e-1,1.5e-1,6e-2]]),
                                'SKA2 low' : np.array([[50,90,200,350],[5e-2,3e-2,1.5e-2,6e-3]]),
                                'GMRT' : np.array([[150],[3]]),
                                'VLA' : np.array([[70,80],[9e1,9e1]]),
                                'UTR-2' : np.array([[10,30],[1e1,1e1]])}

                color_dict = {'NenuFAR' : ['tab:red','solid'],
                                'LOFAR low' : ['tab:green','dashed'],
                                'LOFAR high' : ['tab:green','solid'],
                                'SKA1 low' : ['tab:orange','solid'],
                                'SKA2 low' : ['tab:orange', 'dashed'],
                                'GMRT' : ['tab:red','solid'],
                                'VLA' : ['tab:blue','solid'],
                                'UTR-2' : ['tab:blue','dashed']}
            else :
                raise ValueError("Invalid instrument name. Available instruments are : 'NenuFAR', 'LOFAR','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'.")

        if interaction == 'MS' :
            flux_to_plot = np.array(self.data_base['flux_received_magnetic'][1:],dtype='float')
            frequencies_to_plot = np.array(self.data_base['fc_max_planet'][1:],dtype='float')
            xlabel = '$f_{c,p}^{max}$  [MHz]' ;  ylabel = '$\Phi_{radio}^{MS}$  [mJy]'

            if test_escaping :
                fp_planet = np.array(self.data_base['fp_planet'][1:][~np.isnan(flux_to_plot)],dtype='float')
                fc_planet = frequencies_to_plot[~np.isnan(flux_to_plot)]

                escaping = fc_planet > fp_planet
                not_escaping = ~escaping

        elif interaction == 'SPI' :
            flux_to_plot = np.array(self.data_base['flux_received_spi'][1:],dtype='float')
            frequencies_to_plot = np.array(self.data_base['fc_max_star'][1:],dtype='float')
            xlabel = '$f_{c,*}^{max}$  [MHz]' ;  ylabel = '$\Phi_{radio}^{SPI}$  [mJy]'

            if test_escaping :
                fp_star = np.array(self.data_base['fp_star'][1:][~np.isnan(flux_to_plot)],dtype='float')
                fc_star = frequencies_to_plot[~np.isnan(flux_to_plot)]

                escaping = fc_star > 10*fp_star
                not_escaping = ~escaping
        else :
            raise ValueError("Invalid interaction name. Available names are 'MS' or 'SPI'.")

        if test_alfven_velocity :
            alfven_velocity = np.array(self.data_base['alfven_velocity'][1:][~np.isnan(flux_to_plot)],dtype='float')
            sw_velocity = np.array(self.data_base['sw_velocity'][1:][~np.isnan(flux_to_plot)],dtype='float')
            sw_velocity_eff = np.array(self.data_base['sw_effective_velocity'][1:][~np.isnan(flux_to_plot)],dtype='float')

            sub_alfvenic_sw = alfven_velocity > sw_velocity
            sub_alfvenic_eff = alfven_velocity > sw_velocity_eff
            super_alfvenic_eff = ~sub_alfvenic_eff
            super_alfvenic_sw = ~sub_alfvenic_sw

        frequencies_to_plot = frequencies_to_plot[~np.isnan(flux_to_plot)]
        flux_to_plot = flux_to_plot[~np.isnan(flux_to_plot)]

        xmin = kwargs.get('xmin',0.8*np.min(frequencies_to_plot)) ; xmax = kwargs.get('xmax',1.2*np.max(frequencies_to_plot))
        ymin = kwargs.get('ymin',0.8*np.min(flux_to_plot)) ; ymax = kwargs.get('ymax',1.2*np.max(flux_to_plot))
        print(len(frequencies_to_plot))

        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        if test_alfven_velocity :
            ax.scatter(frequencies_to_plot[sub_alfvenic_eff & sub_alfvenic_sw], flux_to_plot[sub_alfvenic_eff & sub_alfvenic_sw], marker = 'v', color='tab:orange', label='$v_{SW}$ and $v_{SW,eff}$ < $v_A$', alpha=0.6)
            ax.scatter(frequencies_to_plot[super_alfvenic_eff & sub_alfvenic_sw], flux_to_plot[super_alfvenic_eff & sub_alfvenic_sw], marker = 'd', color='tab:green', label='$v_{SW}$ < $v_A$ and $v_{SW,eff}$ > $v_A$', alpha=0.6)
            ax.scatter(frequencies_to_plot[super_alfvenic_eff & super_alfvenic_sw], flux_to_plot[super_alfvenic_eff & super_alfvenic_sw], marker ='^', color='tab:purple', label='$v_{SW}$ and $v_{SW,eff}$ > $v_A$', alpha=0.6)

        elif test_escaping:
            ax.scatter(frequencies_to_plot[not_escaping], flux_to_plot[not_escaping], marker = '^', color='tab:purple', label="$f_{ce}$ < $f_{pe}$" if interaction=="MS" else "$f_{ce}$ < 10$f_{pe}$", alpha=0.6)
            ax.scatter(frequencies_to_plot[escaping], flux_to_plot[escaping], marker ='v', color='tab:orange', label="$f_{ce}$ > $f_{pe}$" if interaction == "MS" else "$f_{ce}$ > 10$f_{pe}$", alpha=0.6)
        else :
            ax.scatter(frequencies_to_plot,flux_to_plot, 
                marker='^', 
                color='gray',
                alpha=0.8)
            
        if instruments is not None :
            for name,sensitivity in sensitivity_dict.items():
                if name == 'GMRT':
                    ax.scatter(sensitivity[0,:],sensitivity[1,:], color=color_dict[name][0], label = name)
                else :
                    ax.plot(sensitivity[0,:],sensitivity[1,:], color=color_dict[name][0], linestyle=color_dict[name][1], linewidth=3, label = name)
            
        ax.plot([10.0,10.0],[ymin,ymax], linestyle = 'dashed', color='black',label='ionospheric cut-off')
        rect = plt.Rectangle((xmin,ymin),10-xmin,ymax-ymin,facecolor='black',alpha=0.1)
        ax.add_patch(rect)
        ax.set_xlabel(xlabel, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        ax.tick_params(axis='both',labelsize=14)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(kwargs.get('title',""))
        plt.grid()
        plt.legend(fontsize=12, ncol = 1)
        plt.tight_layout()

        figname = kwargs.get("figname","")
        if figname != "" :
            plt.savefig(
                figname,
                transparent=kwargs.get('transparent',True),
                bbox_inches='tight', 
                dpi=kwargs.get('dpi',150)
                )
        
        plt.show()
        plt.close('all')

    def plot_quantities(self,
        x : str,
        y : str,
        z : str = None,
        test_alfven_velocity : bool = False,
        **kwargs
        ):

        """
        Function to plot any quantity with respect to any other one.
        """

        if (x not in self.data_base.keys()) or (y not in self.data_base.keys()) :
            raise ValueError("Invalid field name. Make sure you chose keys that exits in the database. Field names can be : {}".format(self.data_base.keys()))
        
        xdata = np.array(self.data_base[x][1:],dtype='float')
        ydata = np.array(self.data_base[y][1:],dtype='float')
        zdata = np.array(self.data_base[z][1:],dtype='float') if z is not None else None

        if test_alfven_velocity :
            alfven_velocity = np.array(self.data_base['alfven_velocity'][1:][~np.isnan(ydata)], dtype='float')
            sw_velocity = np.array(self.data_base['sw_velocity'][1:][~np.isnan(ydata)], dtype='float')
            sw_velocity_eff = np.array(self.data_base['sw_effective_velocity'][1:][~np.isnan(ydata)],dtype='float')

            sub_alfvenic_sw = alfven_velocity > sw_velocity
            sub_alfvenic_eff = alfven_velocity > sw_velocity_eff
            super_alfvenic_eff = ~sub_alfvenic_eff
            super_alfvenic_sw = ~sub_alfvenic_sw

        xmin = kwargs.get('xmin',0.9*np.nanmin(xdata)) ; xmax = kwargs.get('xmax',1.1*np.nanmax(xdata))
        ymin = kwargs.get('ymin',0.9*np.nanmin(ydata)) ; ymax = kwargs.get('ymax',1.1*np.nanmax(ydata))

        xlabel = kwargs.get('xlabel', self.dict_axis_label[x][0] + " " + self.dict_axis_label[x][1])
        ylabel = kwargs.get('ylabel', self.dict_axis_label[y][0] + " " + self.dict_axis_label[y][1])
        zlabel = kwargs.get('zlabel', self.dict_axis_label[z][0] + " " + self.dict_axis_label[z][1]) if z is not None else None

        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        if z is None :
            xdata = xdata[~np.isnan(ydata)] ; ydata = ydata[~np.isnan(ydata)]
            if test_alfven_velocity :
                ax.scatter(xdata[sub_alfvenic_eff & sub_alfvenic_sw], ydata[sub_alfvenic_eff & sub_alfvenic_sw], marker = 'v', color='tab:orange', label='$v_{SW}$ and $v_{SW,eff}$ < $v_A$', alpha=0.6)
                ax.scatter(xdata[super_alfvenic_eff & sub_alfvenic_sw], ydata[super_alfvenic_eff & sub_alfvenic_sw], marker = 'd', color='tab:green', label='$v_{SW}$ < $v_A$ and $v_{SW,eff}$ > $v_A$', alpha=0.6)
                ax.scatter(xdata[super_alfvenic_eff & super_alfvenic_sw], ydata[super_alfvenic_eff & super_alfvenic_sw], marker ='^', color='tab:purple', label='$v_{SW}$ and $v_{SW,eff}$ > $v_A$', alpha=0.6)
                ax.legend(fontsize=12)
            else :
                ax.scatter(xdata[~np.isnan(ydata)],ydata[~np.isnan(ydata)], 
                    marker='+', 
                    color='tab:blue')
        else :
            zscale = kwargs.get("zscale","linear")
            zdata = 10*np.log10(zdata) if zscale == "log" else zdata
            im = ax.scatter(xdata[~np.isnan(zdata)],ydata[~np.isnan(zdata)], 
                c=zdata[~np.isnan(zdata)], marker='o', cmap = 'viridis',
                vmin = kwargs.get("zmin", 0.9 * np.min(zdata[~np.isnan(zdata)])),
                vmax = kwargs.get("zmax", 1.1 * np.max(zdata[~np.isnan(zdata)])))
            cax = inset_axes(
                    ax,
                    width='3%',
                    height='100%',
                    loc='lower left',
                    bbox_to_anchor=(1.03, 0., 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                )
            cbar = plt.colorbar(im, cax=cax,)
            cbar.set_label(zlabel)
        
        ax.set_xlabel(xlabel, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        ax.tick_params(axis='both',labelsize=14)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xscale(kwargs.get('xscale','linear'))
        ax.set_yscale(kwargs.get('yscale','linear'))
        ax.set_title(kwargs.get('title',""))
        ax.grid()

        figname = kwargs.get("figname","")
        if figname != "" :
            plt.savefig(
                figname,
                transparent=kwargs.get('transparent',True),
                bbox_inches='tight', 
                dpi=kwargs.get('dpi',150)
                )

        plt.show()
        plt.close('all')

    def plot_simple_histogram_parameters(self,
        xfield : str,
        **kwargs) -> None : 
        
        """This functions allows to do a simple histogram any parameter. 

        :param xfield:
            The key of the column to use as the x axis.
        :type xfield:
            str

        .. rubric:: `**kwargs`

        :param xlabel:
            To specify what to put as the label of the x-axis. Default is `xfield`.
        :type xlabel:
            str

        :param title:
            The title to give the plot. Default is no title.
        :type title:
            str

        :param bins:
            The number of bins for the histogram. Default is 25.
        :type bins:
            int

        :param range:
            The range of values for the histogram. Default is (vmin,vmax).
        :type range:
            Tuple[float,float]

        :param color:
            The color of the points. Default is 'tab:blue'.
        :type color:
            str

        :param figname:
            If this parameter is set, the plot will be saved with the given name. Default is None.
        :type figname:
            str

        :param transparent:
            To specify if the background of the plot should be transparent or white. Default is True.
        :type transparent:
            bool

        :param dpi:
            To define the resolution of the plot if saved. Default is 150.
        :type dpi:
            int
        """
        if (xfield not in self.data_base.keys()) :
            raise ValueError("Invalid field name. Make sure you chose keys that exits in the database. Field names can be : {}".format(self.data_base.keys()))
        
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        ax.hist(np.array(self.data_base[xfield][1:],dtype='float'), 
                bins = kwargs.get('bins',25),
                color = kwargs.get('color','tab:blue'),
                range = kwargs.get('range',),
                histtype = 'step'
            )

        ax.set_xlabel(kwargs.get('xlabel', self.dict_axis_label[xfield][0] + " " + self.dict_axis_label[xfield][1]), fontsize=18)
        ax.set_ylabel('Count', fontsize=18)
        ax.tick_params(axis='both',labelsize=14)
        ax.set_title(kwargs.get('title',""))
        plt.grid()
        plt.tight_layout()

        figname = kwargs.get("figname","")
        if figname != "" :
            plt.savefig(
                figname,
                transparent=kwargs.get('transparent',True),
                bbox_inches='tight', 
                dpi=kwargs.get('dpi',150)
                )

        plt.show()
        plt.close('all')

    ##### TARGET SELECTION ####

    def target_selection(self,
        fc_range_MHz : Tuple[float,float] = None,
        flux_min_mJy : float = None,
        declination_range : Tuple[float,float] = None,
        filename : str =  None):

        """ This method allows to select the targets with a maximum cyclotron frequency and a maximum estimated flux that are 
        above the sensitivity of the considered instrument. Setting a declination range allows the possibility to also select the 
        one observable with the telescope using its coordinates and the RA/DEC coordinates of the targets."""

        data_base_filtered = self.data_base.iloc[1:].copy()
        sensitivities = self.instrument_sensitivity
        if (fc_range_MHz is None) :
            fc_range_MHz = (np.min(sensitivities[0,:]), np.max(sensitivities[0,:]))

        if (flux_min_mJy is None) :
            flux_min_mJy = np.min(sensitivities[1,:])

        mask_MS = (data_base_filtered['fc_max_planet'].astype(float) >= fc_min_MHz) & (data_base_filtered['flux_received_magnetic'].astype(float) >= flux_min_mJy)
        mask_SPI = (data_base_filtered['fc_max_star'].astype(float) >= fc_min_MHz) & (data_base_filtered['flux_received_spi'].astype(float) >= flux_min_mJy)
        if declination_range is not None :
            mask_observable = (data_base_filtered['dec'].astype(float) >= declination_range[0]) & (data_base_filtered['dec'].astype(float) <= declination_range[1])
        else :
            mask_observable = np.ones_like(mask_MS).astype(bool)
        df_select = data_base_filtered[(mask_MS | mask_SPI) & mask_observable]

        if filename is not None :
            df_select.to_csv(filename, sep=';')

        return DataManipulation(df_select, instrument_name=self.instrument_name)








