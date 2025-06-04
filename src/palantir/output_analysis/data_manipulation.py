#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  15 10:09:06 2024

@author: Emilie Mauduit
"""

# --------------------------------------------------------- #
# ------------------------ Imports ------------------------ #

import pandas as pd
import numpy as np
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
        self.instrument_name = instrument_name
        self._instrument_coordinates = None
        self._instrument_sensitivity = None

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    #### Setter ####

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
                                'LOFAR low' : ['goldenrod','solid'],
                                'LOFAR high' : ['goldenrod','solid'],
                                'SKA1 low' : ['purple','solid'],
                                'SKA2 low' : ['purple', 'dashed'],
                                'GMRT' : ['y','solid'],
                                'VLA' : ['darkblue','solid'],
                                'UTR-2' : ['forestgreen','solid']}
            else :
                raise ValueError("Invalid instrument name. Available instruments are : 'NenuFAR', 'LOFAR','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'.")

        if interaction == 'MS' :
            flux_to_plot = np.array(self.data_base['pow_received_magnetic'][1:],dtype='float')
            frequencies_to_plot = np.array(self.data_base['freq_max_planet'][1:],dtype='float')
            xlabel = '$f_{c,p}^{max}$  [MHz]' ;  ylabel = '$\Phi_{radio}^{MS}$  [mJy]'

        elif interaction == 'SPI' :
            flux_to_plot = np.array(self.data_base['pow_received_spi'][1:],dtype='float')
            frequencies_to_plot = np.array(self.data_base['freq_max_star'][1:],dtype='float')
            xlabel = '$f_{c,*}^{max}$  [MHz]' ;  ylabel = '$\Phi_{radio}^{SPI}$  [mJy]'
        else :
            raise ValueError("Invalid interaction name. Available names are 'MS' or 'SPI'.")

        xmin = kwargs.get('xmin',0.9*np.min(frequencies_to_plot[~np.isnan(flux_to_plot)])) ; xmax = kwargs.get('xmax',1.1*np.max(frequencies_to_plot[~np.isnan(flux_to_plot)]))
        ymin = kwargs.get('ymin',0.9*np.min(flux_to_plot[~np.isnan(flux_to_plot)])) ; ymax = kwargs.get('ymax',1.1*np.max(flux_to_plot[~np.isnan(flux_to_plot)]))

        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        ax.scatter(frequencies_to_plot[~np.isnan(flux_to_plot)],flux_to_plot[~np.isnan(flux_to_plot)], 
            marker='+', 
            color='tab:blue')
        
        if instruments is not None :
            for name,sensitivity in sensitivity_dict.items():
                if name == 'GMRT':
                    ax.scatter(sensitivity[0,:],sensitivity[1,:], color=color_dict[name][0], label = name)
                else :
                    ax.plot(sensitivity[0,:],sensitivity[1,:], color=color_dict[name][0], linestyle=color_dict[name][1], linewidth=2, label = name)
            
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

        figname = kwargs.get("figname","")
        if figname != "" :
            plt.savefig(
                figname,
                transparent=kwargs.get('transparent',True),
                bbox_inches='tight', 
                dpi=kwargs.get('dpi',150)
                )

        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()
        plt.close('all')

    def plot_quantities(self,
        x : str,
        y : str,
        z : str = None,
        **kwargs
        ):

        """
        Function to plot any quantity with respect to any other one.
        """

        if (x not in self.data_base.keys()) or (y not in self.data_base.keys()) :
            raise ValueError("Invalid field name. Make sure you chose keys that exits in the database. Please refer to the documentation to see available field names.")
        
        xdata = np.array(self.data_base[x][1:],dtype='float')
        ydata = np.array(self.data_base[y][1:],dtype='float')
        zdata = np.array(self.data_base[z][1:],dtype='float') if z is not None else None

        xmin = kwargs.get('xmin',0.9*np.nanmin(xdata)) ; xmax = kwargs.get('xmax',1.1*np.nanmax(xdata))
        ymin = kwargs.get('ymin',0.9*np.nanmin(ydata)) ; ymax = kwargs.get('ymax',1.1*np.nanmax(ydata))

        xlabel = kwargs.get('xlabel',x + " [" + self.data_base[x][0] + "]")
        ylabel = kwargs.get('ylabel',y + " [" + self.data_base[y][0] + "]")
        zlabel = kwargs.get('zlabel',z + " [" + self.data_base[z][0] + "]") if z is not None else None

        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        if z is None :
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

    ##### TARGET SELECTION ####

    def target_selection(self,
        fc_min_MHz : float = None,
        flux_min_mJy : float = None,
        filename : str =  None):

        """ This method allows to select the targets with a maximum cyclotron frequency and a maximum estimated flux that are 
        above the sensitivity of the considered instrument. To be added : possibility to also select the one observable with 
        the telescope using its coordinates and the RA/DEC coordinates of the targets."""

        data_base_filtered = self.data_base.iloc[1:].copy()
        if (fc_min_MHz is None) and (flux_min_mJy is None) :
            sensitivities = self.instrument_sensitivity

            fc_min_MHz = np.min(sensitivities[0,:])
            flux_min_mJy = np.min(sensitivities[1,:])

        elif (fc_min_MHz is not None) and (flux_min_mJy is None) :
            flux_min_mJy = 0.
        
        elif (fc_min_MHz is None) and (flux_min_mJy is not None) :
            fc_min_MHz = 0.

        mask_MS = (data_base_filtered['freq_max_planet'].astype(float) >= fc_min_MHz) & (data_base_filtered['pow_received_magnetic'].astype(float) >= flux_min_mJy)
        mask_SPI = (data_base_filtered['freq_max_star'].astype(float) >= fc_min_MHz) & (data_base_filtered['pow_received_spi'].astype(float) >= flux_min_mJy)
        df_select = data_base_filtered[mask_MS | mask_SPI]

        if filename is not None :
            df_select.to_csv(filename, sep=';')

        return DataManipulation(df_select, instrument_name=self.instrument_name)








