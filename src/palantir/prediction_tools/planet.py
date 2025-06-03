#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:53:23 2021

@author: Emilie Mauduit
"""


import pandas as pd
import numpy as np
from math import pow
from typing import List
from scipy.interpolate import interp1d
from importlib_resources import files

import logging
log = logging.getLogger('palantir.prediction_tools.planet')
# ============================================================= #
# --------------------------- Planet -------------------------- #
# ============================================================= #


class Planet:
    def __init__(
        self,
        name: str,
        mass: float,
        radius: dict,
        distance: float,
        eccentricity: float,
        worb: dict,
        luminosity: dict,
        detection_method: str = None,
        wrot: float = None,
    ):
        """Define a planet object. Every parameter must be normalized at Jupiter.
        :param name:
            Name of the planet
        :type name:
            str
        :param mass:
            Mass of the planet, expected to be normalized to Jupiter's.
        :type mass:
            float
        :param radius:
            Radius of the planet, expected to be normalized to Jupiter's.
        :type radius:
            float
        :param distance:
            Distance between the planet and its host star, in AU.
        :type distance:
            float
        :param eccentricity:
            Eccentricity of the orbit, adimensional. Should be between 0 and 1.
        :type eccentricity:
            float
        :param worb:
            Orbital period of the planet, expected to be normalized to Jupiter's.
        :type worb:
            float
        :param luminosity:
            Equivalent luminosity of the planet, expected to be normalized to the Sun's.
        :type luminosity:
            float
        :param detection_method:
            The method used to detect this planet.
        :type detection_method:
            str
        :param wrot:
            Rotation rate of the planet, expected to be normalized to Jupiter's.
        :type wrot:
            float
        """

        self.name = name
        self.mass = mass
        self.radius = radius
        self.stardist = distance
        self.eccentricity = eccentricity
        self.orbitperiod = worb
        self.rotrate = wrot
        self.detection_method = detection_method
        self.luminosity = luminosity
    
    def __str__(self):
        return("Name : " + self.name + "\n"
            + "Mass : Mp ={} MJ\n".format(self.mass)
            + "Radius : Rp ={} RJ\n".format(self.radius)
            + "Rotation rate : wrot = {} wrotJ\n".format(self.rotrate)
            + "Orbital period : worb = {} worbJ\n".format(self.orbitperiod)
            + "Luminosity : Lp = {} LS\n".format(self.luminosity)
            + "Distance to host star : d_sp = {} AU".format(self.stardist)
        )

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value: dict):
        models = value["models"]
        radius = value["radius"]
        if np.isnan(radius):
            self._radius = self._calculate_radius(models, self.mass)
        else:
            self._radius = radius

    @property
    def orbitperiod(self):
        return self._orbitperiod

    @orbitperiod.setter
    def orbitperiod(self, value: dict):
        star_mass = value["star_mass"]
        orbitperiod = value["worb"]
        dua = 1.49597870700e11  # m
        wJ = 1.77e-4  # s-1
        if np.isnan(orbitperiod):
            d = self.stardist * dua
            G = 6.6725985e-11
            MS = 1.989e30
            self._orbitperiod = pow(star_mass * MS * G / pow(d, 3), 1 / 2.0) / wJ
        else:
            self._orbitperiod = orbitperiod

    @property
    def luminosity(self):
        return self._luminosity

    @luminosity.setter
    def luminosity(self, value: dict):
        models = value["models"]
        luminosity = value["luminosity"]
        star_age = value["star_age"]

        if np.isnan(luminosity):
            self._luminosity = self._calculate_luminosity(
                models, planet_mass=self.mass, star_age=star_age
            )
        else:
            self._luminosity = luminosity

    def unnormalize_mass(self):
        mass_jup = 1.8986e27  # kg
        mass = self.mass * mass_jup
        return mass

    def unnormalize_radius(self):
        radius_jup = 71492e3  # m
        radius = self.radius * radius_jup
        return radius

    def unnormalize_orbital_period(self):
        worb_jup = 1.77e-4  # s-1
        worb = self.orbitperiod * worb_jup
        return worb

    def unnormalize_rotrate(self):
        wrot_jup = 1.77e-4  # s-1
        wrot = self.rotrate * wrot_jup
        return wrot

    def radius_expansion(self,luminosity: float) -> None :
        """
        Compute the factor of expansion of the planetary radius, depending on the equilibrium temperature of the host star.
        Uses a simple polyfit model of R_measure/R_predicted with respect to equilibrium temperature, (in log-lin scale).
        And a linear correction from R_measure vs R_predicted*expension factor.
        :param luminosity:
            Luminosity of the star, either known or computed from (Tout et al, 1996)
        :type luminosity:
            float
        :param eccentricity:
            Eccentricity of the planetary orbit
        :type eccentricity:
            float
        """

        LS = 3.826e26  # W
        A = 0.4
        sigmaSB = 5.670374419e-8

        
        B = (1 + (self.eccentricity**2) / 2) ** 2 if not np.isnan(self.eccentricity) else 1.
        
        d = self.stardist * 1.49597870700e11
        Teq = pow(
            (1 - A) * luminosity * LS / (16 * np.pi * (d**2) * sigmaSB * B), 1.0 / 4
        )

        ##expansion factor depending on Teq
        coeffs_teq = [7.13173118 ,-4.78521151 , 0.93174824]
        expansion_factor = 0
        for i,c in enumerate(coeffs_teq):
            expansion_factor += c*np.power(np.log10(Teq),i)
        
        Rp = self._radius

        ##linear correction
        coeffs_lin = [-0.06070923,  1.39140901, -1.19830686, -0.90527977,  0.655981  ]
        lin_corr = 0
        for i,c in enumerate(coeffs_lin):
            lin_corr += c*np.power(np.log10(Rp * expansion_factor),i)
        
        self._radius =  (10**lin_corr)

    def radius_expansion_jmg2007(self, luminosity: float):
        """Compute the factor of expansion of the planetary radius, depending on the distance to the host star. 
        Uses Griessmeier et al, 2007 formula.
        """

        LS = 3.826e26  # W
        ct1 = 764
        ct2 = 0.28
        cg1 = 0.59
        cg2 = 1.03
        A = 0.4
        sigmaSB = 5.670374419e-8

        B = (1 + (self.eccentricity**2) / 2) ** 2
        d = self.stardist * 1.49597870700e11

        T0 = ct1 * pow(self.mass, ct2)
        gamma = 1.15 + 0.05 * pow(cg1 / self.mass, cg2)
        Teq = pow(
            (1 - A) * luminosity * LS / (16 * np.pi * (d**2) * sigmaSB * B), 1.0 / 4
        )
        Rp = self._radius
        self._radius = Rp * (1 + 0.05 * pow(Teq / T0, gamma))
        return Teq

    def tidal_locking(self, age: float, star_mass: float, Qpp: float = 3.16e5):
        """Computes the rotation rate of the planet depending on a synchronized or free rotation.
        Depending on the age of the host star, a limit distance is evaluated and if the star-planet distance is lower than that we consider a synchronized orbit, else we assume the rotation rate to be Jupiter's.

            :param age:
                Age of the host star.
            :type age:
                float
            :param star_mass:
                Mass of the star, normalized to the Sun's.
            :type star_mass:
                float
        """

        wrot_J = 1.77e-4  # s-1
        d = self.stardist
        dsync = self._calculate_synchro_dist(
            self.rotrate * wrot_J, Qpp, age, self.mass, self.radius, star_mass
        )
        if d <= dsync:
            self.rotrate = self.orbitperiod
        else:
            self.rotrate = self.rotrate

    @staticmethod
    def _calculate_radius(
        models: List[str], mass: float, Rmean: bool = True, Rmax: bool = False
    ):
        R = []
        for model in models:
            if "radius_original" in model:
                Mmax = 3.16 * 1.8986e27  # kg
                rho0 = 394  # kg.m-3
                RJ = 69911e3
                R.append(
                    pow((4 / 3) * np.pi * rho0, -1.0 / 3)
                    * pow(mass * 1.8986e27, 1.0 / 3)
                    / (1 + pow(mass * 1.8986e27 / Mmax, 2.0 / 3))
                    / RJ
                )
            if "polyfit" in model :
                coeffs = [0.08960267, -0.03390545, -0.32544411,  0.35903958, 0.19102357, -0.1489227,
                    -0.08079977,  0.01321372,  0.01195444,  0.00157772]
                Rtemp = 0
                for i,c in enumerate(coeffs):
                    Rtemp += c*np.power(np.log10(mass),i)
                R.append(10**Rtemp)

        if Rmean:
            #print("Rmean : ", np.mean(R))
            return np.mean(R)
        if Rmax:
            return np.max(R)



    @staticmethod
    def _calculate_luminosity(
        models,
        planet_mass: float,
        star_age: float,
        Lmean: bool = True,
        Lmax: bool = False,
    ):

        """Retourne la liste des luminosités en fonction de l'âge selon la valeur de M"""

        MJ = 1.8986e27
        ME = 6e24
        M = planet_mass
        L = []
        maps_dir = files("palantir.scripts.input_files")
        for model in models:
            if model == "Burrows":
                table = pd.read_csv(maps_dir / "Burrows.csv", delimiter=";")
                mass_dict = ["M=1MJ_", "M=5MJ_", "M=10MJ_", "M=20MJ_"]
                luminosities = [
                    np.interp(
                        np.log10(star_age),
                        table[mass + "log(t) (Gyr)"],
                        table[mass + "log(L/Ls)"],
                    )
                    for mass in mass_dict
                ]
                masses = np.log10([1.0, 5.0, 10.0, 20.0])
                L.append(10 ** (np.interp(np.log10(M), masses, luminosities)))
            elif model == "Baraffe_noirrad":
                #print('mass in Baraffe = ', M)
                """Tables taken from Baraffe et al, 2008, link :
                https://perso.ens-lyon.fr/isabelle.baraffe/PLANET08/
                Same for the other model."""
                table = pd.read_csv(maps_dir /"Baraffe_no_irrad.csv", delimiter=";")
                age_dict = [
                    "t=0.01_log(L/Ls)",
                    "t=0.05_log(L/Ls)",
                    "t=0.10_log(L/Ls)",
                    "t=0.50_log(L/Ls)",
                    "t=1.00_log(L/Ls)",
                    "t=3.00_log(L/Ls)",
                    "t=5.00_log(L/Ls)",
                    "t=7.00_log(L/Ls)",
                ]
                luminosities = []
                for age in age_dict :
                    f_noirrad = interp1d(np.log10(table["M/M_E"]),table[age],kind='linear',fill_value='extrapolate')
                    luminosities.append(f_noirrad(np.log10(M * MJ / ME)))
                #print("luminosities : ", luminosities)
                ages = np.array([0.01, 0.05, 0.10, 0.50, 1.0, 3.0, 5.0, 7.0])
                f_ages = interp1d(np.log10(ages),luminosities,kind='linear',fill_value='extrapolate')
                L.append(10 ** f_ages(np.log10(star_age)))
            elif model == "Baraffe_irrad":
                table = pd.read_csv(maps_dir /"Baraffe_irrad.csv", delimiter=";")
                age_dict = [
                    "t=0.01_log(L/Ls)",
                    "t=0.05_log(L/Ls)",
                    "t=0.10_log(L/Ls)",
                    "t=0.50_log(L/Ls)",
                    "t=1.00_log(L/Ls)",
                    "t=3.00_log(L/Ls)",
                    "t=5.00_log(L/Ls)",
                    "t=7.00_log(L/Ls)",
                ]
                luminosities = [
                    np.interp(np.log10(M * MJ / ME), np.log10(table["M/M_E"]), table[age])
                    for age in age_dict
                ]
                ages = np.array([0.01, 0.05, 0.10, 0.50, 1.0, 3.0, 5.0, 7.0])
                L.append(10 ** (np.interp(np.log10(star_age), np.log10(ages), luminosities)))

        if Lmean and not Lmax:
            return pow(np.prod(L), 1 / len(L))
        elif Lmax and not Lmean:
            return np.max(L)
        elif not Lmax and not Lmean:
            return L[0]
        else:
            log.error("Wrong value for Lmean :{Lmean} or Lmax : {Lmax}, only one can be set to True")
            raise ValueError(
                "Wrong value for Lmean :{Lmean} or Lmax : {Lmax}, only one can be set to True"
            )

    @staticmethod
    def _calculate_synchro_dist(
        wrot: float,
        Qp: float,
        tsync: float,
        planet_mass: float,
        planet_radius: float,
        star_mass: float,
    ):
        """Compute the distance at which a planet should be so that its orbit is synchronised, depending on the age of the system.
        :param wrot:
            orbital period
        :type wrot:
            float
        :param Qp:
            factor of dissipation by tidal effect
        :type Qp:
            float
        :param tsync:
            time of synchronization
        :type tsync:
            float
        :param planet_mass:
            Mass of the planet [Mjup]
        :type planet_mass:
            float
        :param planet_radius:
            Radius of the planet [Rjup]
        :type planet_radius:
            float
        :param star_mass:
            Mass of the star [Msun]
        :type star_mass:
            float
        """
        G = 6.6725985e-11
        a = 9 * tsync / (4 * wrot * 0.26 * Qp)
        b = G * planet_mass / pow(planet_radius, 3)
        c = pow(star_mass / planet_mass, 2)
        res = planet_radius * pow(a * b * c, 1 / 6)
        return res
