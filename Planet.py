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

from calc_tools import synchro_dist


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
        worb: dict,
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
        :param wrot:
            Rotation rate of the planet, expected to be normalized to Jupiter's.
        :type wrot:
            float
        :param worb:
            Orbital period of the planet, expected to be normalized to Jupiter's.
        :type worb:
            float
        :param distance:
            Distance between the planet and its host star, in AU.
        :type distance:
            float
        :param axis:
            Semi-major axis, in AU.
        :type axis:
            float
        """

        self.name = name
        self.mass = mass
        self.radius = radius
        self.stardist = distance
        self.orbitperiod = worb
        self.rotrate = wrot
        self.detection_method = detection_method

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
        if np.isnan(orbitperiod):
            d = self.stardist * dua
            G = 6.6725985e-11
            MS = 1.989e30
            wJ = 1.77e-4  # s-1
            self._orbitperiod = pow(star_mass * MS * G / pow(d, 3), 1 / 2.0) / wJ
        else:
            self._orbitperiod = orbitperiod

    def talk(self, talk: bool):
        if talk:
            print("Name : ", self.name)
            print("Mass : ", self.mass, " Mj")
            print("Radius : ", self.radius, " Rj")
            print("Rotation rate : ", self.rotrate, " wj")
            print("Orbital period : ", self.orbitperiod, "woj")
            print("Distance to host star : ", self.stardist, " AU")

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

    def radius_expansion(self, luminosity: float, eccentricity: float):
        """Compute the factor of expansion of the planetary radius, depending on the distance to the host star.
        :param mass:
            Planetary mass, in MJ.
        :type mass:
            float
        """

        LS = 3.826e26  # W
        ct1 = 764
        ct2 = 0.28
        cg1 = 0.59
        cg2 = 1.03
        A = 0.4
        sigmaSB = 5.670374419e-8
        B = (1 + (eccentricity**2) / 2) ** 2
        d = self.stardist * 1.49597870700e11

        T0 = ct1 * pow(self.mass, ct2)
        gamma = 1.15 + 0.05 * pow(cg1 / self.mass, cg2)
        Teq = pow(
            (1 - A) * luminosity * LS / (16 * np.pi * (d**2) * sigmaSB * B), 1.0 / 4
        )

        Rp = self._radius
        self._radius = Rp * (1 + 0.05 * pow(Teq / T0, gamma))

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
        dsync = synchro_dist(
            self.rotrate * wrot_J, Qpp, age, self.mass, self.radius, star_mass
        )
        if d <= dsync:
            self.rotrate = self.orbitperiod
        else:
            self.rotrate = self.rotrate

    def calculate_luminosity(
        self,
        age: float,
        table1: pd.DataFrame,
        table2: pd.DataFrame,
        table3: pd.DataFrame,
    ):
        """Retourne la liste des luminosités en fonction de l'âge selon la valeur de M"""
        # MJ = 1.8986e27
        M = self.mass
        print("Mass ds calc", M)
        L1 = np.interp(
            np.log10(age * 1e-9), table1["log(t) (Gyr)"], table1["log(L/Ls)"]
        )
        L2 = np.interp(
            np.log10(age * 1e-9), table2["log(t) (Gyr)"], table2["log(L/Ls)"]
        )
        L3 = np.interp(
            np.log10(age * 1e-9), table3["log(t) (Gyr)"], table3["log(L/Ls)"]
        )
        luminosities = [L1, L2, L3]
        masses = np.log10([1.0, 5.0, 10.0])
        L = np.interp(np.log10(M), masses, luminosities)
        self.luminosity = 10**L
        return 10**L

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
        if Rmean:
            print("Rmean : ", R)
            return np.mean(R)
        if Rmax:
            return np.max(R)
