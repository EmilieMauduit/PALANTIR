#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:53:23 2021

@author: Emilie Mauduit
"""


import pandas as pd
import numpy as np
from math import pow

from calc_tools import synchro_dist


# ============================================================= #
# --------------------------- Planet -------------------------- #
# ============================================================= #


class Planet:
    def __init__(
        self,
        name: str,
        mass: float,
        radius: float,
        distance: float,
        axis: float,
        worb: float,
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
        self.orbitperiod = worb
        self.rotrate = wrot
        self.stardist = distance
        self.sm_axis = axis

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self, talk: bool):
        if talk:
            print("Name : ", self.name)
            print("Mass : ", self.mass, " Mj")
            print("Radius : ", self.radius, " Rj")
            print("Rotation rate : ", self.rotrate, " wj")
            print("Orbital period : ", self.orbitperiod, "woj")
            print("Distance to host star : ", self.stardist, " AU")
            print("Semi-major axis : ", self.sm_axis, " AU")

    def unnormalize_mass(self):
        mass_jup = 1.8986e27  # kg
        mass = self.mass * mass_jup
        return mass

    def unnormalize_radius(self):
        radius_jup = 71492e3  # m
        radius = self.radius * radius_jup
        return radius

    def unnormalize_rotrate(self):
        wrot_jup = 1.77e-4  # s-1
        wrot = self.rotrate * wrot_jup
        return wrot

    def calculate_orbitalperiod(self, star_mass: float):
        """Computes the orbital period of the planet.
        :param star_mass:
            Mass of the host star of the planet.
        :type star_mass:
            float
        """
        d = self.stardist
        G = 6.6725985e-11
        self.orbitperiod = pow(star_mass * G / pow(d, 3), 1 / 2)

    def tidal_locking(self, age: float, star_mass: float):
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
        if age >= 1e10:
            Qpp = 1e5
            dsync = synchro_dist(wrot_J, Qpp, 1e10, self.mass, self.radius, star_mass)
            if d <= dsync:
                self.rotrate = self.orbitperiod
            else:
                self.rotrate = 1.0
        elif 1e8 <= age <= 1e10:
            Qppmax = 1e6
            Qppmin = 1e5
            dsync_min = synchro_dist(
                wrot_J, Qppmin, 1e10, self.mass, self.radius, star_mass
            )
            dsync_max = synchro_dist(
                wrot_J, Qppmax, 1e8, self.mass, self.radius, star_mass
            )
            ddsync = dsync_max - dsync_min
            if d <= dsync_min + ddsync:
                self.rotrate = self.orbitperiod
            else:
                self.rotrate = 1.0
        else:
            Qpp = 1e6
            dsync = synchro_dist(wrot_J, Qpp, 1e8, self.mass, self.radius, star_mass)
            if d <= dsync:
                self.rotrate = self.orbitperiod
            else:
                self.rotrate = 1.0

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
