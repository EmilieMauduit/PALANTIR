#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:02:05 2021

@author: Emilie Mauduit
"""

from math import pow
import numpy as np
from typing import List

# --------------------------------------------------------- #
# ------------ Useful functions for the class ------------- #
# --------------------------------------------------------- #


def TOUT(mass: float):
    "From Tout et al, 1996"
    theta2 = 1.71535900
    iota2 = 6.59778800
    kappa2 = 10.08855000
    lambda2 = 1.01249500
    mu2 = 0.07490166
    nu2 = 0.01077422
    eta2 = 3.08223400
    omega2 = 17.84778000
    pi2 = 0.00022582
    a = (
        theta2 * pow(mass, 2.5)
        + iota2 * pow(mass, 6.5)
        + kappa2 * pow(mass, 11)
        + lambda2 * pow(mass, 19)
        + mu2 * pow(mass, 19.5)
    )
    b = (
        nu2
        + eta2 * pow(mass, 2)
        + omega2 * pow(mass, 8.5)
        + pow(mass, 18.5)
        + pi2 * pow(mass, 19.5)
    )
    return a / b


# ============================================================= #
# ---------------------------- Star --------------------------- #
# ============================================================= #


class Star:
    def __init__(
        self, name: str, mass: float, radius: dict, age: float, obs_dist: float
    ):

        """Creates a Star object.
        :param name:
            Name of the star.
        :type name:
            str
        :param M:
            Star mass, in Sun masses.
        :type M:
            float
        :param R:
            Star radius, in Sun radiuses.
        :type R:
            float
        :param t:
            Star age, in yr.
        :type t:
            float
        :param s:
            Distance from Earth, in pc.
        :type s:
            float
        :param B:
            Star magnetic field, in T. Either known from litterature or computed.
        :type B:
            float
        """

        self.name = name
        self.mass = mass
        self.radius = radius
        self.age = age * 1e9
        self.obs_dist = obs_dist
        self._rotperiod = None
        self._magfield = None
        self._luminosity = None

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
    def rotperiod(self):
        if self._rotperiod is None:
            self._rotperiod = self._compute_rotperiod(age=self.age)
        return self._rotperiod

    @property
    def magfield(self):
        if self._magfield is None:
            self._magfield = self._compute_magfield(self.rotperiod)
        return self._magfield

    @property
    def luminosity(self):
        a = 0.39704170
        b = 8.52762600
        c = 0.00025546
        d = 5.43288900
        e = 5.56357900
        f = 0.78866060
        g = 0.00586685
        M = self.mass
        res1 = a * pow(M, 5.5) + b * pow(M, 11)
        res2 = (
            c
            + pow(M, 3)
            + d * pow(M, 5)
            + e * pow(M, 7)
            + f * pow(M, 8)
            + g * pow(M, 9.5)
        )
        self._luminosity = res1 / res2
        return self._luminosity

    def talk(self, talk: bool):
        if talk:
            print("Star name : ", self.name)
            print("Star mass : ", self.mass, " MS")
            print("Star radius : ", self.radius, " RS")
            print("Star age : ", self.age * 1e-9, " Gyr")
            print("Star distance to Earth : ", self.obs_dist, " pc")
            print("Star rotational period : ", self.rotperiod, " days")
            print("Star magnetic field :  ", self.magfield, " T")
            print("Star luminosity : ", self.luminosity, "LS")

    # Methods

    def unnormalize_mass(self) -> float:
        MS = 1.989e30  # kg
        M = self.mass * MS
        return M

    def unnormalize_radius(self) -> float:
        RS = 6.96342e8  # m
        R = self.radius * RS
        return R

    def unnormalize_luminosity(self) -> float:
        LS = 3.826e26  # W
        L = self.luminosity * LS
        return L

    def obs_dist_meters(self) -> float:
        pc = 3.08568e16  # m
        return self.obs_dist * pc

    def alfven_radius(self, d: float) -> float:
        """Computes the Alfvén radius of the star.
        :param d:
            Distance between the star and the planet, in m.
        :type d:
            float
        """

        Ra = 1

        return Ra

    @staticmethod
    def _calculate_radius(
        models: List[str], mass: float, Rmean: bool = True, Rmax: bool = False
    ):
        R = []
        for model in models:
            if "Tout" in model:
                R.append(TOUT(mass))
        if Rmean:
            return np.mean(R)
        if Rmax:
            return np.max(R)

    @staticmethod
    def _compute_rotperiod(age: float):
        """Define the rotational period of the star in days."""
        tau = 2.56e7  # yr
        K = 0.6709505359255223
        rotperiod = K * pow(1 + (age / tau), 0.7)
        return rotperiod

    @staticmethod
    def _compute_magfield(rotperiod: float):
        Psun = 25.5  # days
        Bsun = 1.435e-4  # T
        magfield = Bsun * Psun / rotperiod
        return magfield
