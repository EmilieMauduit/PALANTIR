#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:46:05 2021

@author: emauduit
"""

from logging import raiseExceptions
from multiprocessing.dummy import Value
from Planet import Planet
import numpy as np
import pandas as pd

from typing import List

from planet import Planet
from star import Star
from dynamo_region import DynamoRegion




# --------------------------------------------------------- #
# ------------ Useful functions for the class ------------- #
# --------------------------------------------------------- #


# Modèles du moment dipolaire magnétique

# Scale laws


def blackett(rc, w, rhoc):
    """Modèle de Blackett, 1947"""
    return rhoc * w * (rc**5)


def Busse(rc, w, rhoc):
    """Modèle de Busse, 1976"""
    res = pow(rhoc, 1 / 2) * w * pow(rc, 4)
    return res


def Mizu_moderate(rc, w, rhoc, sig):
    """Modèle de Mizutani avec convection modérée"""
    res = pow(rhoc, 1 / 2) * pow(w, 3 / 4) * pow(rc, 7 / 2) * pow(sig, -1 / 4)
    return res


def Mizu_slow(rc, w, rhoc, sig):
    """Modèle de Mizutani avec convection lente"""
    res = pow(rhoc, 1 / 2) * pow(w, 1 / 2) * pow(rc, 3) * pow(sig, -1 / 4)
    return res


def curtis(rc, w, rhoc, E):
    """Modèle de Curtis, 1986"""
    a = pow(rhoc, 1 / 3) * pow(w, 1 / 2) * pow(E, 1 / 6)
    b = pow(rc, 7 / 2)
    return a * b


def sano(rc, w, rhoc):
    """Modèle de Sano, 1993"""
    a = pow(rhoc, 1 / 2) * pow(rc, 7 / 2) * w
    return a


# Simulations


def Reiners_Christensen(planet:Planet, star:Star, jup:Planet, sol:Star, table1:pd.DataFrame, table2:pd.DataFrame, table3:pd.DataFrame):
    """Modèle de Reiners-Christensen, 2010."""
    # if normalize==0:
    Mp = planet.mass
    Rp = planet.radius
    t = star.age
    # -0.5e9
    LJ = jup.luminosity
    if planet.name == "Jupiter":
        L = 1.0
        print("jup")
    else:
        L = planet.luminosity / LJ
        # L=planet.calculate_luminosity(t,table1,table2,table3)/LJ
    print("dans RC : L=", L)
    # else :
    #    Mp=planet.normalize_mass(jup)
    #    Rp=planet.normalize_radius(jup)
    #    t=star.age
    #    L=planet.calculate_luminosity(t,table1,table2,table3)/jup.luminosity

    a = Mp * (L**2) * pow(Rp, 11)
    b = pow(1 - (0.17 / Mp), 3) / pow((1 - 0.17), 3)
    res = b * pow(a, 1.0 / 6)
    return res


# ============================================================= #
# ----------------------- MagneticMoment ---------------------- #
# ============================================================= #


class MagneticMoment:
    def __init__(self, models: List[str], Mm: float, Rs: float):
        """Creates a MgneticMoment object, characterized by a value of the magnetic moment and the corresponding radius of the magnetosphere.

        :param models:
            List of the models used to compute the magnetic moment. List of the names : 'blackett', 'busse', 'mizu_mod', 'mizu_slow', 'sano', 'rein-chris'.
        :type models:
            list[str]
        :param Mm:
            Magnetic moment, in MmJ
        :type Mm:
            float
        :param Rs:
            Radius of the magnetosphere, in Rp.
        :type Rs:
            float
        """

        self.models = models
        self.mag_moment = Mm
        self.mag_radius = Rs

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def magnetic_moment(
        self,
        dynamo: DynamoRegion,
        planet: Planet,
        star: Star,
        Mmean: bool = True,
        Mmax: bool = False,
    ):

        """Creates a MagneticMoment object associated to a given dynamo region.
        :param dynamo:
            Dynamo region from which the magnetic moment willbe computed
        :type dynamo:
            DynamoRegion
        :param planet:
            The planet in consideration.
        :type planet:
            Planet
        :param star:
            The host star of the system.
        :type star:
            Star
        :param Mmean:
            Default is True. If you give more than one model, the mean value of the moments will be returned.
        :type Mmean:
            bool
        :param Mmax:
            Default is False. If you give more than one model, the maximum value of the moments will be returned.
        :type Mmax:
            bool
        """
        M = []

        for model in self.models:
            if model == "blackett":
                M.append(blackett(dynamo.radius, planet.rotrate, dynamo.density))
            if model == "busse":
                M.append(Busse(dynamo.radius, planet.rotrate, dynamo.density))
            if model == "mizu_mod":
                M.append(
                    Mizu_moderate(dynamo.radius, planet.rotrate, dynamo.density, 1.0)
                )
            if model == "mizu_slow":
                M.append(Mizu_slow(dynamo.radius, planet.rotrate, dynamo.density, 1.0))
            if model == "sano":
                M.append(sano(dynamo.radius, planet.rotrate, dynamo.density))
            if model == "rein-chris":
                table_1MJ = pd.read_csv(
                    r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/1MJ.csv",
                    delimiter=";",
                )
                table_5MJ = pd.read_csv(
                    r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/5MJ.csv",
                    delimiter=";",
                )
                table_10MJ = pd.read_csv(
                    r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/10MJ.csv",
                    delimiter=";",
                )
                M.append(
                    Reiners_Christensen(self, star, table_1MJ, table_5MJ, table_10MJ)
                )

        if Mmean and not Mmax:
            self.mag_moment = np.mean(M)
        elif Mmax and not Mmean:
            self.mag_moment = np.max(M)
        elif not Mmean and not Mmax:
            self.mag_moment = M[0]
        else:
            raise ValueError(
                "Wrong value for Mmean :{Mmean} or Mmax : {Mmax} only one can be set to True"
            )
