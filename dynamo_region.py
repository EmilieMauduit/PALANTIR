#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:36:08 2021

@author: Emilie Mauduit
"""

import numpy as np
from scipy import optimize
import math as m

from planet import Planet

# --------------------------------------------------------- #
# ------------ Useful functions for the class ------------- #
# --------------------------------------------------------- #

# Models to compute the radius of the dynamo region


def rlin(a, Rp):
    return a * Rp


def rhoLE(r, Mp, Rp, rhot):
    a = np.pi * Mp / (4 * pow(Rp, 3))
    b = np.pi * r / Rp
    res = a * np.sin(b) / b
    return res - rhot


def rhoLEp(r, Mp, Rp, rhot):
    a = np.pi * Mp / (4 * pow(Rp, 3))
    res = a * (
        (np.cos(np.pi * r / Rp) / r)
        - (np.sin(np.pi * r / Rp) * Rp / (np.pi * pow(r, 2)))
    )
    return res


def LaneEmden(Mp, Rp, rhot):
    if Mp <= 5 * 1.8986e27:
        try:
            res = optimize.newton(
                rhoLE, Rp / 2, fprime=rhoLEp, args=(Mp, Rp, rhot), maxiter=50
            )
        except (ZeroDivisionError, RuntimeError):
            res = 0.0
    else:
        res = optimize.newton(
            rhoLE, 0.93 * Rp, fprime=rhoLEp, args=(Mp, Rp, rhot), tol=1.0, maxiter=50
        )

    return res


# Models to compute the density of the dynamo region


def rhodyn(Mp, Rp, rc):
    rhom = 3 * Mp / (4 * np.pi * pow(Rp, 3))
    a = rc / Rp
    if a == 0.0:
        return 0.0
    else:
        res = rhom * pow(a, -3) * ((m.sin(a * np.pi) / np.pi) - a * m.cos(a * np.pi))
    return res


# ============================================================= #
# ------------------------ DynamoRegion ----------------------- #
# ============================================================= #


class DynamoRegion:
    def __init__(self, rhocrit: float, rhoc: float, rc: float):
        """Creates a DynamoRegion object.
        :param rhocrit:
            Critical density, corresponds to the transition between the molecular phase and the liquid metallic phase.
        :type rhocrit:
            float
        :param rhoc:
            Density of the dynamo region.
        :type rhoc:
            float
        :param Rc:
            Radius of the dynamo region.
        :type Rc:
            float
        """
        self.critical_density = rhocrit
        self.density = rhoc
        self.radius = rc

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self, talk: bool):
        if talk:
            print(
                "Critical transition density, rho_crit=",
                self.critical_density,
                " g.cm-3.",
            )
            print("Density of the dynamo region, rhoc=", self.density, " g.cm-3.")
            print("Radius of the dynamo region, Rc=", self.radius, " m.")

    @classmethod
    def from_planet(cls, planet: Planet, rhocrit: float):
        """Creates a DynamoRegion object associated to a given Planet object.
        :param planet:
            Planet of the system.
        :type planet:
            Planet
        :param rhocrit:
            Critical density, corresponds to the transition between the molecular phase and the liquid metallic phase.
        :type rhocrit:
            float
        """
        Rc = LaneEmden(planet.unnormalize_mass(), planet.unnormalize_radius(), rhocrit)
        rhoc = rhodyn(planet.unnormalize_mass(), planet.unnormalize_radius(), Rc)
        return cls(rhocrit, rhoc, Rc)
