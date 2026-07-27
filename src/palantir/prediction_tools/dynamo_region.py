#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:36:08 2021

@author: Emilie Mauduit
"""

import numpy as np
from scipy import optimize
import math as m

import logging
log = logging.getLogger('palantir.prediction_tools.dynamo_region')

from palantir.prediction_tools.planet import Planet


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
        self.mag_field_dynamo = None
        self.mag_field_equatorial = None

    def __str__(self):
        return(
            "Critical transition density, rho_crit={} g.cm-3.\n".format(self.critical_density)
            + "Density of the dynamo region, rhoc={} rhocJ.\n".format(self.density)
            + "Radius of the dynamo region, Rc={} RcJ \n".format(self.radius)
            + "Dynamo magnetic field Bdyn={} T \n".format(self.mag_field_dynamo)
            + "Equatorial magnetic field Beq={} T".format(self.mag_field_equatorial)
        )


    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def normalize(self, other):
        self.density = self.density / other.density
        self.radius = self.radius / other.radius

    def magnetic_field(self, planet : Planet, rc_dyn : bool = False , jup : bool = False) :
        """
        Dynamo and equatorial magnetic field computed using eq (1) and (2) of Reiners & Christensen, 2010.
        """
        MS = 1.989e30  # kg
        RS = 6.96342e8  # m
        RJ = 69911e3  # m
        rcJ = 0.8487819514093978 #Rj
        Bdyn = 4.8 * pow(
            (planet.unnormalize_mass() / MS)
            * (planet.luminosity ** 2)
            * pow(RS / planet.unnormalize_radius(), 7),
            1.0 / 6,
        ) * 1e-1 #T
        self.mag_field_dynamo = Bdyn 

        if jup :
            if rc_dyn : 
                self.mag_field_equatorial = pow(self.radius / RJ ,3) * Bdyn  / (2 * m.sqrt(2))
            else :
                if planet.mass > 13. :
                    self.mag_field_equatorial = Bdyn / (2 * m.sqrt(2))
                else :
                    self.mag_field_equatorial = pow(1 - (0.17 / planet.mass), 3) * Bdyn / (2 * m.sqrt(2))
        else :
            if rc_dyn :
                self.mag_field_equatorial = pow(self.radius / RJ ,3) * Bdyn / (2 * m.sqrt(2))
            else :
                if planet.mass > 13. :
                    self.mag_field_equatorial = Bdyn / (2 * m.sqrt(2))
                elif planet.mass < 0.17 :
                    self.mag_field_equatorial = np.nan
                    self.mag_field_dynamo = np.nan
                else :
                    self.mag_field_equatorial = pow(1 - (0.17 / planet.mass), 3) * Bdyn / (2 * m.sqrt(2))

    @classmethod
    def from_planet(cls, planet: Planet, rhocrit: float):
        """
        Creates a DynamoRegion object associated to a given Planet object.
        :param planet:
            Planet of the system.
        :type planet:
            Planet
        :param rhocrit:
            Critical density, corresponds to the transition between the molecular phase and the liquid metallic phase.
        :type rhocrit:
            float
        """
        Rc = cls._LaneEmden(planet.unnormalize_mass(), planet.unnormalize_radius(), rhocrit)
        rhoc = cls._rhodyn(planet.unnormalize_mass(), planet.unnormalize_radius(), Rc)
        return cls(rhocrit, rhoc, Rc)

    # --------------------------------------------------------- #
    # -------------------- Static methods --------------------- #

    # Models to compute the radius of the dynamo region

    @staticmethod
    def _rlin(a, Rp):
        return a * Rp
    
    @staticmethod
    def _rhoLE(r, Mp, Rp, rhot):
        a = np.pi * Mp / (4 * np.power(Rp, 3))
        b = np.pi * r / Rp
        res = a * np.sin(b) / b
        return res - rhot
    
    @staticmethod
    def _rhoLEp(r, Mp, Rp, rhot):
        a = np.pi * Mp / (4 * np.power(Rp, 3))
        res = a * (
            (np.cos(np.pi * r / Rp) / r)
            - (np.sin(np.pi * r / Rp) * Rp / (np.pi * pow(r, 2)))
        )
        return res
    

    @staticmethod
    def _LaneEmden(Mp, Rp, rhot):
        """ Compute the radius of the dynamo region, if existing, by using Lane-Emden equation.
        If the density of the planet is smaller than the chosen transition density, value is set to NaN.
        """
        if np.pi*Mp/(4*(Rp**3)) < rhot :
            return np.nan
        
        else : 
            radius = np.linspace(0.,1,1000) * Rp
            rho = DynamoRegion._rhoLE(radius,Mp,Rp,rhot)
            indexes = np.argwhere(rho <= 0.)[0]
            return radius[indexes[0]]

    @staticmethod
    def _LaneEmden_old(Mp, Rp, rhot):
        """ Compute the radius of the dynamo region, if existing, by using Lane-Emden equation.
        If the density of the planet is smaller than the chosen transition density, value is set to zero.
        """
        if np.pi*Mp/(4*(Rp**3)) < rhot :
            return 0.
        else :
            RJ = 69911e3  # m
            radius = np.linspace(0.0,Rp/RJ,100) * RJ
            rho = DynamoRegion._rhoLE(radius,Mp,Rp,rhot)
            indexes = np.argwhere(rho <= 0.)[0]
            try:
                res = optimize.newton(
                    DynamoRegion._rhoLE, radius[indexes[0]], fprime=DynamoRegion._rhoLEp, args=(Mp, Rp, rhot), maxiter=50
                )
            except (ZeroDivisionError, RuntimeError):
                res = np.nan

        return res

    # Models to compute the density of the dynamo region

    @staticmethod
    def _rhodyn(Mp, Rp, rc):
        rhom = 3 * Mp / (4 * np.pi * pow(Rp, 3))
        a = rc / Rp
        if a == 0.0:
            return 0.0
        else:
            res = rhom * pow(a, -3) * ((m.sin(a * np.pi) / np.pi) - a * m.cos(a * np.pi))
        return res
