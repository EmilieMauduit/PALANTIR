#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:46:05 2021

@author: Emilie Mauduit
"""

import numpy as np
from typing import List
from math import sqrt, pow

from planet import Planet
from star import Star
from dynamo_region import DynamoRegion
from stellar_wind import StellarWind

import logging
log = logging.getLogger('palantir.prediction_tools.magnetic_moment')

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
        self.standoff_dist = Rs

    def __str__(self):
        models_str =''
        for m in self.models :
            models_str = models_str + m + ", "

        return (
            "Model(s) used: " + models_str + ". \n"
            + "Magnetic moment, M={} MmJ. \n".format(self.mag_moment)
            + "Standoff distance, Rs={} m.\n".format(self.standoff_dist)
        )

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def normalize_standoff_dist(self, planet: Planet):
        return self.standoff_dist / planet.unnormalize_radius()

    def unormalize_magnetic_moment(self, other):
        return self.mag_moment * other.mag_moment

    def magnetic_moment(
        self,
        dynamo: DynamoRegion,
        planet: Planet,
        jup: Planet,
        dynamo_jup : DynamoRegion,
        Mmean: bool = True,
        Mmax: bool = False,
        normalize: bool = True,
    ):

        """Creates a MagneticMoment object associated to a given dynamo region.
        :param dynamo:
            Dynamo region from which the magnetic moment will be computed.
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
            Default is True. If you give more than one model, the geometric mean value of the moments will be returned.
        :type Mmean:
            bool
        :param Mmax:
            Default is False. If you give more than one model, the maximum value of the moments will be returned.
        :type Mmax:
            bool
        """
        M = []

        if normalize:
            wrot = planet.rotrate
        else:
            wrot = planet.unnormalize_rotrate()

        for model in self.models:
            if model == "blackett":
                M.append(self._Blackett(dynamo.radius, wrot, dynamo.density))
            if model == "busse":
                M.append(self._Busse(dynamo.radius, wrot, dynamo.density))
            if model == "mizu_mod":
                M.append(self._Mizu_moderate(dynamo.radius, wrot, dynamo.density, 1.0))
            if model == "mizu_slow":
                M.append(self._Mizu_slow(dynamo.radius, wrot, dynamo.density, 1.0))
            if model == "sano":
                M.append(self._Sano(dynamo.radius, wrot, dynamo.density))
            if model == "rein-chris":
                if planet.mass >= 0.17:
                    M.append(self._Reiners_Christensen(planet, jup, dyn_region=dynamo, dynamo_jup = dynamo_jup))
                else:
                    log.warning("Warning : Planet mass is lower than 0.17 MJ. Magnetic moment has been set to 0.")
                    print(
                        "Warning : Planet mass is lower than 0.17 MJ. Magnetic moment has been set to 0."
                    )
                    M.append(np.nan)
            if model == "rein-chris-dyn":
                if planet.mass >= 0.17:
                    M.append(self._Reiners_Christensen(planet, jup, dyn_region=dynamo,dynamo_jup = dynamo_jup))
                else:
                    log.warning("Warning : Planet mass is lower than 0.17 MJ. Magnetic moment has been set to 0.")
                    print(
                        "Warning : Planet mass is lower than 0.17 MJ. Magnetic moment has been set to 0."
                    )
                    M.append(np.nan)

        if Mmean and not Mmax:
            self.mag_moment = pow(np.prod(M), 1 / len(M))
        elif Mmax and not Mmean:
            self.mag_moment = np.max(M)
        elif not Mmean and not Mmax:
            self.mag_moment = M[0]
        else:
            raise ValueError(
                "Wrong value for Mmean :{Mmean} or Mmax : {Mmax}, only one can be set to True"
            )

    def magnetosphere_radius(self, other, stellar_wind: StellarWind):
        """Computes the radius of the magnetosphere of a given planet.
        :param stellar_wind:
            Stellar wind parameters associated with the system studied.
        :type stellar_wind:
            StellarWind
        """
        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        res1 = (mp * stellar_wind.density * pow(stellar_wind.effective_velocity, 2)) + (
            2 * stellar_wind.density * kb * stellar_wind.corona_temperature
        )
        self.standoff_dist = pow(
            (np.pi * 4e-7 * (1.16**2) * (self.unormalize_magnetic_moment(other) ** 2))
            / (res1 * 8 * (np.pi**2)),
            1 / 6,
        )

    # --------------------------------------------------------- #
    # --------------------- Static methods -------------------- #

    # Modèles du moment dipolaire magnétique

    # Scale laws

    @staticmethod
    def _Blackett(rc: float, wrot: float, rhoc: float):
        """Compute the magnetic moment using Blackett's scaling law, 1947.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        """
        return rhoc * wrot * (rc**5)

    @staticmethod
    def _Busse(rc: float, wrot: float, rhoc: float):
        """Compute the magnetic moment using Busse's scaling law, 1976.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        """
        res = pow(rhoc, 1 / 2) * wrot * pow(rc, 4)
        return res


    @staticmethod
    def _Mizu_moderate(rc: float, wrot: float, rhoc: float, sigma: float = 1.0):
        """Compute the magnetic moment using moderate Mizutani's scaling law.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        :param sigma:
            Conductivity in the dynamo region of the planet, supposed to be the same as Jupiter's. Default is 1.
        :type sigma:
            float

        """
        res = pow(rhoc, 1 / 2) * pow(wrot, 3 / 4) * pow(rc, 7 / 2) * pow(sigma, -1 / 4)
        return res

    @staticmethod
    def _Mizu_slow(rc: float, wrot: float, rhoc: float, sigma: float):
        """Compute the magnetic moment using slow Mizutani's scaling law.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        :param sigma:
            Conductivity in the dynamo region of the planet, supposed to be the same as Jupiter's. Default is 1.
        :type sigma:
            float

        """
        res = pow(rhoc, 1 / 2) * pow(wrot, 1 / 2) * pow(rc, 3) * pow(sigma, -1 / 4)
        return res

    @staticmethod
    def _Curtis(rc: float, wrot: float, rhoc: float, heat_flux: float):
        """Compute the magnetic moment using curtis's scaling law, 1986.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        :param heat_flux:
            Heat flux of the planet, difficult to estimate.
        :type heat_flux:
            float

        """
        a = pow(rhoc, 1 / 3) * pow(wrot, 1 / 2) * pow(heat_flux, 1 / 6)
        b = pow(rc, 7 / 2)
        return a * b

    @staticmethod
    def _Sano(rc: float, wrot: float, rhoc: float):
        """Compute the magnetic moment using Sano's scaling law, 1993.
        :param rc:
            Radius of the dynamo region of the planet, normalized to Jupiter's.
        :type rc:
            float
        :param wrot:
            Rotation rate of the planet, normalized to Jupiter's.
        :type wrot:
            float
        :param rhoc:
            Density of the dynamo region of the planet, normalized to Jupiter's.
        :type rhoc:
            float
        """
        a = pow(rhoc, 1 / 2) * pow(rc, 7 / 2) * wrot
        return a


    # Simulations

    @staticmethod
    def _Reiners_Christensen(planet: Planet, jup: Planet, dyn_region : DynamoRegion, dynamo_jup : DynamoRegion):

        """Computing the magnetic moment based on Reiners-Christensen's simulations, 2010.
        :param planet:
            The planet for which the magnetic moment will be computed.
        :type planet:
            Planet
        :param dyn_region:
            The dynamo region of the planet
        :type dyn_region:
            DynamoRegion
        :param dynamo_jup:
            Jupiter's dynamo region (the magnetic moment of the planet is computed normalized to Jupiter's)
        :type dynamo_jup:
            DynamoRegion
        """

        #Bdyn = dyn_region.mag_field_dynamo ; Bdyn_jup = dynamo_jup.mag_field_dynamo
        Beq = dyn_region.mag_field_equatorial
        Beq_jup = dynamo_jup.mag_field_equatorial  

        mag_moment = pow(planet.radius,3) * Beq / Beq_jup #eq 4
        
        return mag_moment