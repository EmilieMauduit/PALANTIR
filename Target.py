#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:02:05 2021

@author: Emilie Mauduit
"""

from math import pow
import numpy as np

# from magnetic_moment import MagneticMoment
# from planet import Planet
# from star import Star
# from stellar_wind import StellarWind


# ============================================================= #
# --------------------------- Target -------------------------- #
# ============================================================= #


class Target:
    def __init__(
        self,
        name: str,
        mag_field: dict,
        pow_emission : dict,
        pow_received: dict,
        fmax_star: float,
    ):
        """Creates a Target object.
        :param name:
            Name of the target.
        :type name:
            str
        :param Bmax:
            Maximum magnetic field estimated, in T.
        :type Bmax:
            float
        :param fmax:
            Maximum frequency of the emission at the planet, in Hz.
        :type fmax:
            float
        :param Pem:
            Power of the emission, in W.
        :type Pem:
            float
        :param Prec:
            Power received by the planet from the star, in W.
        :type Prec:
            float
        :param F:
            Flux of the emission at Earth, in W.m-2
        :type F:
            float
        :param fmax_star:
            Maximum frequency of the emission at the star, in Hz.
        :type fmax_star:
            float
        """

        self.name = name
        self.mag_field = mag_field
        self._freq_max = None
        self.pow_emission = pow_emission
        self._flux = None
        self.pow_received = pow_received
        self.freq_max_star = fmax_star

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def mag_field(self):
        return self._mag_field

    @mag_field.setter
    def mag_field(self, value: dict):
        if "planet" not in value or "magnetic_moment" not in value:
            raise KeyError("planet or magnetic_moment not in value")
        planet = value["planet"]
        magnetic_moment = value["magnetic_moment"]
        mag_moment_jup = 1.56e27  # A.m-2
        self._mag_field = (
            1e-7
            * 2
            * magnetic_moment.mag_moment
            * mag_moment_jup
            / pow(planet.unnormalize_radius(), 3)
        )

    @property
    def freq_max(self):
        if self._freq_max is None:
            me = 9.1093897e-31  # kg
            e = 1.60217733e-19  # C
            self._freq_max = e * self.mag_field / (2 * np.pi * me)
        return self._freq_max

    @property
    def pow_emission(self):
        return self._pow_emission

    @pow_emission.setter
    def pow_emission(self, value: dict):
        if (
            "planet" not in value
            or "star" not in value
            or "magnetic_moment" not in value
            or "stellar_wind" not in value
        ):
            raise KeyError(
                "planet or star or magnetic_moment or stellar_wind not in value"
            )
        prad_jup = 2.1e11  # W
        standoff_dist_jup = 40.1  # RJ
        density_jup = 1.98e5  # m-3
        veff_jup = 523e3  # m/s

        self._pow_emission = (
            prad_jup
            * pow(value["planet"].radius,2)
            * pow((value["magnetic_moment"].normalize_standoff_dist(planet=value["planet"]) / standoff_dist_jup), 2)
            * (value["stellar_wind"].density / density_jup)
            * pow(value["stellar_wind"].effective_velocity / veff_jup, 3)
        )

    @property
    def flux(self):
        if self._flux is None :
            dua = 1.49597870700e11  # m
            self._flux = self._pow_emission / (1.6 * self.freq_max * (dua ** 2))
        return self._flux

    @property
    def pow_received(self):
        return self._pow_received

    @pow_received.setter
    def pow_received(self, value: dict):
        if ("star" not in value):
            raise KeyError(
                "star not in value"
            )

        pc = 3.08568e16  # m
        self._pow_received = self._pow_emission / (1.6 * self.freq_max * pow(value['star'].obs_dist * pc, 2))


    def talk(self, talk: bool):
        if talk:
            print("Name of the system : ", self.name)
            print("Maximum magnetic field : ", self.mag_field, " T")
            print(
                "Maximum frequency emission at the planet : ",
                self.freq_max * 1e-6,
                " MHz",
            )
            print("Power of the emission at the star : ", self.pow_emission/1e14, ".10^14 W")
            print( "Flux of the emission emitted at 1 AU : ",
                self.flux / 1e-26 / 1e10,
                ".10^10 Jy ",
            )
            print(
                "Flux of the emission received by the instrument : ",
                self._pow_received * 1e3 / 1e-26 ,
                " mJy ",
            )
            print(
                "Maximum frequency of the emission at th star : ",
                self.freq_max_star * 1e-6,
                " MHz",
            )

    def select_target(self) -> bool:
        """Select or not a target according to various criterions set by the user."""
        if ():
            return True
        else:
            return False
