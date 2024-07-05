#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Emilie Mauduit
"""

from scipy import optimize
import numpy as np
from math import sqrt, log, atan, sin

from planet import Planet
from star import Star

# ============================================================= #
# ------------------------ StellarWind ------------------------ #
# ============================================================= #


class StellarWind:
    def __init__(self, ne: float, ve: float, Tcor: float, Bsw: dict):
        """Creates a stellar wind object.

        :param ne:
            Electron density in the stellar wind, in m-3.
        :type ne:
            float
        :param ve:
            Effective velocity of the stellar wind, in m.s-1.
        :type ve:
            float
        :param Tcor:
            Coronal temperature, in K.
        :type Tcor:
            float
        :param Bsw:
            Magnetic field of the stellar wind, in T.
        :type Bsw:
            float
        """

        self.density = ne
        self.effective_velocity = ve
        self.corona_temperature = Tcor
        self.mag_field = Bsw
        self._alfven_velocity = None
    
    def __str__(self):
        return ("Electron density : {} m-3 \n".format(self.density)
        + "Effective velocity of the stellar wind : {} m.s-1 \n".format(self.effective_velocity)
        + "Temperature of the corona : {} MK \n".format(self.corona_temperature * 1e-6)
        + "Stellar wind magnetic field : {} T \n".format(self.mag_field)
        + "Electron Alfven velocity : {} m.s-1 \n".format(self.alfven_velocity))


    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def mag_field(self):
        return self._mag_field

    @mag_field.setter
    def mag_field(self, value: dict):
        if ("planet" or "star" or "vsw") not in value:
            raise KeyError("Planet or Star or SW velocity not in value.")
        G = 6.6725985e-11  # N.m^2/kg^2
        dua = 1.49597870700e11  # m
        Psun = 25.5  # days
        vorb = sqrt(
            G * value["star"].unnormalize_mass() / (value["planet"].stardist * dua)
        )
        Bimf_r, Bimf_p = self._calc_Bimf(stardist=value["planet"].stardist,magfield_surf = value["star"].magfield)
        alpha = atan(Bimf_p / Bimf_r)
        beta = atan(vorb / value["vsw"])
        #Bimf_r *= Psun / value["star"].rotperiod
        #Bimf_p *= Psun / value["star"].rotperiod
        self._mag_field = sqrt(Bimf_r**2 + Bimf_p**2) * abs(sin(alpha - beta))

    @property
    def alfven_velocity(self):
        """ Compute the electronic Alfven velocity. """
        me = 9.1e-31 #kg
        if self._alfven_velocity is None :
            mu0 = 4*np.pi*1e-7
            self._alfven_velocity = self._mag_field / np.sqrt(mu0 * self.density * me)
        return self._alfven_velocity

    @classmethod
    def from_system(cls, star: Star, planet: Planet):
        """Creates a StellarWind object corresponding to the given Star-Planet system.
        :param star:
            The star of the system.
        :type star:
            Star
        :param planet:
            The planet of the system.
        :type planet:
            Planet
        """

        v, veff, ne, T = cls._Parker(star=star, planet=planet)

        return cls(ne, veff, T, Bsw={"planet": planet, "star": star, "vsw": v})

    # --------------------------------------------------------- #
    # -------------------- Static methods --------------------- #

    @staticmethod
    def _calc_vsun(age: float):
        """Compute the velocity of the solar wind at 1 AU.
        :param age:
            stellar age [yr]
        :type age:
            float
        """
        tau = 2.56e7  # yr
        vo = 3971e3  # m/s
        vsun = vo / pow((1 + (age / tau)), 0.43)
        return vsun

    @staticmethod
    def _calc_nsun(age: float):
        """Compute the electron density of the solar wind at 1 AU.
        :param age:
            stellar age [yr]
        :type age:
            float
        """
        tau = 2.56e7  # yr
        no = 1.04e11  # m-3
        nsun = no / pow((1 + (age / tau)), 1.86)
        return nsun

    @staticmethod
    def _calc_Bimf(stardist: float, magfield_surf : float = None):

        if magfield_surf is None :
            Br0 = 2.6e-9
            Bp0 = 2.4e-9  # T
        else :
            Br0 = magfield_surf
            Bp0 = magfield_surf
        Br = Br0 * pow(stardist, -2)
        Bp = Bp0 * pow(stardist, -1)
        return (Br, Bp)

    @staticmethod
    def _mass_lossrate(t: float, R: float) -> float:
        """Compute the stellar mass loss rate normalized to the sun
        :param t:
            stellar age [yr]
        :type t:
            float
        :param R:
            radius of the star [R_sun]
        :type R:
            float
        """
        mp = 1.660540210e-27  # kg
        do = 1.49597870e11  # m

        vsun = StellarWind._calc_vsun(t)
        nsun = StellarWind._calc_nsun(t)
        M = 4 * np.pi * (do**2) * nsun * vsun * mp * (R**2)
        return M

    @staticmethod
    def _parker_velocity(v: float, d: float, rc: float, vc: float) -> float:
        """Expression of f(v(d))=0 based on Parker's model for SW.
        :param v:
            velocity [m/s]
        :type v:
            float
        :param d:
            distance between the star and the planet [m]
        :type d:
            float
        :param rc:
            critical radius above which the solar wind become super-sonic
        :type rc:
            float
        :param vc:
            critical velocity associated to rc
        :type vc:
            float
        """
        res = (
            ((v / vc) ** 2)
            - log((v / vc) ** 2)
            - log((d / rc) ** 4)
            - (4.0 / (d / rc))
            + 3.0
        )
        return res

    @staticmethod
    def _parker_velocity_derivative(v: float, d: float, rc: float, vc: float) -> float:
        """Expression of df(v(d))/dv based on Parker's model for SW.
        :param v:
            velocity [m/s]
        :type v:
            float
        :param d:
            distance between the star and the planet [m]
        :type d:
            float
        :param rc:
            critical radius above which the solar wind become super-sonic
        :type rc:
            float
        :param vc:
            critical velocity associated to rc
        :type vc:
            float
        """
        res = (2 * v / (vc**2)) - (2 / v)
        return res

    @staticmethod
    def _calc_temperature(M: float, t: float) -> float:
        """Adjust the temperature of the corona for a star of given age.

        :param M:
            Mass of the star in kg
        :type M:
            float
        :param t:
            Age of the star in yr
        :type t:
            float

        """

        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        d = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2

        vsun = StellarWind._calc_vsun(t)
        # nsun = calc_nsun(t)

        Tini = (1e6 * 3.6e9) / (t + 1e9)
        vc = sqrt(2 * kb * Tini / mp)
        rc = mp * G * M / (4 * kb * Tini)
        v = optimize.newton(
            StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
        )
        if v <= vsun:
            while (v >= (vsun + 2e3)) or (v <= (vsun - 2e3)):
                Tini = Tini + 0.005e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
                )
            if v <= vsun:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini + 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
            else:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini - 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini

        else:
            while (v >= (vsun + 2e3)) or (v <= (vsun - 2e3)):
                Tini = Tini - 0.005e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
                )
            if v <= vsun:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini + 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
            else:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini - 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
        return T

    @staticmethod
    def _Parker(star: Star, planet: Planet, T: float = None):
        """Compute the velocity and the density of the SW using the Parker model.
        :param planet:
            The planet studied
        :type planet:
            Planet
        :param star:
            The associated star
        :type star:
            Star
        """

        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        dua = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2
        d = planet.stardist
        t = star.age  # yr
        if T is None:
            T = StellarWind._calc_temperature(star.unnormalize_mass(), t)
        vc = sqrt(2 * kb * T / mp)
        rc = mp * G * star.unnormalize_mass() / (4 * kb * T)
        vorb = sqrt(G * star.unnormalize_mass() / (d * dua))
        # Warning messages
        if d < 0.01:
            print(
                "Warning : The Parker model has not been verified for such star-planet distances"
            )
        if t < 0.7e9:
            print("Warning : The Parker model is not precise for stars with t<0.7 Gyr")

        if d >= 1.0:
            v = optimize.newton(
                StellarWind._parker_velocity,
                350.0e3,
                StellarWind._parker_velocity_derivative,
                args=(d * dua, rc, vc),
                maxiter=50,
            )
            print("planet dist > 1UA")
        else:
            d_temp = 1.0
            v_temp_ini = optimize.newton(
                StellarWind._parker_velocity,
                350.0e3,
                StellarWind._parker_velocity_derivative,
                args=(d_temp * dua, rc, vc),
                maxiter=50,
            )
            while abs(d_temp - d) >= 1e-5:
                d_temp = (9 * d_temp + d) / 10.0
                v_temp = optimize.newton(
                    StellarWind._parker_velocity,
                    v_temp_ini,
                    StellarWind._parker_velocity_derivative,
                    args=(d_temp * dua, rc, vc),
                    maxiter=50,
                )
                if (v_temp / vc > 1.0) and (d_temp / (rc / dua) <= 1.0):
                    v_temp = optimize.newton(
                        StellarWind._parker_velocity,
                        0.9 * v_temp_ini,
                        StellarWind._parker_velocity_derivative,
                        args=(d_temp * dua, rc, vc),
                        maxiter=50,
                    )
                v_temp_ini = v_temp
            v = v_temp

        veff = sqrt((v**2) + (vorb**2))
        Mls = StellarWind._mass_lossrate(t, star.radius)
        n = Mls / (4 * np.pi * ((d * dua) ** 2) * v * mp)
        if n < 0:
            raise ValueError("Negative stellar wind density is not physical.")
        return (v, veff, n, T)

    @staticmethod
    def _CME(star: Star, planet: Planet):
        """Evaluates if it is necessary to take into account CME in the stellar wind model.
        Returns a tuple containing the speed, the effective speed and the density of the stellar wind.

        :param star:
            Star of the system considered
        :type star:
            Star
        :param pla:
            Planet of the system considered
        :type pla:
            Planet
        """
        T = 2e6  # K
        dua = 1.49597870700e11  # m
        d = planet.stardist * dua
        G = 6.6725985e-11  # N.m^2/kg^2
        # nwo = 4.88e6
        nso = 7.1e6
        # nw = nwo * pow((d / dua), -2.31)
        ns = nso * pow((d / dua), -2.99)
        v = 5.26e5
        veff = sqrt((v**2) + (star.unnormalize_mass() * G / d))
        return (veff, ns, T)

