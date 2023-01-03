#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Emilie Mauduit
"""

from scipy import optimize
import numpy as np
from math import sqrt, log

from calc_tools import calc_vsun, calc_nsun

from planet import Planet
from star import Star

# --------------------------------------------------------- #
# ------------ Useful functions for the class ------------- #
# --------------------------------------------------------- #


def masslossrate(t: float, R: float) -> float:
    """Compute the stellar mass loss rate normalized to the sun
    :param t:
        stellar age [yr]
    :type t:
        float
    :param R:
        radius of the planet [R_sun]
    :type R:
        float
    """
    mp = 1.660540210e-27  # kg
    do = 1.49597870e11  # m

    vsun = calc_vsun(t)
    nsun = calc_nsun(t)
    M = 4 * np.pi * (do**2) * nsun * vsun * mp * (R**2)
    return M


def parker_velocity(v: float, d: float, rc: float, vc: float) -> float:
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


def parker_velocity_p(v: float, d: float, rc: float, vc: float) -> float:
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


def calc_temperature(M: float, t: float) -> float:
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

    vsun = calc_vsun(t)
    # nsun = calc_nsun(t)

    Tini = (1e6 * 4.6e9) / (t + 1e9)
    vc = sqrt(2 * kb * Tini / mp)
    rc = mp * G * M / (4 * kb * Tini)
    v = optimize.newton(
        parker_velocity, vsun, parker_velocity_p, args=(d, rc, vc), maxiter=50
    )
    if v <= vsun:
        while (v >= (vsun + 2e3)) or (v <= (vsun - 2e3)):
            Tini = Tini + 0.005e6
            vc = sqrt(2 * kb * Tini / mp)
            rc = mp * G * M / (4 * kb * Tini)
            v = optimize.newton(
                parker_velocity, vsun, parker_velocity_p, args=(d, rc, vc), maxiter=50
            )
        if v <= vsun:
            while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                Tini = Tini + 0.001e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    parker_velocity,
                    vsun,
                    parker_velocity_p,
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
                    parker_velocity,
                    vsun,
                    parker_velocity_p,
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
                parker_velocity, vsun, parker_velocity_p, args=(d, rc, vc), maxiter=50
            )
        if v <= vsun:
            while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                Tini = Tini + 0.001e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    parker_velocity,
                    vsun,
                    parker_velocity_p,
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
                    parker_velocity,
                    vsun,
                    parker_velocity_p,
                    args=(d, rc, vc),
                    maxiter=50,
                )
            T = Tini
    return T


def parker(star: Star, planet: Planet):
    """Compute the velocity and the density of the SW using the Parker model.
    :param pla:
        The planet studied
    :type pla:
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
    t = star.age # yr
    T = calc_temperature(star.unnormalize_mass(), t)
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

    # Finding the rigth case

    if 0.7e9 <= t <= 1.6e9:
        dlim = 0.01
    elif 1.6e9 < t < 3.5e9:
        dlim = 0.02
    elif 3.5e9 <= t:
        dlim = 0.03
    else:
        dlim = 0.0
        print("Star age of " + star.name + " is too small")

    if d > dlim:
        try:
            v = optimize.newton(
                parker_velocity,
                350.0e3,
                parker_velocity_p,
                args=(d * dua, rc, vc),
                maxiter=50,
            )
        except (ZeroDivisionError, RuntimeError):
            v = 0.0
            print("")
    else:
        try:
            v = optimize.newton(
                parker_velocity,
                10e3,
                parker_velocity_p,
                args=(d * dua, rc, vc),
                maxiter=50,
            )
        except (ZeroDivisionError, RuntimeError):
            v = 0.0

    veff = sqrt((v**2) + (vorb**2))
    Mls = masslossrate(t, star.radius)
    n = Mls / (4 * np.pi * ((d * dua) ** 2) * v * mp)
    return (veff, n, T)


def CME(star: Star, planet: Planet):
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


# ============================================================= #
# ------------------------ StellarWind ------------------------ #
# ============================================================= #


class StellarWind:
    def __init__(self, ne: float, ve: float, Tcor: float, Bsw: float = None):
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

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self, talk: bool):
        if talk:
            print("Electron density : ", self.density, " m-3")
            print(
                "Effective velocity of the stellar wind : ",
                self.effective_velocity,
                " m.s-1",
            )
            print("Temperature of the corona : ", self.corona_temperature * 1e-6, " MK")
            print("Stellar wind magnetic field : ", self.mag_field, " T")

    # @property
    # def mag_field(self):
    #   return self.mag_field

    # @mag_field.setter
    # def mag_field(self):
    # self.mag_field = 1.0

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

        ve, ne, T = parker(star=star, planet=planet)
        #if planet.stardist <= 0.1:
            #ve_cme, ne_cme, T_cme = CME(star=star, planet=planet)
            #return cls((0.7*ne + 0.3*ne_cme) , (0.7*ve + 0.3*ve_cme), (0.7*T + 0.3*T_cme))
        #else :
        return cls(ne,ve,T)
