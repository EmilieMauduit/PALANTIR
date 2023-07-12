#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: emauduit
"""


from array import array
import numpy as np
from math import pow


def find_value(value: float, list: array, res: bool = False):
    """Find the index of the value closest to a in a list l.
    :param a:
        Value to find in the list
    :type a:
        float
    :param l:
        List to search into
    :type l:
        array
    :param res:
        If not specified, the function returns the index. If res=True, returns the value associated to the index.
    """
    min = np.abs(list[0] - value)
    indmin = 0
    for i in range(1, len(list)):
        nmin = np.abs(list[i] - value)
        if nmin <= min:
            min = nmin
            indmin = i
    if not res:
        return indmin
    else:
        return list[indmin]


def synchro_dist(
    wrot: float,
    Qp: float,
    tsync: float,
    planet_mass: float,
    planet_radius: float,
    star_mass: float,
):
    """Compute the distance at which a planet should be so that its orbit is synchronised, depending on the age of the system.
    :param wrot:
        orbital period
    :type wrot:
        float
    :param Qp:
        factor of dissipation by tidal effect
    :type Qp:
        float
    :param tsync:
        time of synchronization
    :type tsync:
        float
    :param planet_mass:
        Mass of the planet [Mjup]
    :type planet_mass:
        float
    :param planet_radius:
        Radius of the planet [Rjup]
    :type planet_radius:
        float
    :param star_mass:
        Mass of the star [Msun]
    :type star_mass:
        float
    """
    G = 6.6725985e-11
    a = 9 * tsync / (4 * wrot * 0.26 * Qp)
    b = G * planet_mass / pow(planet_radius, 3)
    c = pow(star_mass / planet_mass, 2)
    res = planet_radius * pow(a * b * c, 1 / 6)
    return res


def calc_vsun(age: float):
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


def calc_nsun(age: float):
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


def calc_Bimf(stardist: float):
    Br0 = 2.6e-9
    Bp0 = 2.4e-9  # T
    Br = Br0 * pow(stardist, -2)
    Bp = Bp0 * pow(stardist, -1)
    return (Br, Bp)
