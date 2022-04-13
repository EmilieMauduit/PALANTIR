#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: emauduit
"""

import os
import scipy
import scipy.io as sio
from scipy.io import readsav
from scipy import optimize
import numpy as np
from math import *
from sympy.solvers import solve
from sympy import Symbol

from classPlanet import *
from classStar import *

from calc_tools import *

def masslossrate(t,R):
    """ Compute the stellar mass loss rate normalized to the sun"""
    tau=2.56e7 #yr
    mp=1.660540210e-27 #kg
    vo=3397e3 #m/s
    no=1.6e10 #m-3
    do=1.49597870e11 #m

    vsun=vo/pow((1+(t/tau)),0.4)
    nsun=no/pow((1+(t/tau)),1.5)
    M=4*np.pi*(do**2)*nsun*vsun*mp*(R**2)
    return(M)

def parker_velocity(v,d,rc, vc):
    """Expression of f(v(d))=0 based on Parker's model for SW"""
    res=((v/vc)**2)-log((v/vc)**2)-log((d/rc)**4)-(4.0/(d/rc))+3.0 
    return(res)

def parker_velocity_p(v,d,rc,vc):
    res=(2*v/(vc**2))-(2/v)
    return(res)

def parker(pla,star,T,vo,limd):
    """Compute the velocity and the density of the SW using the Parker model"""
    kb=1.380658e-23 #J/K
    mp=1.660540210e-27 #kg
    G=6.6725985e-11 #N.m^2/kg^2
    M=star.mass
    R=star.radius
    d=pla.stardist
    t=star.age
    vc=sqrt(2*kb*T/mp)
    rc=mp*G*M/(4*kb*T)
    vorb=sqrt(G*M/d)
    if (d >= limd):
        try:
            v=optimize.newton(parker_velocity,350.0e3,parker_velocity_p, args=(d,rc,vc), maxiter=50)
        except (ZeroDivisionError, RuntimeError):
            v=0.
    else :
        try:
            v=optimize.newton(parker_velocity,10.0e3,parker_velocity_p, args=(d,rc,vc), maxiter=100)
        except (ZeroDivisionError, RuntimeError):
            v=0.

    veff=sqrt((v**2)+(vorb**2))
    Mls=masslossrate(t,R)
    n=Mls/(4*np.pi*(d**2)*v*mp)
    return(v,veff,n)