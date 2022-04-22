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

def masslossrate(t:float,R:float):
    """ Compute the stellar mass loss rate normalized to the sun
        :param t:
            stellar age [yr]
        :type t:
            float
        :param R:
            radius of the planet [m]
        :type R:
            float
    """
    mp=1.660540210e-27 #kg
    do=1.49597870e11 #m

    vsun=calc_vsun(t)
    nsun=calc_nsun(t)
    M=4*np.pi*(do**2)*nsun*vsun*mp*(R**2)
    return(M)


def parker_velocity(v:float,d:float,rc:float, vc:float):
    """Expression of f(v(d))=0 based on Parker's model for SW
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
    res=((v/vc)**2)-log((v/vc)**2)-log((d/rc)**4)-(4.0/(d/rc))+3.0 
    return(res)

def parker_velocity_p(v:float,d:float,rc:float,vc:float):
    """Expression of df(v(d))/dv based on Parker's model for SW
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
    res=(2*v/(vc**2))-(2/v)
    return(res)

def calc_temperature(M:float,t:float):
    """ Adjust the temperature of the corona for a star of given age.
    
        :param M:
            Mass of the star in kg
        :type M:
            float
        :param t:
            Age of the star in yr
        :type t:
            float

    """

    kb=1.380658e-23 #J/K
    mp=1.660540210e-27 #kg
    d=1.49597870700e11 #m
    G=6.6725985e-11 #N.m^2/kg^2

    vsun=calc_vsun(t)
    nsun=calc_nsun(t)

    Tini=(1e6*4.6e9)/(t+1e9)
    vc=sqrt(2*kb*Tini/mp)
    rc=mp*G*M/(4*kb*Tini)
    v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
    if (v <= vsun):
        while ((v >= (vsun+2e3)) or (v <= (vsun-2e3))): 
            Tini=Tini+0.005e6
            vc=sqrt(2*kb*Tini/mp)
            rc=mp*G*M/(4*kb*Tini)
            v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
        if (v <= vsun):
            while ((v >= (vsun+0.5e3)) or (v <= (vsun-0.5e3))): 
                Tini=Tini+0.001e6
                vc=sqrt(2*kb*Tini/mp)
                rc=mp*G*M/(4*kb*Tini)
                v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
            T=Tini
        else:
            while ((v >= (vsun+0.5e3)) or (v <= (vsun-0.5e3))): 
                Tini=Tini-0.001e6
                vc=sqrt(2*kb*Tini/mp)
                rc=mp*G*M/(4*kb*Tini)
                v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
            T=Tini

    else:
        while ((v >= (vsun+2e3)) or (v <= (vsun-2e3))): 
            Tini=Tini-0.005e6
            vc=sqrt(2*kb*Tini/mp)
            rc=mp*G*M/(4*kb*Tini)
            v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
        if (v <= vsun):
            while ((v >= (vsun+0.5e3)) or (v <= (vsun-0.5e3))): 
                Tini=Tini+0.001e6
                vc=sqrt(2*kb*Tini/mp)
                rc=mp*G*M/(4*kb*Tini)
                v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
            T=Tini
        else:
            while ((v >= (vsun+0.5e3)) or (v <= (vsun-0.5e3))): 
                Tini=Tini-0.001e6
                vc=sqrt(2*kb*Tini/mp)
                rc=mp*G*M/(4*kb*Tini)
                v=optimize.newton(parker_velocity,vsun,parker_velocity_p, args=(d,rc,vc), maxiter=50)
            T=Tini
    return(T)






def parker(pla:Planet,star:Star):
    """Compute the velocity and the density of the SW using the Parker model
        :param pla:
            The planet studied
        :type pla:
            Planet
        :param star:
            The associated star
        :type star:
            Star
    """

    kb=1.380658e-23 #J/K
    mp=1.660540210e-27 #kg
    dua=1.49597870700e11 #m
    G=6.6725985e-11 #N.m^2/kg^2
    M=star.mass
    R=star.radius
    d=pla.stardist
    t=star.age
    T=calc_temperature(M,t)
    vc=sqrt(2*kb*T/mp)
    rc=mp*G*M/(4*kb*T)
    vorb=sqrt(G*M/d)
    
    #Warning messages
    if (d<0.01*dua ) :
        print('Warning : The Parker model has not been verified for such star-planet distances')
    if (t <0.7e9):
        print('Warning : The Parker model is not precise for stars with t<0.7 Gyr')
    if (0.7e9 <= t <= 1.6e9):
        dlim=0.01*dua
    elif ( 1.6e9 < t < 3.5e9):
        dlim=0.02*dua
    elif (3.5e9 <= t ):
        dlim=0.03*dua
    else :
        dlim=0.
        print ('Star age is too small')

    if (d > dlim):
        try:
            v=optimize.newton(parker_velocity,350.0e3,parker_velocity_p, args=(d,rc,vc), maxiter=50)
        except (ZeroDivisionError, RuntimeError):
            v=0.
    else :
        try:
            v=optimize.newton(parker_velocity,10e3,parker_velocity_p, args=(d,rc,vc), maxiter=50)
        except (ZeroDivisionError, RuntimeError):
            v=0.
    
    veff=sqrt((v**2)+(vorb**2))
    Mls=masslossrate(t,R)
    n=Mls/(4*np.pi*(d**2)*v*mp)
    return(v,veff,n)