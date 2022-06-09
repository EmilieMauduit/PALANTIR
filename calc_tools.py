#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: emauduit
"""


from array import array
import os
import scipy.io as sio
from scipy.io import readsav
import numpy as np

from math import *


def find_value(a:float,l:array, res:bool =False):
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
    min=np.abs(l[0]-a)
    indmin=0
    for i in range(1,len(l)):
        nmin=np.abs(l[i]-a)
        if (nmin <= min):
            min=nmin
            indmin=i
    if not res :
        return(indmin)
    else :
        return (l[indmin])

def synchro_dist(w:float,Qp:float,tsync:float,Mp:float,Rp:float,Ms:float):
    """ Compute the distance at which a planet should be so that its orbit is synchronised, depending on the age of the system.
        :param w:
            orbital period 
        :type w:
            float
        :param Qp:
            factor of dissipation by tidal effect 
        :type Qp:
            float
        :param tsync:
            time of synchronization
        :type tsync:
            float
        :param Mp:
            Mass of the planet [Mjup]
        :type Mp:
            float
        :param Rp:
            Radius of the planet [Rjup]
        :type Rp:
            float
        :param Ms:
            Mass of the star [Msun]
        :type Ms:
            float
    """
    G=6.6725985e-11 
    a=9*tsync/(4*w*0.26*Qp)
    b=G*Mp/pow(Rp,3)
    c=pow(Ms/Mp,2)
    res=Rp*pow(a*b*c,1/6)
    return(res)

def calc_vsun(t:float):
    """ Compute the velocity of the solar wind at 1 AU.
        :param t:
            stellar age [yr]
        :type t:
            float
    """
    tau=2.56e7 #yr
    vo=3397e3 #m/s
    vsun=vo/pow((1+(t/tau)),0.4)
    return(vsun)

def calc_nsun(t:float):
    """ Compute the electron density of the solar wind at 1 AU.
        :param t:
            stellar age [yr]
        :type t:
            float
    """
    tau=2.56e7 #yr
    no=1.6e10 #m-3
    nsun=no/pow((1+(t/tau)),1.5)
    return(nsun)

def magnetosphere_radius(Mm,n,veff,T):
    kb=1.380658e-23 #J/K
    mp=1.660540210e-27 #kg
    res1=(mp*n*(veff**2))+(2*n*kb*T)
    Rs=pow((np.pi*4e-7*(1.16**2)*(Mm**2))/(res1*8*(np.pi**2)),1/6)
    return(Rs)