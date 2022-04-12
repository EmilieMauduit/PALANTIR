#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:36:08 2021

@author: emauduit
"""

import numpy as np
from math import *
from scipy import optimize
import math as m
from calc_tools import *



# Modèles de calcul du rayon de la région dynamo

def rlin(a,Rp):
    return(a*Rp)

def rhoLE(r,Mp, Rp, rhot):
    a=np.pi*Mp/(4*pow(Rp, 3))
    b=np.pi*r/Rp
    res=a*np.sin(b)/b
    return(res-rhot)

def rhoLEp(r,Mp,Rp,rhot):
    a=np.pi*Mp/(4*pow(Rp, 3))
    res=a*((np.cos(np.pi*r/Rp)/r)-(np.sin(np.pi*r/Rp)*Rp/(np.pi*pow(r,2))))
    return(res)

def LaneEmden(Mp, Rp,rhot):
    try:
        res=optimize.newton(rhoLE,Rp/2, fprime=rhoLEp , args=(Mp, Rp, rhot), maxiter=50)
    except(ZeroDivisionError, RuntimeError):
        res=0.0
    return(res)

# Modèles de calcul de la densité de masse dans la région dynamo

def rhodyn(Mp, Rp,rc):
    rhom=3*Mp/(4*np.pi*pow(Rp,3))
    a=rc/Rp
    if a==0.0 :
        return(0.0)
    else:
        res=rhom*pow(a,-3)*((m.sin(a*np.pi)/np.pi) - a*m.cos(a*np.pi))
    return(res)
