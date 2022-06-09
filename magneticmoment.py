#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:46:05 2021

@author: emauduit
"""

import numpy as np
from math import *

from Planet import *
from Star import *

from calc_tools import *



#### Modèles du moment dipolaire magnétique

#Scale laws

def blackett(rc,w,rhoc):
    """Modèle de Blackett, 1947"""
    return(rhoc*w*(rc**5))

def Busse(rc,w,rhoc):
    """Modèle de Busse, 1976"""
    res=pow(rhoc,1/2)*w*pow(rc,4)
    return(res)

def Mizu_moderate(rc,w,rhoc,sig):
    """Modèle de Mizutani avec convection modérée"""
    res=pow(rhoc,1/2)*pow(w,3/4)*pow(rc,7/2)*pow(sig,-1/4)
    return(res)

def Mizu_slow(rc,w,rhoc,sig):
    """Modèle de Mizutani avec convection lente"""
    res=pow(rhoc,1/2)*pow(w,1/2)*pow(rc,3)*pow(sig,-1/4)
    return(res)

def curtis(rc,w,rhoc,E):
    """Modèle de Curtis, 1986"""
    a=pow(rhoc, 1/3)*pow(w,1/2)*pow(E,1/6)
    b=pow(rc,7/2)
    return(a*b)

def sano(rc,w,rhoc):
    """Modèle de Sano, 1993"""
    a=pow(rhoc,1/2)*pow(rc,7/2)*w
    return(a)

#Simulations
    

def Reiners_Christensen(planet,star,jup,sol,table1,table2,table3, normalize=0):
    """Modèle de Reiners-Christensen, 2010. Si normalize n'est pas précisé on ne normalise pas à Jupiter et au Soleil par défaut. Pour normaliser : normalize=1 
         Retourne le moment magnétique en A.m2 avec un champ B en Gauss."""
    #if normalize==0:
    Mp=planet.mass
    Rp=planet.radius
    t=star.age
    #-0.5e9
    LJ=jup.luminosity
    if( planet.name == 'Jupiter') :
        L=1.0
        print('jup')
    else :
        L=planet.luminosity/LJ
        #L=planet.calculate_luminosity(t,table1,table2,table3)/LJ
    print('dans RC : L=',L)
    #else :
    #    Mp=planet.normalize_mass(jup)
    #    Rp=planet.normalize_radius(jup)
    #    t=star.age
    #    L=planet.calculate_luminosity(t,table1,table2,table3)/jup.luminosity
    
    a=Mp*(L**2)*pow(Rp,11)
    b=pow(1-(0.17/Mp),3)/pow((1-0.17),3)
    res=b*pow(a,1./6)
    return(res)


def magnetic_moment(self, model:str, rc:float, rhoc:float, star:Star):
    """ Compute the value of the magnetic moment of the planet using the specified model
        
    :param model:
        Name of the model : blackett, busse, mizu_mod, mizu_slow, sano, rein-chris
    :type model:
        str
    :param rc:
        Radius of the dynamo region [m]
    :type rc:
        float
    :param rhoc:
        Mass density in the dynamo region [m-3]
    :type rhoc:
        float
    """
    if (model=='blackett'):
        Mm=blackett(rc, self.rotrate, rhoc)
        return(Mm)
    if (model=='busse'):
        Mm=Busse(rc, self.rotrate, rhoc)
        return(Mm)
    if (model=='mizu_mod'):
        Mm=Mizu_moderate(rc,self.rotrate,rhoc,1.0)
        return(Mm)
    if (model=='mizu_slow'):
        Mm=Mizu_slow(rc, self.rotrate,rhoc,1.0)
    if (model=='sano') :
        Mm=sano(rc, self.rotrate,rhoc)
        return(Mm)
    if (model=='rein-chris'):
        table_1MJ=pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/1MJ.csv', delimiter=';')
        table_5MJ=pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/5MJ.csv', delimiter=';')
        table_10MJ=pd.read_csv(r'/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/10MJ.csv', delimiter=';')
        Mm=Reiners_Christensen(self,star,table_1MJ,table_5MJ,table_10MJ)
        return(Mm)


# ============================================================= #
# ----------------------- MagneticMoment ---------------------- #
# ============================================================= #

