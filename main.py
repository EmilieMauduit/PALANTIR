#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:09:06 2021

@author: emauduit
"""

import os
import scipy.io as sio
from scipy.io import readsav
import numpy as np
from math import *
import pandas as pd

# Class
from Planet import *
from Star import *
from Target import *
from DynamoRegion import *
from MagneticMoment import *
from StellarWind import *

###Physical constant


MS=1.989e30 #kg
RS=6.96342e8 #m
AS=4.6e9 #yo
BSsw= 1 #T
LS=3.826e26 #W


MJ=1.8986e27 #kg
RJ=69911e3   #m
wJ=1.77e-4   #s-1

ME=5.97237e24 #kg
RE=6371.0e3   #m
wE=7.27e-5  #s-1

rc=0.85*RJ
rhoc=1800 #kg/m3


# Configuration settings input

config=pd.read_csv(r"/Users/emauduit/Documents/Thèse/Sélection des cibles/Programmes/parametres.csv", delimiter=';')


# Data input


# Main

sol=Star('Soleil', MS,RS,AS,BSsw,LS)
#sol.affiche()
jup=Planet('Jupiter',MJ,RJ,wJ)
#jup.affiche()


# Calcul du rayon critique de la région dynamo

rc1,rc2=LaneEmden(jup.mass,jup.radius)

rhoc1=rhodyn(jup.mass, jup.radius, rc1)
rhoc2=rhodyn(jup.mass, jup.radius, rc2)

print("Rho dynamo : rhoc1= ",rhoc1," ,rhoc2= " ,rhoc2)


# Calcul du moment magnétique

if config.Value[0]==1:
    M1_bla=blackett(rc1, wJ, rhoc1)
    M2_bla=blackett(rc2, wJ, rhoc2)
    print("Blackett : M1= ", M1_bla, " T, M2= " ,M2_bla, ' T')
else :
    M1_bla=0
    M2_bla=0
    
if config.Value[1]==1:
    M1_curt=curtis(rc1, wJ, rhoc1, E)
    M2_curt=curtis(rc2,wJ, rhoc2, E)
    print("Curtis : M1= ", M1_curt, " T, M2= " ,M2_curt, ' T')
else :
    M1_curt=0
    M2_curt=0

if config.Value[2]==1:
    M_RC=Reiners_Christensen(sol, jup, sol, jup)
    print("Reiners-Christensen : M=", M_RC, " T" )
else :
    M_RC=0

if config.Value[3]==1:
    M1_sa=sano(rc1,wJ,rhoc1)
    M2_sa=sano(rc2,wJ,rhoc2)
    print("Sano : M1= ", M1_sa, " T, M2= " ,M2_sa, ' T')
else :
    M1_sa=0
    M2_sa=0



# print(rc1/jup.radius, rc2/jup.radius)

# earth=Planet('Earth', ME, RE, wE)
# earth.affiche()
