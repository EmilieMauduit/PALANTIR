#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:53:23 2021

@author: emauduit
"""
import os
import Star
import scipy.io as sio
import pandas as pd
from scipy.io import readsav
import numpy as np
from math import *


from magneticmoment import *
from calc_tools import *


# ============================================================= #
# --------------------------- Planet -------------------------- #
# ============================================================= #



class Planet:
    
    def __init__(self,name:str,Mp:float,Rp:float,d:float,a:float,wo:float=None,w:float=None):
        """ Define a planet object. Every parameter must be normalized at Jupiter.
            :param name:
                Name of the planet
            :type name:
                str
            :param Mp:
                Mass of the planet, expected to be normalized to Jupiter's.
            :type Mp:
                float
            :param Rp:
                Radius of the planet, expected to be normalized to Jupiter's.
            :type Rp:
                float
            :param w:
                Rotation rate of the planet, expected to be normalized to Jupiter's.
            :type w:
                float
            :param wo:
                Orbital period of the planet, expected to be normalized to Jupiter's.
            :type wo:
                float
            :param d:
                Distance between the planet and its host star, in AU.
            :type d:
                float
            :param a:
                Semi-major axis, in AU.
            :type a:
                float
        """

        self.name=name
        self.mass=Mp
        self.radius=Rp
        self.orbitperiod=wo
        self.rotrate= tidal_locking
        self.stardist=d
        self.sm_axis=a    

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " Mj")
        print("Radius : ", self.radius, " Rj")
        print("Rotation rate : ", self.rotrate, " wj")
        print("Orbital period : ", self.orbitperiod, "woj")
        print("Distance to host star : ", self.stardist, " AU")
        print("Semi-major axis : ", self.sm_axis, " AU")
        
    def unnormalize_mass(self):
        MJ=1.8986e27 #kg
        M=self.mass*MJ
        return(M)
    def unnormalize_radius(self):
        RJ=71492e3 #m
        R=self.radius*RJ
        return(R)
    def unnormalize_rotrate(self):
        wJ=1.77e-4   #s-1
        w=self.rotrate*wJ
        return(w)


    def calculate_orbitalperiod(self,Ms:float):
        """ Computes the orbital period of the planet.
            :param Ms:
                Mass of the host star of the planet.
            :type Ms:
                float
        """
        d=self.stardist
        G=6.6725985e-11 
        self.orbitperiod=pow(Ms*G/pow(d,3),1/2)


    def tidal_locking(self,t:float,Ms:float):
        """Computes the rotation rate of the planet depending on a synchronized or free rotation. 
        Depending on the age of the host star, a limit distance is evaluated and if the star-planet distance is lower than that we consider a synchronized orbit, else we assume the rotation rate to be Jupiter's.

            :param t:
                Age of the host star.
            :type t:
                float
            :param Ms:
                Mass of the star, normalized to the Sun's.
            :type Ms:
                float
        """

        wJ=1.77e-4 #s-1
        d=self.stardist
        if (t>=1e10):
            Qpp=1e5
            dsync=synchro_dist(wJ,Qpp,1e10,self.mass,self.radius,Ms)
            if (d<=dsync):
                self.rotrate=self.orbitperiod
            else: 
                self.rotrate=1.0
        elif (1e8 <= t <= 1e10):
            Qppmax=1e6
            Qppmin=1e5
            dsync_min=synchro_dist(wJ,Qppmin,1e10,self.mass,self.radius,Ms)
            dsync_max=synchro_dist(wJ,Qppmax,1e8,self.mass,self.radius,Ms)
            ddsync=dsync_max-dsync_min
            if (d <= dsync_min+ddsync):
                self.rotrate=self.orbitperiod
            else :
                self.rotrate=1.0
        else :
            Qpp=1e6
            dsync=synchro_dist(wJ,Qpp,1e8,self.mass,self.radius,Ms)
            if (d<=dsync):
                self.rotrate=self.orbitperiod
            else: 
                self.rotrate=1.0

    def calculate_luminosity(self,t:float,table1:pd.DataFrame,table2:pd.DataFrame,table3:pd.DataFrame):
        """Retourne la liste des luminosités en fonction de l'âge selon la valeur de M"""
        MJ=1.8986e27
        M=self.mass
        print("Mass ds calc",M)
        if (0.0 <= M <= 5.0):
            dM=5.0-M
            if ( dM >= 2.5):
                L_ind=find_value(t,table1['t_Gyr'])
                L=table1['L_Ls'][L_ind]
            else :
                L_ind=find_value(t,table2['t_Gyr'])
                L=table2['L_Ls'][L_ind]
        elif (5.0 <= M <= 10.0) :
            dM=10.0-M
            if (dM >= 2.5):
                L_ind=find_value(t,table2['t_Gyr'])
                L=table2['L_Ls'][L_ind]
            else :
                L_ind=find_value(t,table3['t_Gyr'])
                L=table3['L_Ls'][L_ind]
        else :
            L_ind=find_value(t,table3['t_Gyr'])
            L=table3['L_Ls'][L_ind]
        self.luminosity=L
        return(L)

    def magnetosphere_radius(self,Mm:float,n:float,veff:float,T:float)->float:
        """Computes the radius of the magnetosphere in planetary radius.
        If lower than 1, the value is set to 1. A magnetosphere can not be inside the planet.
        
            :param Mm:
                Magnetic moment of the planet.
            :type Mm:
                float
            :param n:
                Electrons density of the stellar wind at the planet, in m-3.
            :type n:
                float
            :param veff:
                effective velocity of the electrons of the stellar wind at the planet, in m/s.
            :type veff:
                float
            :param T:
                Temperature of the corona of the host star, in K.
            :type T:
                float"""
        kb=1.380658e-23 #J/K
        mp=1.660540210e-27 #kg
        res1=(mp*n*(veff**2))+(2*n*kb*T)
        Rs=pow((np.pi*4e-7*(1.16**2)*(Mm**2))/(res1*8*(np.pi**2)),1/6)
        if ((Rs/self.radius) < 1.) :
            return(1.0)
        else:
            return(Rs)
    


        




    
