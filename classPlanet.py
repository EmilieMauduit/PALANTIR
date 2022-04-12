#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:53:23 2021

@author: emauduit
"""
import os
import scipy.io as sio
from scipy.io import readsav
import numpy as np
from math import *
from calc_tools import *


class Planet:
    
    def __init__(self,name,mp,rp,w,wo,L,d):
        self.name=name
        self.mass=mp
        self.radius=rp
        self.rotrate=w
        self.orbitperiod=wo
        self.luminosity=L
        self.stardist=d        
        
    def affiche(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " kg")
        print("Radius : ", self.radius, " m")
        print("Rotation rate : ", self.rotrate, " s-1")
        print("Orbital period : ", self.orbitperiod, "s-1")
        print("Luminosity : ", self.luminosity, "W")
        print("Distance to host star : ", self.stardist, "m")
        
    def affiche_norm(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " MJ")
        print("Radius : ", self.radius, " RJ")
        print("Rotation rate : ", self.rotrate, " wJ")

    
    #Methods
    
    def normalize_all(self,planorm):
        #Normalize the parameters to the values of the ones of the input planet (usually Jupiter)
        M=self.mass/planorm.mass
        R=self.radius/planorm.radius
        w=self.rotrate/planorm.rotrate
        return(M,R,w)
    def normalize_mass(self,planorm):
        M=self.mass/planorm.mass
        return(M)
    def normalize_radius(self, planorm):
        R=self.radius/planorm.radius
        return(R)
    def normalize_rotrate(self, planorm):
        w=self.rotrate/planorm.rotrate
        return(w)


    def calculate_luminosity(self,t,table1,table2,table3):
        """Retourne la liste des luminosités en fonction de l'âge selon la valeur de M"""
        MJ=1.8986e27
        M=self.mass/MJ
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

    def calculate_orbitalperiod(self,star):
        M=star.mass
        d=self.stardist
        G=6.6725985e-11 
        self.rotrate=pow(M*G/pow(d,3),1/2)

    
    def tidal_locking(self,t,Ms):
        """Calcul le période de rotation de la planète selon si l'on considère une rotation synchrone ou libre"""
        wJ=1.77e-4
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





    
