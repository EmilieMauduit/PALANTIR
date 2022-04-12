#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:02:05 2021

@author: emauduit
"""

import os
import scipy.io as sio
from scipy.io import readsav
import numpy as np
from math import *
from calc_tools import *


class Star:
    
    def __init__(self,name,M,R,t,Bsw,L):
        self.name=name
        self.mass=M
        self.radius=R
        self.age=t
        self.Bsw=1
        self.luminosity=L
        
    def affiche(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " kg")
        print("Radius : ", self.radius, " m")
        print("Age : ", self.age, " y")
        print("SW magnetic field : ", self.Bsw, 'T')
        print("Luminosity : ", self.luminosity, 'W')
        
    def affiche_norm(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " MS")
        print("Radius : ", self.radius, " RS")
        print("Age : ", self.age, " yS")
        print("SW magnetic field ", self.Bsw, 'BswS')
        print("Luminosity : ", self.luminosity, 'LS')

    def calculate_luminosity(self):
        a=0.39704170
        b=8.52762600
        c=0.00025546
        d=5.43288900
        e=5.56357900
        f=0.78866060
        g=0.00586685
        M=self.mass
        res1=a*pow(M,5.5)+b*pow(M,11)
        res2=c+pow(M,3)+d*pow(M,5)+e*pow(M,7)+f*pow(M,8)+g*pow(M,9.5)
        self.luminosity=res1/res2
    
    # Methods
    
    def normalize_all(self,starnorm):
        #Normalize the parameters to the values of the ones of the input star (usually the Sun)
        M=self.mass/starnorm.mass
        R=self.radius/starnorm.radius
        A=self.age/starnorm.age
        Bsw=self.Bsw/starnorm.Bsw
        L=self.luminosity/starnorm.luminosity
        return(M,R,A,Bsw,L)
    
    def normalize_mass(self, starnorm):
        M=self.mass/starnorm.mass
        return(M)
    def normalize_radius(self, starnorm):
        R=self.radius/starnorm.radius
        return(R)
    def normalize_age(self, starnorm):
        A=self.age/starnorm.age
        return(A)
    def normalize_Bsw(self, starnorm):
        Bsw=self.Bsw/starnorm.Bsw
        return(Bsw)
    def normalize_luminosity(self, starnorm):
        L=self.luminosity/starnorm.luminosity
        return(L)