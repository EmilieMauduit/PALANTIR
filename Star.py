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


# ============================================================= #
# ---------------------------- Star --------------------------- #
# ============================================================= #


class Star:
    """ Class to handle Star object.

    :param name:
        Name of the star used in the SIMBAD catalog.
    :type name:
        str
    :param M:
        Star mass in Solar masses
    :type M:
        float
    :param R:
        Star radius in Solar radius
    :type R:
        float
    :param t:
        Star age in yr
    :type t:
        float
    :param L:
        Star luminosity in Solar luminosity
    :type L:
        float
    """
    def __init__(self,name:str,M:float,R:float,t:float,s:float,B:float,L:float):

        """ Creates a Star object.
            :param name:
                Name of the star.
            :type name:
                str
            :param M:
                Star mass, in Sun masses.
            :type M:
                float
            :param R:
                Star radius, in Sun radiuses.
            :type R:
                float
            :param t:
                Star age, in yr.
            :type t:
                float
            :param s:
                Distance from Earth, in pc.
            :type s:
                float
            :param B:
                Star magnetic field, in T. Either known from litterature or computed.
            :type B:
                float
            :param L:
                Star luminosity, normalized to the Sun.
            :type L:
                float
        """

        self.name=name
        self.mass=M
        self.radius=R
        self.age=t
        self.obs_dist=s
        self.magfield=B
        self.luminosity=L
        
    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self):
        print("Name : ", self.name)
        print("Mass : ", self.mass, " MS")
        print("Radius : ", self.radius, " RS")
        print("Age : ", self.age*1e-9, " Gyr")
        print("Distance to Earth : ", self.obs_dist, " pc")
        print("SW magnetic field :  ", self.magfield, ' T')
        print("Luminosity : ", self.luminosity, 'LS')

    @property
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

    @property
    def calculate_B(self):     
        self.B=1
    
    # Methods

    def unnormalize_mass(self)->float:
        MS=1.989e30 #kg
        M=self.mass*MS
        return(M)
    def unnormalize_radius(self)->float:
        RS=6.96342e8 #m
        R=self.radius*RS
        return(R)
    def unnormalize_luminosity(self)->float:
        LS=3.826e26 #W
        L=self.luminosity*LS
        return(L)

    def alfven_radius(self,d:float)->float:
        """ Computes the Alfvén radius of the star.
            :param d:
                Distance between the star and the planet, in m.
            :type d:
                float
        """

        return(Ra)