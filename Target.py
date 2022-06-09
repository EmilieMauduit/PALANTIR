#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:02:05 2021

@author: emauduit
"""

import os
import numpy as np
from math import *

from calc_tools import *


# ============================================================= #
# --------------------------- Target -------------------------- #
# ============================================================= #


class Target:

    def __init__(self, name:str,Bmax:float,fmax:float,Pem:float,Prec:float,F:float,fmax_star:float):
        """ Creates a Target object.
            :param name:
                Name of the target.
            :type name:
                str
            :param Bmax:
                Maximum magnetic field estimated, in T.
            :type Bmax:
                float
            :param fmax:
                Maximum frequency of the emission at the planet, in Hz.
            :type fmax:
                float
            :param Pem:
                Power of the emission, in W.
            :type Pem:
                float
            :param Prec:
                Power received by the planet from the star, in W.
            :type Prec:
                float
            :param F:
                Flux of the emission at Earth, in W.m-2
            :type F:
                float
            :param fmax_star:
                Maximum frequency of the emission at the star, in Hz.
            :type fmax_star:
                float
        """

        self.name=name
        self.mag_field=Bmax
        self.freq_max=fmax
        self.pow_emission=Pem
        self.pow_received=Prec
        self.flux=F
        self.freq_max_star=fmax_star

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    def talk(self):
        print("Name of the system : ", self.name)
        print("Maximum magnetic field : ", self.mag_field, " T")
        print("Maximum frequency emission at the planet : ", self.freq_max*1e-6, " MHz")
        print("Power of the emission at the star : ", self.pow_emission, " W")
        print("Power of the emission received at the planet : ", self.pow_received, " W")
        print("Flux of the emission received by the instrument : ", self.flux, " W.m-2")
        print("Maximum frequency of the emission at th star : ", self.freq_max_star*1e-6, " MHz")
    
    
    def select_target(self)->bool:
        """Select or not a target according to various criterions set by the user.
        """
        if():
            return(True)
        else:
            return(False)
