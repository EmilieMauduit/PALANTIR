#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:02:05 2021

@author: Emilie Mauduit
"""

from math import pow
import numpy as np
from typing import List
from astroquery.simbad import Simbad

import logging
log = logging.getLogger('palantir.prediction_tools.star')

# --------------------------------------------------------- #
# ------------ Useful functions for the class ------------- #
# --------------------------------------------------------- #


def TOUT(mass: float):
    "From Tout et al, 1996"
    theta2 = 1.71535900
    iota2 = 6.59778800
    kappa2 = 10.08855000
    lambda2 = 1.01249500
    mu2 = 0.07490166
    nu2 = 0.01077422
    eta2 = 3.08223400
    omega2 = 17.84778000
    pi2 = 0.00022582
    a = (
        theta2 * pow(mass, 2.5)
        + iota2 * pow(mass, 6.5)
        + kappa2 * pow(mass, 11)
        + lambda2 * pow(mass, 19)
        + mu2 * pow(mass, 19.5)
    )
    b = (
        nu2
        + eta2 * pow(mass, 2)
        + omega2 * pow(mass, 8.5)
        + pow(mass, 18.5)
        + pi2 * pow(mass, 19.5)
    )
    return a / b


# ============================================================= #
# ---------------------------- Star --------------------------- #
# ============================================================= #


class Star:
    def __init__(self, 
        name: str,
        mass: float, 
        radius: dict, 
        age: float, 
        obs_dist: float, 
        sp_type : str, 
    ):

        """Creates a Star object.
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
            Star age, in Gyr.
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
        """

        self.name = name
        self.main_id = name
        self.mass = mass
        self.radius = radius
        self.age = age * 1e9
        self.luminosity = mass
        self.obs_dist = obs_dist
        self.sp_type = sp_type
        self._sp_type_code = None
        self._rotperiod = None
        self.effective_temperature = None
        self.magfield = None

    
    def __str__(self):
        return("Star name : " + self.name + "\n"
            +"Star mass : M* = {} MS\n".format(self.mass)
            +"Star radius : R* = {} RS\n".format(self.radius)
            +"Star age : t* = {} Gyr\n".format(self.age * 1e-9)
            +"Star distance to Earth : s = {} pc\n".format(self.obs_dist)
            +"Star rotational period : wrot* = {} days\n".format(self.rotperiod)
            +"Star magnetic field :  B* = {} T\n".format(self.magfield)
            +"Star luminosity : L* = {} LS\n".format(self.luminosity)
            +"Spectral type : " + self.sp_type
        )

    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def main_id(self):
        return self._main_id
    
    @main_id.setter
    def main_id(self, value:str):
        if value  == 'Sun' :
            self._main_id = value
        else :
            query_id = Simbad.query_object(value)
            if query_id is not None :
                self._main_id = query_id['MAIN_ID'][0]
            else :
                self._main_id = value


    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value: dict):
        models = value["models"]
        radius = value["radius"]
        if np.isnan(radius):
            self._radius = self._calculate_radius(models, self.mass)
        else:
            self._radius = radius

    @property
    def rotperiod(self):
        self._rotperiod = self._compute_rotperiod(age=self.age)
        return self._rotperiod
    
    @property
    def luminosity(self):
        return self._luminosity
    
    @luminosity.setter
    def luminosity(self, value : float):
        a = 0.39704170
        b = 8.52762600
        c = 0.00025546
        d = 5.43288900
        e = 5.56357900
        f = 0.78866060
        g = 0.00586685
        M = value
        res1 = a * pow(M, 5.5) + b * pow(M, 11)
        res2 = (
            c
            + pow(M, 3)
            + d * pow(M, 5)
            + e * pow(M, 7)
            + f * pow(M, 8)
            + g * pow(M, 9.5)
        )
        self._luminosity = res1 / res2

    @property
    def sp_type_code(self):
        if self._sp_type_code is None :
            self._sp_type_code = self._decode_sp_type(self.sp_type)
        return self._sp_type_code

    def unnormalize_mass(self) -> float:
        MS = 1.989e30  # kg
        M = self.mass * MS
        return M

    def unnormalize_radius(self) -> float:
        RS = 6.96342e8  # m
        R = self.radius * RS
        return R

    def unnormalize_luminosity(self) -> float:
        LS = 3.826e26  # W
        L = self._luminosity * LS
        return L

    def obs_dist_meters(self) -> float:
        pc = 3.08568e16  # m
        return self.obs_dist * pc

    def alfven_radius(self, d: float) -> float:
        """Computes the Alfvén radius of the star.
        :param d:
            Distance between the star and the planet, in m.
        :type d:
            float
        """

        Ra = 1

        return Ra

    def compute_effective_temperature(self, value):
        """ Effective temperature of the star, computed with the Stefan-Boltzman law."""
        Teff = value
        LS = 3.826e26  # W
        sigma_sb = 5.670374e-8 # W.m-2.K-4
        if self.effective_temperature is None :
            if np.isnan(Teff) :
                self.effective_temperature = np.power(self.luminosity*LS/(4*np.pi*sigma_sb* (self.unnormalize_radius()**2)),1/4)
            else :
                self.effective_temperature = Teff


    def compute_magnetic_field(self, value :dict):
        model = ['Bstar_original'] if (self.mass > 1.6) else value['model']
        mag_field = value['mag_field']
        if self.magfield is None :
            if np.isnan(mag_field):
                if len(model) > 1 :
                    log.error('ValueError:Two models were selected for computing stellar magnetic field, only one can be chosen.')
                    raise ValueError('Two models were selected for computing stellar magnetic field, only one can be chosen.')
                else :
                    self.magfield = self._compute_magfield(model = model[0], 
                                        rotperiod = self.rotperiod, 
                                        Teff=self.effective_temperature,
                                        mass = self.mass,
                                        sp_type= self.sp_type
                                        ) * 1e-4
            else : 
                self.magfield = mag_field * 1e-4

    @staticmethod
    def _calculate_radius(
        models: List[str], mass: float, Rmean: bool = True, Rmax: bool = False
    ):
        R = []
        for model in models:
            if "Tout" in model:
                R.append(TOUT(mass))
        if Rmean:
            return np.mean(R)
        if Rmax:
            return np.max(R)

    @staticmethod
    def _compute_rotperiod(age: float):
        """Define the rotational period of the star in days."""
        tau = 2.56e7  # yr
        K = 0.6709505359255223
        rotperiod = K * pow(1 + (age / tau), 0.7)
        return rotperiod

    @staticmethod
    def _compute_magfield(model : str, rotperiod: float, Teff : float, mass : float, sp_type :str):

        if model == 'Bstar_original' :
            Psun = 25.5  # days
            Bsun = 1.435e-4  # T
            magfield = Bsun * Psun / rotperiod

        elif model == 'Bstar_polyfit' :
            sptype_dict = {'D': 374, 'M':120.4,'K':6.53,'G':1.32,'F':2.68}
            coeffs_polyfit = { 'B1' : [5.9252252,-46.444525,218.10474,-481.21961,520.51852,-271.33952,54.7305532],
                                'B4' : [290242.68,-321299.4,133345.3,-24588.728,1699.7549],
                                'B6' : [5.707718,-28.489377,120.03282,-247.65387,255.10662,-128.3299,25.233201]}
            B1=0 ; B4 = 0 ; B6 = 0
            for i,coeff in enumerate(coeffs_polyfit['B1']):
                B1 += coeff * (mass**i)
            B1 = 10**B1

            T = np.log10(Teff)
            for i,coeff in enumerate(coeffs_polyfit['B4']):
                B4 += coeff * (T**i)
            B4 = 10**B4

            for i,coeff in enumerate(coeffs_polyfit['B6']):
                B6 += coeff * (mass**i)

            B6 = (10**B6) / rotperiod

            if ((Teff > 7000) or (Teff < 2000)) and (mass > 1.6) :
                for i,sptype in enumerate(sptype_dict) :
                    for letter,B in sptype_dict.items():
                        if letter in sptype :
                            magfield = B
            elif ((Teff > 7000) or (Teff < 2000)) and not (mass > 1.6) :
                magfield = B6
            elif (mass > 1.6) and not ((Teff > 7000) or (Teff < 2000)) :
                magfield = B4
            else :
                magfield = np.sqrt(B4*B6)

        return magfield

    @staticmethod
    def _decode_sp_type(sp_type:str):
        types_roman = ["III", "IV", "VI", "VII", "II", "I"]
        sp_type_code = 3

        if sp_type == 'nan' :
            sp_type_code = 2 
            return sp_type_code 
        else :
            sp_type = 'G0V+pul' if (sp_type == 'G0Vpul') else sp_type
            sp_type = 'M3.5V' if (sp_type == 'M(3.5+/-0.5) V') else sp_type

        #Special type that we do not consider
        special_type = ['AM Her', 'AM', 'Am', 'Catac. var.']

        for spe in special_type:
            if spe in sp_type :
                sp_type_code = 3
                return sp_type_code
        
        if ('wd' in sp_type) or ('psr' in sp_type) or ('pul' in sp_type) :
            sp_type_code = 3
            return sp_type_code

        if ('L' in sp_type) or ('T' in sp_type) or ('Y' in sp_type) or ('W' in sp_type) :
            sp_type_code=3
            return sp_type_code
        
        if ('V' in sp_type) and ('IV' not in sp_type) and ('VI' not in sp_type):
            sp_type_code = 1
            return sp_type_code
        elif ('IV-V' in sp_type) or ('IV/V' in sp_type) :
            sp_type_code = 1
            return sp_type_code
        elif ("III" not in sp_type) and ("IV" not in sp_type) and ("VI" not in sp_type) and ("VII" not in sp_type) and ("II" not in sp_type) and ("I"not in sp_type) :
            if ('D' in sp_type) or ('M' in sp_type) or ('K' in sp_type) or ('G' in sp_type) or ('F' in sp_type) :
                sp_type_code=1
                return sp_type_code
        return sp_type_code
        
    
    @staticmethod
    def _decode_sp_type_old(sp_type : str):
        types_roman = ["III", "IV", "VI", "VII", "II", "I"]

        if sp_type == 'nan' :
            sp_type_code = 2 
            return sp_type_code 
        else :
            sp_type = 'G0V+pul' if (sp_type == 'G0Vpul') else sp_type
            sp_type = 'M3.5V' if (sp_type == 'M(3.5+/-0.5) V') else sp_type

        #Special type that we do not consider
            special_type = ['AM Her', 'AM', 'Am', 'Catac. var.']

            for spe in special_type:
                if spe in sp_type :
                    sp_type_code = 3
                    return sp_type_code 
                    
            #Checking for binary systems
            sp_type = sp_type.split('+')
            if len(sp_type) <= 1 :
                sp_type = sp_type[0].split('or')		

            sp_type_code = []
            for sp in sp_type :
                code = 0
                #Removing unnecessary spaces
                names = sp.split(' ')
                sp=''
                for i in range (len(names)):
                    sp += names[i]

                #Discard pulsars
                if ('wd' in sp.lower()) or ('psr' in sp.lower()) or ('pul' in sp.lower()) :
                    code = 3
                    sp_type_code.append(code)
                    continue

                #Check for post-fixes
                postfixes_1 = ['A','B','C','D','E','p','e','m']
                for pf in postfixes_1 :
                    if sp[-1] == pf :
                        sp = sp[0:-1]
                        print('pf1')
                        break
                
                #postfixes_2 = ['ep','pe']
                #for pf in postfixes_2 :
                #	if sp[-2:] == pf :
                #		sp = sp[0:-2]
                #		print('pf2')
                #		break

                sp = sp[0:-5] if (sp[-5:] == 'pecul') else sp
                sp = sp[0:-6] if (sp[-6:] == 'pecul.') else sp

                #Case of dwarves

                sp = sp[0:-6]+'V' if (sp[-6:] == '-dwarf') else sp
                sp = sp[0:-5]+'V' if (sp[-5:] == 'dwarf') else sp

                #First check for types in roman numbers

                for rom in types_roman :
                    if rom in sp :
                        code = 3
                        
                if code != 0 :
                    sp_type_code.append(3)
                    continue

                #Prefixes 'd' and 'sd'
                sp = sp[1:]+'VII' if (sp[0:1].lower() == 'd') else sp
                sp = sp[2:]+'VI' if (sp[0:2].lower() == 'sd') else sp

                #Postfixes 'Va','Vb','Ia', 'Ib'
                postfixes_2 = ['Va','Vb','Ia', 'Ib']
                for pf in postfixes_2 :
                    if sp[-2:] == pf :
                        sp = sp[0:-1]
                        break

                postfixes_3 = ['Vab','Iab']
                for pf in postfixes_3 :
                    if sp[-3:] == pf :
                        sp = sp[0:-2]
                        break

                if sp == '' :
                    code = 2 
                if code != 0 :
                    sp_type_code.append(code)
                    continue

                code = 3 if (sp[-1:] == 'I' or sp[-2:] == 'IV') else 0
                if code != 0 :
                    sp_type_code.append(code)
                    continue

                code = 1 if (sp[-1:] == 'V') else 0
                if code != 0 :
                    sp_type_code.append(code)
                    continue

                if code == 0 :
                    code = 2
                    sp_type_code.append(code)
                    continue
            sp_type_code = np.array(sp_type_code)

            return np.min(sp_type_code)			
