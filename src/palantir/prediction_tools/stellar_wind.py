#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Emilie Mauduit
"""

from scipy import optimize
import numpy as np
from math import sqrt, log, atan, sin
from importlib_resources import files
from scipy.io import readsav
from scipy.interpolate import RegularGridInterpolator

from palantir.prediction_tools import Planet, Star

import logging
logger = logging.getLogger('palantir.prediction_tools.stellar_wind')


# ============================================================ #
# ------------------------ ParkerGrid ------------------------ #
# ============================================================ #

class ParkerGrid:
    def __init__(self):
        maps_dir = files("palantir.scripts.input_files")
        parker_grid = readsav(maps_dir / "parker.sav")
        self.mass_ramp = np.array(parker_grid['m'],dtype='float')
        self.tcor_ramp = np.array(parker_grid['t'], dtype = 'float')
        self.dist_ramp = np.array(parker_grid['d'], dtype = 'float')
        self.velocity_grid = np.array(parker_grid['v'],dtype="float")
        self.interpolation_function_parker = RegularGridInterpolator((self.dist_ramp,self.tcor_ramp,self.mass_ramp),self.velocity_grid)


# ============================================================= #
# ------------------------ StellarWind ------------------------ #
# ============================================================= #


class StellarWind:
    def __init__(self, 
            Tcor : dict,
            d_alfven_point : dict,
        ):
        """Creates a stellar wind object.

        :param Tcor:
            Coronal temperature, in K.
        :type Tcor:
            float
        :param Bsw:
            Magnetic field of the stellar wind, in T.
        :type Bsw:
            float
        """

        self.density_planet = None
        self.density_star = None
        self.velocity_sw = None
        self.effective_velocity = None
        self.mass_loss_rate = None
        self.corona_temperature = Tcor
        self.distance_alfven_point = d_alfven_point
        self._alfven_velocity = None
        self.radial_mag_field = None
        self.azimuthal_mag_field = None
        self._total_mag_field = None
        self.perp_mag_field = None

    
    def __str__(self):
        return ("Electron density at planet: {} m-3 \n".format(self.density_planet)
        + "Electron density at star: {} m-3 \n".format(self.density_star)
        + "Effective velocity of the stellar wind : {} m.s-1 \n".format(self.effective_velocity)
        + "Temperature of the corona : {} MK \n".format(self.corona_temperature * 1e-6)
        + "Stellar wind magnetic field : {} T \n".format(self.perp_mag_field)
        + "Electron Alfven velocity : {} m.s-1 \n".format(self.alfven_velocity)
        + "Stellar mass-loss rate : {} kg.s-1\n".format(self.mass_loss_rate))


    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def corona_temperature(self):
        return self._corona_temperature
    @corona_temperature.setter
    def corona_temperature(self, value : dict):
        """Computes the value of the coronal temperature using X-ray measurements of the ROSAT catalog."""
        if np.isnan(value["Tcor"]) :
            star = value["star"]
            #T = StellarWind._calc_temperature_jmg2007(star.unnormalize_mass(), t) jmg 2007
            self._corona_temperature = StellarWind._calc_temperature(star.mass, star.rotperiod, star.Xray_flux)
        else :
            self._corona_temperature = value["Tcor"]

    @property
    def distance_alfven_point(self):
        return self._distance_alfven_point

    @distance_alfven_point.setter
    def distance_alfven_point(self, value:dict):
        """Finds the position of the Alfven point where Btot^2/2mu0 = n*mp*v^2"""
        parker_grid = value['parker_grid']
        d_values_coarse = np.logspace(np.log10(0.01),np.log10(100),60,endpoint=True)
        f_values_coarse = np.array([self._alfven_point_equation(
                                d,
                                value['star'], 
                                value['eccentricity'], 
                                self.corona_temperature,
                                dist_ramp=parker_grid.dist_ramp,
                                mass_ramp=parker_grid.mass_ramp,
                                tcor_ramp = parker_grid.tcor_ramp,
                                interpolation=parker_grid.interpolation_function_parker
                                ) 
                            for d in d_values_coarse],dtype='float')
        
        mask_coarse = ~np.isnan(f_values_coarse)
        idx_min_coarse = np.argmin(f_values_coarse[mask_coarse])

        if idx_min_coarse == 0 :
            d_min = d_values_coarse[mask_coarse][0] ; d_max = d_values_coarse[mask_coarse][20]
        elif idx_min_coarse == -1 :
            d_min = d_values_coarse[mask_coarse][-20], d_max = d_values_coarse[mask_coarse][-1]
        else :
            idx_min_temp = idx_min_coarse-10 if idx_min_coarse >= 10 else 0
            idx_max_temp = idx_min_coarse+10 if idx_min_coarse <= 49 else -1
            d_min = d_values_coarse[mask_coarse][idx_min_temp] ; d_max = d_values_coarse[mask_coarse][idx_max_temp]

        d_values = np.logspace(np.log10(d_min),np.log10(d_max),40, endpoint=True)
        f_values = np.array([self._alfven_point_equation(
                                d,
                                value['star'], 
                                value['eccentricity'], 
                                self.corona_temperature,
                                dist_ramp=parker_grid.dist_ramp,
                                mass_ramp=parker_grid.mass_ramp,
                                tcor_ramp = parker_grid.tcor_ramp,
                                interpolation=parker_grid.interpolation_function_parker
                                ) 
                            for d in d_values],dtype='float')
        mask = ~np.isnan(f_values)
        idx_min = np.argmin(f_values[mask])
        
        self._distance_alfven_point = d_values[mask][idx_min]

    @property
    def total_mag_field(self):
        if self._total_mag_field is None :
            self._total_mag_field = sqrt((self.radial_mag_field **2) + (self.azimuthal_mag_field**2))
        return self._total_mag_field

    @property
    def alfven_velocity(self):
        """ Compute the electronic Alfven velocity. """
        if self._alfven_velocity is None :
            mu0 = 4*np.pi*1e-7
            mp = 1.67e-27 #kg 
            c = 3e8 #m.s-1
            a = (c**2) * mu0 * self.density_planet * mp / (self.total_mag_field **2)
            self._alfven_velocity = np.sqrt((c**2)/(1 + a))
            #self.total_mag_field / np.sqrt(mu0 * self.density_planet * mp)
        return self._alfven_velocity

    ########## Methods ##########
    def compute_Parker_solution(self, planet : Planet, star : Star, parker_grid : ParkerGrid = None ):
        if parker_grid is not None:
            v, veff, ne_planet, ne_star, Mls = self._Parker_interpolation_grid(
                                                    star = star, 
                                                    distance = planet.stardist, 
                                                    T = self.corona_temperature,
                                                    dist_ramp=parker_grid.dist_ramp,
                                                    mass_ramp=parker_grid.mass_ramp,
                                                    tcor_ramp = parker_grid.tcor_ramp,
                                                    interpolation=parker_grid.interpolation_function_parker)
        else : 
            v, veff, ne_planet, ne_star, Mls = self._Parker(star=star, distance=planet.stardist, eccentricity=planet.eccentricity, T = self.corona_temperature)
        self.density_planet = ne_planet
        self.density_star = ne_star
        self.velocity_sw = v
        self.effective_velocity = veff
        self.mass_loss_rate = Mls 
    
    def compute_B_imf_components(self, planet: Planet, star: Star):
        self.radial_mag_field = self._calc_B_radial(
                d = planet.stardist,
                d_alfven_point=self._distance_alfven_point,
                R_star =star.unnormalize_radius(),
                magfield_surf=star.magfield
            )
        self.azimuthal_mag_field = self._calc_B_azimuthal(
                d = planet.stardist,
                T_star=star.rotperiod,
                vsw = self.velocity_sw,
                Br = self.radial_mag_field
            )
        self.perp_mag_field = self._calc_Bperp(
                d=planet.stardist,
                star=star,
                planet=planet,
                vsw = self.velocity_sw,
                d_alfven_point= self._distance_alfven_point
            )

    # --------------------------------------------------------- #
    # -------------------- Static methods --------------------- #

    @staticmethod
    def _calc_vsun(age: float):
        """Compute the velocity of the solar wind at 1 AU.
        :param age:
            stellar age [yr]
        :type age:
            float
        """
        tau = 2.56e7  # yr
        vo = 3971e3  # m/s
        vsun = vo / pow((1 + (age / tau)), 0.43)
        return vsun

    @staticmethod
    def _calc_nsun(age: float):
        """Compute the electron density of the solar wind at 1 AU.
        :param age:
            stellar age [yr]
        :type age:
            float
        """
        tau = 2.56e7  # yr
        no = 1.04e11  # m-3
        nsun = no / pow((1 + (age / tau)), 1.86)
        return nsun

    @staticmethod
    def _alfven_point_equation(
        d_alfven_point : float,
        star : Star,
        eccentricity : float,
        T_corona :float,
        dist_ramp : np.ndarray,
        mass_ramp : np.ndarray,
        tcor_ramp : np.ndarray,
        interpolation : RegularGridInterpolator
        ):

        mp = 1.660540210e-27  # kg
        mu0 = 4*np.pi * 1e-7
        try : 
            v, veff, ne_planet, ne_star, Mls = StellarWind._Parker_interpolation_grid(
                                                    star=star, 
                                                    distance=d_alfven_point, 
                                                    T = T_corona,
                                                    dist_ramp=dist_ramp,
                                                    mass_ramp=mass_ramp,
                                                    tcor_ramp = tcor_ramp,
                                                    interpolation=interpolation)

            Br = StellarWind._calc_B_radial(
                d = d_alfven_point, 
                d_alfven_point = d_alfven_point,
                R_star = star.unnormalize_radius(),
                magfield_surf= star.magfield)
            
            B_phi = StellarWind._calc_B_azimuthal(
                d = d_alfven_point,
                T_star=star.rotperiod,
                vsw = v,
                Br = Br
            )
            
            res1 = (Br**2 + B_phi**2) / mu0
            res2 = ne_planet * (v**2) * mp
            return np.abs((res1 - res2)/max(res1,res2))
        
        except ValueError : 
            return np.nan

    @staticmethod
    def _calc_B_radial(
        d : float,
        d_alfven_point : float,
        R_star : float,
        magfield_surf : float
        ):

        dua = 1.49597870700e11  # m
        #if d <= d_alfven_point :
        #    Br = magfield_surf * pow(R_star/(d * dua), 3)
        #else :
        #    Br_ini = magfield_surf * pow(R_star/(d_alfven_point * dua), 3)
        #    Br = Br_ini * pow(d_alfven_point/d, 2)
        
        dist = d * dua / R_star
        f = d_alfven_point * dua/ R_star
        Br = magfield_surf * (1 + ((f-1)/pow(dist, 1.5))) / (f * pow(dist,2))
        return Br
    
    @staticmethod
    def _calc_B_azimuthal(d : float,
        T_star : float,
        vsw : float,
        Br : float
    ) -> float:
        """Computes $B_{\phi}$, the azimuthal component of the interplanetary magnetic field, at the disance d from the star.
        
            :param d:
                The distance from the star where to evaluate $B_{\phi}$. In astronomical units.
            :type d:
                float
            :param T_star:
                The rotation rate of the star, in days.
            :param vsw:
                The stellar wind velocity, in the frame of the star, at the distance d. In $m.s^{-1}$.
            :type vsw:
                float
            :param Br:
                The radial component, $B_r$ of the interplanetary magnetic field, evaluated in d. In Tesla.
            :type Br:
                float

            :returns:
                $B_{\phi}(d)$ in Tesla.
            :rtype:
                float
        """
        dua = 1.49597870700e11  # m
        omega_star = 2 * np.pi / (T_star * 24 * 3600)

        Bp = Br * omega_star * d * dua / vsw
        return Bp
    
    @staticmethod
    def _calc_B_total(d : float,
        star : Star,
        vsw : float,
        d_alfven_point : float,
        ) -> float :
        """ Computes $B_{tot}$ (T) at given distance d (AU) from the star. It is computed as presented in Zarka et al, 2007, PSS.

            :param d:
                The distance from the star where to evaluate the magnetic field, in astronomical units.
            :type d:
                float
            :param star:
                The star to consider.
            :type star:
                :class:`~palantir.prediction_tools.star.Star`
            :param vsw:
                The stellar wind velocity, in the frame of the star, at the distance d. In $m.s^{-1}$.
            :type vsw:
                float
            :param d_alfven_point:
                The distance of the Alfvèn point, from the star, in astronomical units.
            :type d_alfven_point:
                float

            :returns:
                The perpendicular component of the interplanetary magnetic field, evaluated at the planet, in Tesla.
            :rtype:
                float
        """
        if np.isnan(star.magfield) :
            Psun = 25.5  # days
            Br0 = 2.6e-9 * Psun / star.rotperiod # T
            Bp0 = 2.4e-9 * Psun / star.rotperiod # T
            Br = Br0 * pow(d, -2)
            Bp = Bp0 * pow(d, -1)
        else :
            Br = StellarWind._calc_B_radial( 
                d = d, 
                d_alfven_point = d_alfven_point,
                R_star = star.unnormalize_radius(),
                magfield_surf = star.magfield)
            Bp = StellarWind._calc_B_azimuthal(
                d = d,
                T_star = star.rotperiod,
                vsw = vsw,
                Br = Br
            )
        Btot = sqrt(Br**2 + Bp**2) 
        return Btot

    @staticmethod
    def _calc_Bperp(
        d : float,
        star : Star,
        planet : Planet, 
        vsw : float,
        d_alfven_point : float,
        ) -> float :
        """ Compute the perpendicular component of the interplanetary magnetic field at the planet.
            It is computed as presented in Zarka et al, 2007, PSS.

            :param d:
                The distance from the star where to evaluate the magnetic field, in astronomical units.
            :type d:
                float
            :param star:
                The star to consider.
            :type star:
                :class:`~palantir.prediction_tools.star.Star`
            :param planet:
                The planet to consider.
            :type planet:
                :class:`~palantir.prediction_tools.planet.Planet`
            :param vsw:
                The stellar wind velocity, in the frame of the star, at the distance d. In $m.s^{-1}$.
            :type vsw:
                float
            :param d_alfven_point:
                The distance of the Alfvèn point, from the star, in astronomical units.
            :type d_alfven_point:
                float

            :returns:
                The perpendicular component of the interplanetary magnetic field, evaluated at the planet, in Tesla.
            :rtype:
                float
        """
        Psun = 25.5  # days
        dua = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2
        vorb = sqrt(
            G * star.unnormalize_mass() / (planet.stardist * dua)
        )

        if np.isnan(star.magfield) :
            Br0 = 2.6e-9 * Psun / star.rotperiod # T
            Bp0 = 2.4e-9 * Psun / star.rotperiod # T
            Br = Br0 * pow(d, -2)
            Bp = Bp0 * pow(d, -1)
        else :
            Br = StellarWind._calc_B_radial( 
                d = d, 
                d_alfven_point = d_alfven_point,
                R_star = star.unnormalize_radius(),
                magfield_surf = star.magfield)
            Bp = StellarWind._calc_B_azimuthal(
                d = d,
                T_star = star.rotperiod,
                vsw = vsw,
                Br = Br
            )
        alpha = atan(Bp / Br)
        beta = atan(vorb / vsw)

        Bperp = sqrt(Br**2 + Bp**2) * abs(sin(alpha - beta))
        return Bperp

    @staticmethod
    def _mass_lossrate(t: float, R: float, M:float, P:float) -> float:
        """Compute the stellar mass loss rate in kg.s-1
        :param t:
            stellar age [yr]
        :type t:
            float
        :param R:
            radius of the star [R_sun]
        :type R:
            float
        :param M:
            Stellar mass [M_sun]
        :type M:
            float
        :param P:
            Stellar rotation period [days]
        :type P:
            float

        """
        mp = 1.660540210e-27  # kg
        do = 1.49597870e11  # m
        MS = 1.989e30  # kg

        if (t >= 0.1e6) and (t <= 1e6):
            #Cranmer et al 2008 ApJ
            Mls = 1e-9 * MS / (365 * 86400)

        if (M > 0.25) and (M < 1.2) :
            if (t>1e6) and (t<=100e6) :
                # Bouvier & Gallet 2013 AA
                Mls = 1.8e-14 * np.power(25.5/P,1.58) * MS / (365 * 86400)
            elif (t > 100e6) and (t <= 5.2e9) : 
                #Johnstone et al 2015 AA (II)
                Mls = np.power(R,2) * np.power(25.5/P,1.33) * np.power(M,-3.36) * 2e-14 * MS / (365 * 86400)
            else :
                #Newkirk 1980
                vsun = StellarWind._calc_vsun(t)
                nsun = StellarWind._calc_nsun(t)
                Mls = 4 * np.pi * (do**2) * nsun * vsun * mp * (R**2)
        else :
            #Newkirk 1980
            vsun = StellarWind._calc_vsun(t)
            nsun = StellarWind._calc_nsun(t)
            Mls = 4 * np.pi * (do**2) * nsun * vsun * mp * (R**2)
        return Mls

    @staticmethod
    def _parker_velocity(v: float, d: float, rc: float, vc: float) -> float:
        """Expression of f(v(d))=0 based on Parker's model for SW.
        :param v:
            velocity [m/s]
        :type v:
            float
        :param d:
            distance between the star and the planet [m]
        :type d:
            float
        :param rc:
            critical radius above which the solar wind become super-sonic
        :type rc:
            float
        :param vc:
            critical velocity associated to rc
        :type vc:
            float
        """
        res = (
            ((v / vc) ** 2)
            - log((v / vc) ** 2)
            - log((d / rc) ** 4)
            - (4.0 / (d / rc))
            + 3.0
        )
        return res

    @staticmethod
    def _parker_velocity_derivative(v: float, d: float, rc: float, vc: float) -> float:
        """Expression of df(v(d))/dv based on Parker's model for SW.
        :param v:
            velocity [m/s]
        :type v:
            float
        :param d:
            distance between the star and the planet [m]
        :type d:
            float
        :param rc:
            critical radius above which the solar wind become super-sonic
        :type rc:
            float
        :param vc:
            critical velocity associated to rc
        :type vc:
            float
        """
        res = (2 * v / (vc**2)) - (2 / v)
        return res

    @staticmethod
    def _calc_temperature(M : float, P : float, Fx : float) -> float :
        """Compute the coronal temperature from 'Johnstone & Güdel, AA, 2015'.
            Gives Tcor in K.
        
        :param M:
            Stellar mass in Msun
        :type M:
            float

        :param P:
            Stellar rotation period.
        :type P:
            float

        :param Fx:
            X-ray flux of the star in erg.cm-2.s-1
        :type Fx:
            float
        """

        Psun = 25.5
        if not np.isnan(Fx) :
            Tcor = 0.11 * np.power(Fx,0.26)
        elif (not np.isnan(P)) and (not np.isnan(M)) :
            Tcor = 1.77 * np.power(M,-0.42) * np.power(Psun/P,0.52)
        elif not np.isnan(M) :
            Tcor = 1.77 * np.power(M,0.6)
        else :
            Tcor = 1.77
        return Tcor * 1e6

    @staticmethod
    def _calc_temperature_jmg2007(M: float, t: float) -> float:
        """Adjust the temperature of the corona for a star of given age.

        :param M:
            Mass of the star in kg
        :type M:
            float
        :param t:
            Age of the star in yr
        :type t:
            float

        """

        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        d = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2

        vsun = StellarWind._calc_vsun(t)
        # nsun = calc_nsun(t)

        Tini = (1e6 * 3.6e9) / (t + 1e9)
        vc = sqrt(2 * kb * Tini / mp)
        rc = mp * G * M / (4 * kb * Tini)
        v = optimize.newton(
            StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
        )
        if v <= vsun:
            while (v >= (vsun + 2e3)) or (v <= (vsun - 2e3)):
                Tini = Tini + 0.005e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
                )
            if v <= vsun:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini + 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
            else:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini - 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini

        else:
            while (v >= (vsun + 2e3)) or (v <= (vsun - 2e3)):
                Tini = Tini - 0.005e6
                vc = sqrt(2 * kb * Tini / mp)
                rc = mp * G * M / (4 * kb * Tini)
                v = optimize.newton(
                    StellarWind._parker_velocity, vsun, StellarWind._parker_velocity_derivative, args=(d, rc, vc), maxiter=50
                )
            if v <= vsun:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini + 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
            else:
                while (v >= (vsun + 0.5e3)) or (v <= (vsun - 0.5e3)):
                    Tini = Tini - 0.001e6
                    vc = sqrt(2 * kb * Tini / mp)
                    rc = mp * G * M / (4 * kb * Tini)
                    v = optimize.newton(
                        StellarWind._parker_velocity,
                        vsun,
                        StellarWind._parker_velocity_derivative,
                        args=(d, rc, vc),
                        maxiter=50,
                    )
                T = Tini
        return T

    @staticmethod
    def _Parker_interpolation_grid(
            star: Star, 
            distance: float, 
            T: float, 
            dist_ramp : np.ndarray,
            mass_ramp : np.ndarray,
            tcor_ramp : np.ndarray,
            interpolation : RegularGridInterpolator
            ) :
        """Compute the velocity and the density of the SW using a 3D-grid (M,Tcor,d_*p) of precomputed solutions of the Parker model.
        :param distance:
            The distance from the star where the parameters should be evaluated in astronomical units (AU).
        :type distance:
            float
        :param star:
            The associated star
        :type star:
            Star
        """
        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        dua = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2
        d = distance #AU
        t = star.age  # yr
        vc = sqrt(2 * kb * T / mp)
        rc = mp * G * star.unnormalize_mass() / (4 * kb * T)
        vorb = sqrt(G * star.unnormalize_mass() / (d * dua))


        if (distance * dua > np.max(dist_ramp) ) :
            dist = np.max(dist_ramp)
        elif (distance * dua < np.min(dist_ramp) ) :
            dist = np.min(dist_ramp)
        else :
            dist = d*dua
        
        if (star.unnormalize_mass() > np.max(mass_ramp)):
            mass = np.max(mass_ramp)
        elif (star.unnormalize_mass() < np.min(mass_ramp)):
            mass = np.min(mass_ramp)
        else :
            mass = star.unnormalize_mass()

        if (T > np.max(tcor_ramp)):
            tcor = np.max(tcor_ramp)
        elif (T < np.min(tcor_ramp)):
            tcor = np.min(tcor_ramp)
        else :
            tcor = T

        point = np.array([dist,tcor,mass])
        v = interpolation(point)[0]

        veff = sqrt((v**2) + (vorb**2))
        Mls = StellarWind._mass_lossrate(t, star.radius, star.mass, star.rotperiod)
        ne_planet = Mls / (4 * np.pi * ((d * dua) ** 2) * v * mp)
        ne_star = ne_planet * ((d * dua / star.unnormalize_radius()) ** 2)
        if (ne_planet <= 0) or (ne_star <= 0):
            logger.error("Negative stellar wind density is not physical.")
            raise ValueError("Negative or null stellar wind density is not physical.")
        return (v, veff, ne_planet, ne_star, Mls)

    @staticmethod
    def _Parker(star: Star, distance: float, eccentricity : float, T: float):
        """Compute the velocity and the density of the SW using the Parker model.
        :param distance:
            The distance from the star where the parameters should be evaluated in astronomical units (AU).
        :type distance:
            float
        :param star:
            The associated star
        :type star:
            Star
        """

        kb = 1.380658e-23  # J/K
        mp = 1.660540210e-27  # kg
        dua = 1.49597870700e11  # m
        G = 6.6725985e-11  # N.m^2/kg^2
        d = distance #AU
        t = star.age  # yr
        vc = sqrt(2 * kb * T / mp)
        rc = mp * G * star.unnormalize_mass() / (4 * kb * T)
        vorb = sqrt(G * star.unnormalize_mass() / (d * dua))
        # Warning messages
        if d < 0.01:
            logger.warning("Warning : The Parker model has not been verified for such star-planet distances")
            print(
                "Warning : The Parker model has not been verified for such star-planet distances"
            )
        if t < 0.7e9:
            logger.warning("Warning : The Parker model is not precise for stars with t<0.7 Gyr")
            print("Warning : The Parker model is not precise for stars with t<0.7 Gyr")

        if d >= 1.0 : #(d * dua / rc) >= 1.0:
            v= optimize.newton(
                StellarWind._parker_velocity,
                350e3,#1.1 * vc,
                StellarWind._parker_velocity_derivative,
                args=(d * dua, rc, vc),
                maxiter=50,
            )
                
        else:
            d_temp = 1.0 #0.9 * rc
            v_temp_ini = optimize.newton(
                StellarWind._parker_velocity,
                350.0e3,
                StellarWind._parker_velocity_derivative,
                args=(d_temp * dua, rc, vc),
                maxiter=50,
            )
            while abs(d_temp - d) >= 1e-5:
                d_temp = (9 * d_temp + d) / 10.0
                v_temp = optimize.newton(
                    StellarWind._parker_velocity,
                    v_temp_ini,
                    StellarWind._parker_velocity_derivative,
                    args=(d_temp * dua, rc, vc),
                    maxiter=50,
                )
                if (v_temp / vc > 1.0) and (d_temp / (rc / dua) <= 1.0):
                    v_temp = optimize.newton(
                        StellarWind._parker_velocity,
                        0.9 * v_temp_ini,
                        StellarWind._parker_velocity_derivative,
                        args=(d_temp * dua, rc, vc),
                        maxiter=50,
                    )
                v_temp_ini = v_temp
            v = v_temp

        veff = sqrt((v**2) + (vorb**2))
        Mls = StellarWind._mass_lossrate(t, star.radius, star.mass, star.rotperiod)
        ne_planet = Mls / (4 * np.pi * ((d * dua) ** 2) * v * mp)
        ne_star = ne_planet * ((d * dua / star.unnormalize_radius()) ** 2)
        if (ne_planet <= 0) or (ne_star <= 0):
            logger.error("Negative stellar wind density is not physical.")
            raise ValueError("Negative or null stellar wind density is not physical.")
        return (v, veff, ne_planet, ne_star,Mls)

    @staticmethod
    def _CME(star: Star, planet: Planet):
        """Evaluates if it is necessary to take into account CME in the stellar wind model.
        Returns a tuple containing the speed, the effective speed and the density of the stellar wind.

        :param star:
            Star of the system considered
        :type star:
            Star
        :param pla:
            Planet of the system considered
        :type pla:
            Planet
        """
        T = 2e6  # K
        dua = 1.49597870700e11  # m
        d = planet.stardist * dua
        G = 6.6725985e-11  # N.m^2/kg^2
        # nwo = 4.88e6
        nso = 7.1e6
        # nw = nwo * pow((d / dua), -2.31)
        ns = nso * pow((d / dua), -2.99)
        v = 5.26e5
        veff = sqrt((v**2) + (star.unnormalize_mass() * G / d))
        return (veff, ns, T)

