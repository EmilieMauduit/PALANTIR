#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Emilie Mauduit
"""

from scipy import optimize
import numpy as np
from math import sqrt, log, atan, sin

from palantir.prediction_tools.planet import Planet
from palantir.prediction_tools.star import Star

import logging
logger = logging.getLogger('palantir.prediction_tools.stellar_wind_old')
# ============================================================= #
# ------------------------ StellarWind ------------------------ #
# ============================================================= #


class StellarWindOld:
    def __init__(self, 
            ne_planet: float, 
            ne_star : float, 
            v : float, 
            ve: float, 
            Bsw: dict,
            Tcor: float = None, 
        ):
        """Creates a stellar wind object.

        :param ne_planet:
            Electron density in the stellar wind, in m-3. At the distance of the planet from the star.
        :type ne_planet:
            float
        :param ne_star:
            Electron density in the stellar wind, in m-3. At the surface of the star.
        :type ne_star:
            float
        :param ve:
            Effective velocity of the stellar wind, in m.s-1.
        :type ve:
            float
        :param Tcor:
            Coronal temperature, in K.
        :type Tcor:
            float
        :param Bsw:
            Magnetic field of the stellar wind, in T.
        :type Bsw:
            float
        """

        self.density_planet = ne_planet
        self.density_star = ne_star
        self.velocity_sw = v
        self.effective_velocity = ve
        self._corona_temperature = Tcor
        self.distance_alfven_point = Bsw
        self.radial_mag_field = Bsw
        self.azimuthal_mag_field = Bsw
        self._total_mag_field = None
        self.perp_mag_field = Bsw
        self._alfven_velocity = None
        
    
    def __str__(self):
        return ("Electron density : {} m-3 \n".format(self.density)
        + "Effective velocity of the stellar wind : {} m.s-1 \n".format(self.effective_velocity)
        + "Temperature of the corona : {} MK \n".format(self.corona_temperature * 1e-6)
        + "Stellar wind magnetic field : {} T \n".format(self.perp_mag_field)
        + "Electron Alfven velocity : {} m.s-1 \n".format(self.alfven_velocity))


    # --------------------------------------------------------- #
    # ------------------------ Methods ------------------------ #

    @property
    def corona_temperature(self):
        return self._corona_temperature
    @corona_temperature.setter
    def corona_temperature(self, value :  Star):
        """Computes the value of the coronal temperature using X-ray measurements of the ROSAT catalog."""
        star = value["star"]
        #T = StellarWind._calc_temperature_jmg2007(star.unnormalize_mass(), t) jmg 2007
        self._corona_temperature = StellarWind._calc_temperature(star.mass, star.rotperiod, star.Xray_flux)

    @property
    def distance_alfven_point(self):
        return self._distance_alfven_point

    @distance_alfven_point.setter
    def distance_alfven_point(self, value:dict):
        """Finds the position of the Alfven point where Btot^2/2mu0 = n*mp*v^2"""

        d_values_coarse = np.logspace(np.log10(0.01),np.log10(100),60,endpoint=True)
        f_values_coarse = np.array([self._alfven_point_equation(d,value['star'], value['eccentricity']) for d in d_values_coarse],dtype='float')
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
        f_values = np.array([self._alfven_point_equation(d,value['star'], value['eccentricity']) for d in d_values],dtype='float')
        mask = ~np.isnan(f_values)
        idx_min = np.argmin(f_values[mask])
        
        self._distance_alfven_point = d_values[mask][idx_min]

    @property
    def radial_mag_field(self):
        return self._radial_mag_field
    @radial_mag_field.setter
    def radial_mag_field(self, value: dict) :
        if ("star" or "planet") not in value.keys() :
            logger.error("KeyError : Star or Planet not in value.")
            raise KeyError("Star or Planet not in value.")
        
        Br = self._calc_B_radial(
            d = value['planet'].stardist,
            d_alfven_point=self._distance_alfven_point,
            R_star = value['star'].unnormalize_radius(),
            magfield_surf=value['star'].magfield
        )
        self._radial_mag_field = Br
    
    @property
    def azimuthal_mag_field(self):
        return self._azimuthal_mag_field
    @azimuthal_mag_field.setter
    def azimuthal_mag_field(self, value : dict) :
        if ("star" or "planet" or "vsw") not in value.keys() :
            logger.error("KeyError : Star or Planet or vsw not in value.")
            raise KeyError("Star or Planet or vsw not in value.")
        
        B_phi = self._calc_B_azimuthal(
            d = value['planet'].stardist,
            T_star=value['star'].rotperiod,
            vsw = value['vsw'],
            Br = self._radial_mag_field
        )
        self._azimuthal_mag_field = B_phi
    
    @property
    def total_mag_field(self):
        if self._total_mag_field is None :
            self._total_mag_field = sqrt((self._radial_mag_field **2) + (self._azimuthal_mag_field**2))
        return self._total_mag_field

    @property
    def perp_mag_field(self):
        return self._perp_mag_field

    @perp_mag_field.setter
    def perp_mag_field(self, value: dict):
        if ("planet" or "star" or "vsw") not in value:
            logger.error("KeyError : Planet or Star or Vsw not in value.")
            raise KeyError("Planet or Star or Vsw not in value.")
        
        Bperp = self._calc_Bperp(
            d=value['planet'].stardist,
            star=value['star'],
            planet=value['planet'],
            vsw = value['vsw'],
            d_alfven_point= self._distance_alfven_point
            )
        self._perp_mag_field = Bperp

    @property
    def alfven_velocity(self):
        """ Compute the electronic Alfven velocity. """
        if self._alfven_velocity is None :
            mu0 = 4*np.pi*1e-7
            mp = 1.67e-27 #kg 
            self._alfven_velocity = self.total_mag_field / np.sqrt(mu0 * self.density * mp)
        return self._alfven_velocity

    @classmethod
    def from_system(cls, star: Star, planet: Planet, T_cor : float = None):
        """Creates a StellarWind object corresponding to the given Star-Planet system.
        :param star:
            The star of the system.
        :type star:
            Star
        :param planet:
            The planet of the system.
        :type planet:
            Planet
        :param T_cor:
            Temperature of the corona in K, if already known.
        :type T_cor
        """

        v, veff, ne_planet, ne_star = cls._Parker(star=star, distance=planet.stardist, eccentricity=planet.eccentricity, T = T_cor)

        return cls(ne_planet, ne_star, v, veff, T_cor, Bsw={"planet": planet, "star": star, "vsw" : v, "eccentricity" : planet.eccentricity})

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
        ):

        mp = 1.660540210e-27  # kg
        mu0 = 4*np.pi * 1e-7
        try : 
            v, veff, ne, T = StellarWind._Parker(star=star, distance=d_alfven_point, eccentricity= eccentricity)

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
            
            res1 = (Br**2 + B_phi**2) / (2*mu0)
            res2 = ne * (v**2) * mp

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
        if d <= d_alfven_point :
            Br = magfield_surf * pow(R_star/(d * dua), 3)
        else :
            Br_ini = magfield_surf * pow(R_star/(d_alfven_point * dua), 3)
            Br = Br_ini * pow(d_alfven_point/d, 2)
        # dist = d * dua / Rstar
        # Br = magfield_surf * (1 + ((f-1)/pow(dist, 1.5))) / (f * pow(dist,2))
        
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
    def _mass_lossrate(t: float, R: float) -> float:
        """Compute the stellar mass loss rate normalized to the sun
        :param t:
            stellar age [yr]
        :type t:
            float
        :param R:
            radius of the star [R_sun]
        :type R:
            float
        """
        mp = 1.660540210e-27  # kg
        do = 1.49597870e11  # m

        vsun = StellarWind._calc_vsun(t)
        nsun = StellarWind._calc_nsun(t)
        M = 4 * np.pi * (do**2) * nsun * vsun * mp * (R**2)
        return M

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
            Tcor = 1.77 * np.power(M,-0.42) * np.power(P/Psun,0.52)
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
        e = eccentricity 
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
        Mls = StellarWind._mass_lossrate(t, star.radius)
        ne_planet = Mls / (4 * np.pi * ((d * dua) ** 2) * v * mp)
        ne_star = Mls / (4* np.pi * v * mp)
        if (ne_planet <= 0) or (ne_star <= 0):
            logger.error("Negative stellar wind density is not physical.")
            raise ValueError("Negative or null stellar wind density is not physical.")
        return (v, veff, ne_planet, ne_star)

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

