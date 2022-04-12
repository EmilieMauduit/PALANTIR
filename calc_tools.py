import os
import scipy.io as sio
from scipy.io import readsav
import numpy as np
from math import *

def find_value(a,l, res=0):
    """Retourne la valeur ou l'indice de la valeur la plus proche de a dans la liste l, par défaut l'indice."""
    min=np.abs(l[0]-a)
    indmin=0
    for i in range(1,len(l)):
        nmin=np.abs(l[i]-a)
        if (nmin <= min):
            min=nmin
            indmin=i
    if res==0 :
        return(indmin)
    else :
        return (l[indmin])

def synchro_dist(w,Qp,tsync,Mp,Rp,Ms):
    """ Calcule la distance à laquelle doit se trouver la planète pour qu'il y ait synchronisationd de l'orbite selon l'âge du système"""
    G=6.6725985e-11 
    a=9*tsync/(4*w*0.26*Qp)
    b=G*Mp/pow(Rp,3)
    c=pow(Ms/Mp,2)
    res=Rp*pow(a*b*c,1/6)
    return(res)