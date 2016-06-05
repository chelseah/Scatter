#!/usr/bin/env python
import numpy as np
import scipy as sp
from setting import *
from scipy.stats import rayleigh

def init_orbit(randomstat=1):
    #initial the basic param array for a system
    
    mass_pl=np.zeros(N_pl)
    a_pl=np.zeros(N_pl)
    r_pl=np.zeros(N_pl) #changed to radius rather than density

    #a_inner,a_pl[:3]=set_hill(mass_pl[:3])

    #initial semi-major axis and masses of super earths,
    #in solar units
    M_earth=1./300./1000.
    R_earth = 6371./1.49e8
    density = 3*M_earth/(4*np.pi*R_earth**3)
    loop = True

    while (loop==True):
    #set up other orbital elements
        e_pl=rayleigh.rvs(scale=0.01,size=N_pl,random_state=randomstat)
        i_pl=rayleigh.rvs(scale=2*(np.pi/180.),size=N_pl,random_state=randomstat+1000)
        np.random.seed(randomstat+2000)
        omega_pl=2.*np.pi*np.random.rand(N_pl)
        Omega_pl=2.*np.pi*np.random.rand(N_pl)
        M_pl=2.*np.pi*np.random.rand(N_pl)
        mass_pl = np.random.normal(loc=5*M_earth,scale = 2*M_earth,size=N_pl)
        
        r_pl =2*R_earth+r_pl#(3*mass_pl/(4*np.pi*density))**(1./3.)
        K = np.random.normal(loc=k_Hill,scale=2)
        P_min =rayleigh.rvs(scale=P_inner/(365),random_state=randomstat+3000)
        a_min = P_min**(2./3.)

        in_range=[]
        a_pl[0] = a_min
        loop = False
        for i in range(1,N_pl):
            m_ratio = 2*(3/(mass_pl[i-1]+mass_pl[i]))**(1./3.)
            a_pl[i] =(K/m_ratio + 1)/(1 - K/m_ratio)*a_pl[i-1]
            if a_pl[i]>0.1 and a_pl[i]<0.4:
                in_range.append(i)
        #if len(in_range)>0:
        #    loop=False
        #    accreting_planet = np.random.choice(in_range)
        #    mass_pl[accreting_planet]=12*M_earth
        #    for i in range(accreting_planet,N_pl):
        #        if i>0:
        #            m_ratio = 2*(3/(mass_pl[i-1]+mass_pl[i]))**(1./3.)
        #            a_pl[i] =(K/m_ratio + 1)/(1 - K/m_ratio)*a_pl[i-1]
    #print mass_pl, a_pl, r_pl, e_pl, i_pl, omega_pl, Omega_pl, M_pl
    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]

def read_init(infile):
    #need to reload orbit elements from end result of a file.
    data=np.loadtxt(infile, skiprows=1)
    t=data[:,0][0]
    mass_pl=data[:,3]
    r_pl=data[:,4]
    a_pl=data[:,5]
    e_pl=data[:,6]
    i_pl=data[:,7]
    Omega_pl=data[:,8]
    omega_pl=data[:,9]
    M_pl=data[:,10]
    return [t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]


