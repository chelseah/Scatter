#!/usr/bin/env python
import numpy as np
import scipy as sp
from setting import *
from scipy.stats import rayleigh
from period_dist import semi_major


def set_hill(mass_pl):
    #set up the giant planets so that they will be instable at one point
    i_hill=0
    a_pl=np.zeros(3)
    while (a_pl[0]<1):
        while (i_hill<1):
            spacing=np.random.rand(3)
            print 'spacing,',spacing
            a_inner=a_min+(a_max-a_min)*np.random.rand()
            a_max_i=a_max
            print 'amax_i,',a_max_i
            a_aux=a_inner+(a_max_i-a_inner)*spacing
            print 'a_aux,',a_aux
            a_pl=np.sort(a_aux)

            R_hill_12=((mass_pl[0]+mass_pl[1])/3)**(1./3.)*(a_pl[0]+a_pl[1])/2.
            R_hill_23=((mass_pl[1]+mass_pl[2])/3)**(1./3.)*(a_pl[1]+a_pl[2])/2.
            diff_a_12=abs(a_pl[1]-a_pl[0])
            diff_a_23=abs(a_pl[2]-a_pl[1])

            if(diff_a_12>k_Hill*R_hill_12 and diff_a_23>k_Hill*R_hill_23):
                break
    return [a_inner,a_pl]

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



def init_orbit(randomstat=1):
    #initial the basic param array for a system
    mass_pl=np.zeros(N_pl)
    a_pl=np.zeros(N_pl)
    r_pl=np.zeros(N_pl) #changed to radius rather than density


    #initial semi-major axis and masses of gas giant,
    #in solar units
    gpmass=3
    sesemi=1
    semass=2
    #fixed
    if gpmass==1:
        mass_pl[0]=1.e-3
        mass_pl[1]=1.e-3
        mass_pl[2]=1.e-3
    elif gpmass==2:
        #uniform in a mass range
        mass_pl[:3]=(0.5+np.random.random(3)*1.5)*1.e-3
    elif gpmass==3:
        mass_pl[:3]=(0.3+np.random.random(3)*2.7)*1.e-3

    a_inner,a_pl[:3]=set_hill(mass_pl[:3])
     ##need to move to setting file in the future
    if sesemi==1:
        a_pl[3]=0.1
        a_pl[4]=0.25
        a_pl[5]=0.5
    elif sesemi==2:
        a_pl[3:]=semi_major()

    #initial semi-major axis and masses of super earths,
    #in solar units
    M_earth=1./300./1000.

    if semass==1:
        mass_pl[3]=5*M_earth
        mass_pl[4]=10*M_earth
        mass_pl[5]=15*M_earth
    elif semass==2:
        marr=np.array([5.,10.,15.])
        sortindex=np.random.choice(3,3,False)
        print sortindex
        for i in xrange(3):
	        mass_pl[3+i]=marr[sortindex[i]]*M_earth

    r_pl[:3]=5e-4
    r_pl[3:]=1e-4

    #set up other orbital elements
    e_pl=rayleigh.rvs(scale=sigma_e,size=N_pl,random_state=randomstat)
    i_pl=rayleigh.rvs(scale=sigma_i,size=N_pl,random_state=randomstat+1000)
    np.random.seed(randomstat+2000)
    omega_pl=2.*np.pi*np.random.rand(N_pl)
    Omega_pl=2.*np.pi*np.random.rand(N_pl)
    M_pl=2.*np.pi*np.random.rand(N_pl)

    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]



