#!/usr/bin/env python
import numpy as np
import scipy as sp

def mass_growth(time, initial_mass,doubling_time):
    new_mass = initial_mass*(2.0)**(time/doubling_time)
    return new_mass

def add_mass(sim,time):
    if time<= 4.5*doubling_time*2*np.pi:

        doubling_time = 1.0e4
        acc_time = doubling_time/1000
        t_max = 5e4#doubling_time*4.5 + 1.e5
        N_acc = t_max/acc_time

        sp = sim.particles
        mass_list = []
        M_earth = 1./300./1000.
        threshold = 12*M_earth

        for i in range(1,sim.N):
            mass_list.append(sp[i].m)
        accreting_planet = np.where(np.array(mass_list)>=threshold)[0]+1
        for i in accreting_planet:
            planet_density = (3*sp[i].m/(4*np.pi*sp[i].r**3))
            sp[i].m = mass_growth(time,sp[i].m,doubling_time)
            sp[i].r = (3*sp[i].m/(4*np.pi*sp[i].r))**(1./3.)
        sim.move_to_com()
    return

