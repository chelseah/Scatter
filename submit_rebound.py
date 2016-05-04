#!/usr/bin/env python
from setting import *
import os
import numpy as np
import scipy as sp
import rebound
import sys
from scipy.stats import rayleigh
def set_hill(mass_pl):
    i_hill=0
    a_pl=np.zeros(3)
    while (a_pl[0]<1):
        while (i_hill<1):
            spacing=np.random.randn(3)
            a_inner=a_min+(a_max-a_min)*np.random.rand()
            a_max_i=a_max
            a_aux=a_inner+(a_max_i-a_inner)*spacing
            a_pl=np.sort(a_aux)

            R_hill_12=((mass_pl[0]+mass_pl[1])/3)**(1./3.)*(a_pl[0]+a_pl[1])/2.
            R_hill_23=((mass_pl[1]+mass_pl[2])/3)**(1./3.)*(a_pl[1]+a_pl[2])/2.
            diff_a_12=abs(a_pl[1]-a_pl[0])
            diff_a_23=abs(a_pl[2]-a_pl[1])

            if(diff_a_12>k_Hill*R_hill_12 and diff_a_23>k_Hill*R_hill_23):
                break
    return [a_inner,a_pl]

def callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl):
    sim = rebound.Simulation()
    sim.add(m=1., r=0.005)
    for i in range(len(mass_pl)):
        sim.add(m = mass_pl[i], r = r_pl[i], a=a_pl[i], e=e_pl[i], inc=i_pl[i], Omega=Omega_pl[i], omega=omega_pl[i], M = M_pl[i])
    return sim

def submit(abspath,subfile):
    fout=open(subfile,mode='w')
    fout.write("#!/bin/csh\n")
    fout.write('#PBS -l nodes=1:ppn=1\n')
    fout.write('#PBS -q workq\n')
    fout.write('#PBS -r n \n')
    fout.write('#PBS -l walltime=48:00:00\n')
    fout.write('#PBS -N rebound_kepler\n')
    fout.write('cd %s\n' % abspath)
    fout.write('./submit_rebound')
    #fout.write('./mercury6 >& /dev/null \n')
    #fout.write('./element6 \n')
    fout.close()
    os.system('qsub %s' % subfile)
    return

def orbit2str(particle):
    orbit=particle.orbit
    string="%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"%(particle.m,particle.r,orbit.a,orbit.e,orbit.inc,orbit.Omega,orbit.omega,orbit.f)
    return string

def saveorbit(outfile,sim):
    fout=open(outfile,mode="a")
    for p in range(len(sim.particles)):
        if p==0:
            continue
        line=orbit2str(sim.particles[p])
        fout.write("%f %d %s\n" % (sim.t,p,line))
    fout.close()
    return

def main():
    for i in xrange(1,max_runs+1):
        #rundir="run%d" % i
        #os.chdir(basename)
        #os.system("mkdir %s" % rundir)

        #initial the basic param array for a system
        mass_pl=np.zeros(N_pl)
        a_pl=np.zeros(N_pl)
        r_pl=np.zeros(N_pl) #changed to radius rather than density

        #initial semi-major axis and masses of gas giant,
        #in solar units
        mass_pl[0]=1.e-3
        mass_pl[1]=1.e-3
        mass_pl[2]=1.e-3

        a_inner,a_pl[:3]=set_hill(mass_pl[:3])
        #initial semi-major axis and masses of super earths,
        #in solar units
        M_earth=1./300./1000.
        mass_pl[3]=5*M_earth
        mass_pl[4]=10*M_earth
        mass_pl[5]=15*M_earth

        #need to move to setting file in the future
        a_pl[3]=0.1
        a_pl[4]=0.25
        a_pl[5]=0.5

        #need to modify in the future how this is set up
        r_pl[:2]=5e-4
        r_pl[3:]=1e-4

        #set up other orbital elements
        e_pl=rayleigh.rvs(scale=sigma_e,size=6)
        i_pl=rayleigh.rvs(scale=sigma_i,size=6)
        omega_pl=360.*np.random.randn(N_pl)
        Omega_pl=360.*np.random.randn(N_pl)
        M_pl=360.*np.random.randn(N_pl)
        #t_max=t_orb*365.25*(a_inner)**1.5
        t_max=1.e3
        Noutputs=1000.


        sim = callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl)

        init=[]
        for p in range(len(sim.particles)):
            if p==0:
                continue
            parr=np.array(list(orbit2str(sim.particles[p]).split()),dtype='f8')
            init.append(parr)
        print init
        fig = rebound.OrbitPlot(sim)



        # set up integrator (TO BE EDITED)
        sim.integrator="hybrid"
        sim.ri_hybarid.switch_ratio = 2  #units of Hill radii
        sim.ri_hybarid.CE_radius = 15.  #X*radius
        sim.testparticle_type = 1

        #set up time step
        sim.dt = 0.001 #time step in units of yr/2pi
        #set up collision options
        #by default, the end result of the collision always
        #keep the small id number.
        sim.collision="direct"
        sim.collision_resolve = "merge"
        sim.collisions_track_dE = 1

        #set up escape options
        sim.exit_max_distance = 100.
        #sim.exit_min_distance = 0.01
        #print sim.collisions[0]

        statuscode={"eject":2,"star":3,"collision":1,"survive":0}

        times = np.linspace(0,t_max,Noutputs)
        Ncurrent = sim.N
        finalstatus=np.zeros(Ncurrent-1)
        nstep=np.zeros(Ncurrent-1)
        end=np.zeros([Ncurrent-1,8])
        outfile="runrebound%.4d.txt" % i
        fout=open(outfile,"w") # looks like we open it but never write anything to it?
        fout.close()

        for j,time in enumerate(times):
            try:
                sim.integrate(time)
            #deal with Escape
            except rebound.Escape as error:

                #print error
                max_d2 = 0.
                peject=None
                for p in range(len(sim.particles)):
                    if p==0:
                        continue
                    d2 = sim.particles[p].x*sim.particles[p].x + sim.particles[p].y*sim.particles[p].y + sim.particles[p].z*sim.particles[p].z
                    if d2>max_d2:
                        max_d2 = d2
                        mid = p
                        peject=sim.particles[p]
                end[mid-1,:]=np.array(list(orbit2str(peject).split()),dtype='f8')
                sim.remove(id=mid)
                nstep[mid-1]=int(sim.t/sim.dt)
                Ncurrent-=1
                finalstatus[mid-1]=statuscode['eject']
                #print "final status",mid,"eject"
            #deal with collision
            if Ncurrent>sim.N:
                #print "collision"
                for l in xrange(len(finalstatus)):

                    if finalstatus[l]==0:
                        cflag=True
                        for p in range(len(sim.particles)):
                            #print p,i+1,cflag
                            if p==0:
                                continue
                            if (p)==(l+1):
                                cflag=False
                                break
                        #print i,cflag
                        if cflag:
                            finalstatus[l]=statuscode['collision']
                            nstep[l]=int(sim.t/sim.dt)
                            #print "final status",i+1,'collision'
                Ncurrent=sim.N
            #print orbit2str(sim.particles[1].orbit)
            for p in range(len(sim.particles)):
                if p==0:
                    continue
                end[p-1,:]=np.array(list(orbit2str(sim.particles[p]).split()),dtype='f8')
            if not outfile is None:

                saveorbit(outfile,sim)#end)#sim)
        #TBD:how to handel the collision with star

        import pickle
        #total number of planets left
        npcount=len(sim.particles)-1
        #total number of earth type planet left
        necount=0
        for p in range(len(sim.particles)):
            if p==0:
                continue
            if sim.particles[p].m<0.5e-3:
                necount+=1
            nstep[p-1]=int(sim.t/sim.dt)

        #print init
        #print end
        #print nstep
        #print finalstatus
        #print npcount
        #print necount
        datadump=[init,end,nstep,finalstatus,npcount,necount]
        print 'init=', init
        print 'end=', end
        print 'nstep=', nstep
        print 'finalstatus=', finalstatus
        print 'npcount=', npcount
        print 'necount=', necount
        #print datadump
        infofile="runrebound%.4d.pkl" % i
        def write_outcome(infofile,datadump):
            pickle.dump(datadump,open(infofile,"w"))
            return
        write_outcome(infofile, datadump)

        

    return
if __name__=='__main__':
    if len(sys.argv)==1:
        main()
    elif sys.argv[1]=='submit':
        abspath=basename
        subfile="qsubrebound"
        submit(abspath,subfile)