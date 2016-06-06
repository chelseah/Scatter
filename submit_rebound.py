#!/usr/bin/env python
import time as timing
import setting as st
from setting import *
import os
import numpy as np
import scipy as sp
import rebound
import sys
import pickle
from util import check_for_bad_dt,set_hill
from init import init_orbit, read_init
from submit_nodes import submit
from WJgrow import add_mass
def callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=0):
    #initialize a rebound run
    sim = rebound.Simulation()
    sim.t=t
    sim.add(m=1., r=0.005)
    for i in range(len(mass_pl)):
        sim.add(m = mass_pl[i], r = r_pl[i], a=a_pl[i], e=e_pl[i], inc=i_pl[i], Omega=Omega_pl[i], omega=omega_pl[i], M = M_pl[i],id=(i+1))
    sim.move_to_com()
    return sim

def orbit2str(particle):
    #write the orbit elements to a string
    orbit=particle.orbit
    string="%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"%(particle.m,particle.r,orbit.a,orbit.e,orbit.inc,orbit.Omega,orbit.omega,orbit.M)
    #print string
    return string

def saveorbit(outfile,sim):
    #save the orbits elements to a file
    fout=open(outfile,mode="a")
    E=sim.calculate_energy()
    for p in sim.particles:
        if p.id==0:
            continue
        line=orbit2str(p)
        fout.write("%f %15.12f %d %s\n" % (sim.t,E, p.id,line))
    fout.close()
    return


def integrate(sim,times,outfile):
    #the main integration routine
    finalstatus=np.zeros(N_pl)
    nstep=np.zeros(N_pl)
    end=np.zeros([N_pl,8])
    bad_dts=np.zeros(len(times))
    Ncurrent = sim.N
    dEs=np.zeros(len(times))
    for j,time in enumerate(times):
        try:
            if addmass: 
                add_mass(sim,time)

            sim.integrate(time)
            #deal with Escape
        except rebound.Escape as error:
            #print error
            max_d2 = 0.
            peject=None
            #check distance to be >1000, or (e>1 and distance>100)
            for p in sim.particles:
                if p.id==0:
                    continue
                d2 = p.x*p.x + p.y*p.y + p.z*p.z
                if d2>max_d2:
                    max_d2 = d2
                    mid = p.id
                    peject=p
       
            if not peject is None:
                if max_d2>1000:
                    print mid
                    end[mid-1,:]=np.array(list(orbit2str(peject).split()),dtype='f8')
                    sim.remove(id=mid)
                    nstep[mid-1]=int(sim.t/sim.dt)
                    Ncurrent-=1
                    finalstatus[mid-1]=statuscode['eject']
                    #print "final status",mid,"eject"
                elif max_d2>100:
                    orbit=particle.orbit
                    if orbit.e>1:
                        print mid
                        end[mid-1,:]=np.array(list(orbit2str(peject).split()),dtype='f8')
                        sim.remove(id=mid)
                        nstep[mid-1]=int(sim.t/sim.dt)
                        Ncurrent-=1
                        finalstatus[mid-1]=statuscode['eject']

        #deal with collision
        if Ncurrent>sim.N:
            #print "collision"
            for l in xrange(len(finalstatus)):

                if finalstatus[l]==0:
                    cflag=True
                    for p in sim.particles:
                        #print p,i+1,cflag
                        if p.id==0:
                            continue
                        if (p.id)==(l+1):
                            cflag=False
                            break
                    #print i,cflag
                    if cflag:
                        finalstatus[l]=statuscode['collision']
                        nstep[l]=int(sim.t/sim.dt)
                        #print "final status",i+1,'collision'
            Ncurrent=sim.N
        #print orbit2str(sim.particles[1].orbit)
        for p in sim.particles:
            print p.id
            if p.id==0:
                continue

            end[p.id-1,:]=np.array(list(orbit2str(p).split()),dtype='f8')
        if not outfile is None:

            saveorbit(outfile,sim)#end)#sim)
        if checkpoint:
            checkpointfile=os.path.splitext(outfile)[0]+'.bin'
            sim.save(checkpointfile)
    

        #bad_dts[j] = check_for_bad_dt(sim)
    	dEs[j] = sim.calculate_energy()
    return [finalstatus,end,nstep,bad_dts,dEs]

def one_run(runnumber,infile="",HSR=None,dt=None):
    #this function controls the flow for one run

    #set up the output files and directories (need further modify)

    #initialize the run
    t=0
    outfile=rundir+"rebound%.4d.txt" % runnumber
    fout=open(outfile,"w")
    fout.close()
    infofile=rundir+"runrebound%.4d.pkl" % runnumber
    if not frombin:
        if infile=="":
            mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=init_orbit(runnumber)
        else:
            t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=read_init(infile)

        sim = callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=t)
    else:
        if st.binfile=="":
            binfile=rundir+"rebound%.4d.bin" % runnumber

        sim=rebound.Simulation.from_file(binfile)
        t=sim.t
    #return
    saveorbit(outfile,sim)#save the initial orbits to output file file
    init=[]
    for p in sim.particles:
        if p.id==0:
            continue
        parr=np.array(list(orbit2str(p).split()),dtype='f8')
        init.append(parr)
    if debug:
        print init
        print HSR
        return


    # set up integrator (TO BE EDITED)

    if integrator=="hybarid":
    	sim.integrator="hybarid"
    	if not HSR is None:
    	    sim.ri_hybarid.switch_radius = HSR  #units of Hill radii
    	else:
    	    sim.ri_hybarid.switch_radius = 8  #units of Hill radii
    	sim.ri_hybarid.CE_radius = 20.  #X*radius

    	#set up time step
    	if not dt is None:
    	    sim.dt = dt #time step in units of yr/2pi
    	else:
    	    sim.dt = 0.001 #time step in units of yr/2pi
    sim.t=t
    #set up collision options
    #by default, the end result of the collision always
    #keep the small id number.
    sim.testparticle_type = 1
    sim.collision="direct"
    sim.collision_resolve = "merge"
    sim.collisions_track_dE = 1

    #set up escape options
    sim.exit_max_distance = 100.
    #sim.exit_min_distance = 0.01
    #print sim.collisions[0]

    ####create linear time step arr#
    N_acc = t_max/acc_time
    times=np.linspace(0,t_max*2*np.pi,N_acc)
    #################################
    #times = np.logspace(np.log10(t+1000),np.log10(t+t_max),Noutputs)
    E0 = sim.calculate_energy()
    start_t = timing.time()
    #call integration
    finalstatus,end,nstep,bad_dt,dEs=integrate(sim,times,outfile)

    if checkpoint:
        checkpointfile=rundir+"rebound%.4d.bin" % runnumber
        sim.save(checkpointfile)


    #Final processing
    #bad_dt = sim.ri_hybarid.timestep_too_large_warning
    dE = np.abs((dEs - E0)/E0)
    time = timing.time() - start_t
    #total number of planets left
    npcount=len(sim.particles)-1
    #total number of earth type planet left
    necount=0
    for p in sim.particles:
        if p.id==0:
            continue
        if p.m<0.5e-3:
            necount+=1
        nstep[p.id-1]=int(sim.t/sim.dt)


    datadump=[init,end,nstep,finalstatus,npcount,necount,HSR,dt,dE,bad_dt,time]

    def write_outcome(infofile,datadump):
        pickle.dump(datadump,open(infofile,"w"))
        return
    write_outcome(infofile, datadump)
    print 'exit run', runnumber
    return

def main(start=1):
    if gridsearch:
        #minpowdt,maxpowdt,numdt = [-4,-1,1] #min/max limits in logspace, i.e. 10**min - 10**max.
        minpowHSR,maxpowHSR,numHSR = [2,5,4]      #min/max limits in logspace, i.e. 10**min - 10**
        dt = P_inner/10./365.*2.*np.pi 
        HSRarr = np.linspace(minpowHSR,maxpowHSR,numHSR)
        repeat=int(nnodes/numHSR)
        #print nnodes,repeat
        one_run(start,HSR=HSRarr[int((start)/repeat)],dt=dt)
    else:
        one_run(start)
    return
if __name__=='__main__':
    if len(sys.argv)==1:
        
        main()
    elif sys.argv[1]=='submit':
        abspath=basename
        submit(abspath)
    elif sys.argv[1]=='restart':
        one_run(1,restartinput)
    else:
        start=eval(sys.argv[1])
        main(start)
