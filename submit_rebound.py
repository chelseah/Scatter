#!/usr/bin/env python
import time as timing
from setting import *
import os
import numpy as np
import scipy as sp
import rebound
import sys
from scipy.stats import rayleigh
import multiprocessing as mp
import pickle

def check_for_bad_dt(sim):
    bad_dt = 0

    par = sim.particles
    p0 = sim.particles[0]                                      #star
    p = sim.particles[1]                            #planet
    dx = p.x - p0.x
    dy = p.y - p0.y
    dz = p.z - p0.z
    mh = (p.m/(3*p0.m))**(1./3.)
    rhill_p = mh*(dx*dx + dy*dy + dz*dz)**(0.5)               #hill radius planet
    vmax = 0                                                  #max relative velocity squared
    for i in xrange(2,sim.N):
        dx_i = par[i].x - p0.x
        dy_i = par[i].y - p0.y
        dz_i = par[i].z - p0.z
        mh_i = (par[i].m/(3*p0.m))**(1./3.)
        rhill_p_i = mh_i*(dx_i*dx_i + dy_i*dy_i + dz_i*dz_i)**(0.5)
        rh_sum = rhill_p + rhill_p_i
        dvx = par[i].vx - p.vx
        dvy = par[i].vy - p.vy
        dvz = par[i].vz - p.vz
        v = (dvx*dvx + dvy*dvy + dvz*dvz)**(0.5)
        print 'v',v,dvx,dvy,dvz
        if v > vmax:
            vmax = v
    HSR = sim.ri_hybarid.switch_radius
    min_dt = HSR*rh_sum / vmax
    print 'vmax',vmax,HSR*rh_sum,min_dt,HSR
    if(min_dt < 4*sim.dt):
        bad_dt = 1                                           #factor of 4 for extra wiggle room
    return bad_dt



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

def callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=0):
    #initialize a rebound run
    sim = rebound.Simulation()
    sim.t=t
    sim.add(m=1., r=0.005)
    for i in range(len(mass_pl)):
        sim.add(m = mass_pl[i], r = r_pl[i], a=a_pl[i], e=e_pl[i], inc=i_pl[i], Omega=Omega_pl[i], omega=omega_pl[i], M = M_pl[i],id=(i+1))
    sim.move_to_com()
    return sim

def submit(abspath,start=1,subfileauto=""):
    #create the submission file on the cluster
    if subfileauto=="":
        subfileauto=subfile
    fout=open(subfileauto,mode='w')
    fout.write("#!/bin/bash\n")
    fout.write('#PBS -l nodes=1:ppn=8\n')
    fout.write('#PBS -q workq\n')
    fout.write('#PBS -r n \n')
    fout.write('#PBS -l walltime=48:00:00\n')
    fout.write('#PBS -N rebound_kepler\n')
    fout.write('source %sbin/activate\n'% pythonpath)
    fout.write('cd %s\n' % abspath)
    if start>1:
        fout.write('%sbin/python submit_rebound.py %d' % (pythonpath,start))
    else:
        fout.write('%sbin/python submit_rebound.py' % pythonpath)
    fout.close()
    #os.system('qsub %s' % subfile)
    return

def orbit2str(particle):
    #write the orbit elements to a string
    orbit=particle.orbit
    string="%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"%(particle.m,particle.r,orbit.a,orbit.e,orbit.inc,orbit.Omega,orbit.omega,orbit.M)
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

def read_init(infile):
    #need to reload orbit elements from end result of a file.
    data=np.loadtxt(infile, skiprows=1)
    t=data[:,0][0]
    mass_pl=data[:,2]
    r_pl=data[:,3]
    a_pl=data[:,4]
    e_pl=data[:,5]
    i_pl=np.rad(data[:,6])
    Omega_pl=np.rad(data[:,7])
    omega_pl=np.rad(data[:,8])
    M_pl=np.rad(data[:,9])
    return [t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]

def init_orbit(randomstat=1):
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
    #a_pl[0]=1.91
    #a_pl[1]=2.31
    #a_pl[2]=1.07
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
    r_pl[:3]=5e-4
    r_pl[3:]=1e-4

    #set up other orbital elements
    e_pl=rayleigh.rvs(scale=sigma_e,size=6,random_state=randomstat)
    i_pl=rayleigh.rvs(scale=sigma_i,size=6,random_state=randomstat+1000)
    np.random.seed(randomstat+2000)
    omega_pl=2.*np.pi*np.random.rand(N_pl)
    Omega_pl=2.*np.pi*np.random.rand(N_pl)
    M_pl=2.*np.pi*np.random.rand(N_pl)

    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]

def integrate(sim,times,outfile):
    #the main integration routine
    Ncurrent = sim.N
    finalstatus=np.zeros(Ncurrent-1)
    nstep=np.zeros(Ncurrent-1)
    end=np.zeros([Ncurrent-1,8])
    bad_dts=np.zeros(len(times))
    dEs=np.zeros(len(times))
    for j,time in enumerate(times):
        sim.integrate(time)
        #for p in sim.particles:
            #print p
            #deal with Escape

        #print error
        max_d2 = 0.
        peject=None


        for p in sim.particles:
            if p.id==0:
                continue
            if p.a>100:
                mid = p.id
                peject=p
        if not peject is None:
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
            if p.id==0:
                continue
            end[p.id-1,:]=np.array(list(orbit2str(p).split()),dtype='f8')
        if not outfile is None:

            saveorbit(outfile,sim)#end)#sim)
        if checkpoint:
            checkpointfile=os.path.splitext(outfile)[0]+'.bin'
            sim.save(checkpointfile)
    

        bad_dts[j] = check_for_bad_dt(sim)
    	dEs[j] = sim.calculate_energy()
    return [finalstatus,end,nstep,bad_dts,dEs]

def one_run(runnumber,infile="",HSR=None,dt=None):
    #flow for one run

    #set up the output files and directories (need further modify)
    np.random.seed()

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
        if binfile=="":
            binfile=rundir+"rebound%.4d.bin" % runnumber

        sim=rebound.Simulation.from_file(binfile)
    #return
    saveorbit(outfile,sim)#save the initial orbits to output file file
    init=[]
    for p in sim.particles:
        if p.id==0:
            continue
        parr=np.array(list(orbit2str(p).split()),dtype='f8')
        init.append(parr)
    print init
    #fig = rebound.OrbitPlot(sim)


    # set up integrator (TO BE EDITED)
    #t_max=t_orb*365.25*(a_inner)**1.5

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
    #sim.exit_max_distance = 100.
    #sim.exit_min_distance = 0.01
    #print sim.collisions[0]



    times = np.logspace(np.log10(t+1),np.log10(t+t_max),Noutputs)
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
    pool = mp.Pool(processes=num_proc)
    pool.map(one_run,range(start,max_runs+start))

    return
if __name__=='__main__':
    if len(sys.argv)==1:
        #for run in xrange(1, num_runs+1):
        #    p = mp.Process(target=main)
        #    p.start()

        main()
    elif sys.argv[1]=='submit':
        abspath=basename
        submit(abspath)
    elif sys.argv[1]=='restart':
        one_run(1,restartinput)
    else:
        start=eval(sys.argv[1])
        main(start)
