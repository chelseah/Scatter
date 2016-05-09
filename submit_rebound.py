#!/usr/bin/env python
from setting import *
import os
import numpy as np
import scipy as sp
import rebound
import sys
from scipy.stats import rayleigh
import multiprocessing as mp
import pickle
def set_hill(mass_pl):
    #set up the giant planets so that they will be instable at one point
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

def callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=0):
    #initialize a rebound run
    sim = rebound.Simulation()
    sim.t=t
    sim.add(m=1., r=0.005)
    for i in range(len(mass_pl)):
        sim.add(m = mass_pl[i], r = r_pl[i], a=a_pl[i], e=e_pl[i], inc=i_pl[i], Omega=Omega_pl[i], omega=omega_pl[i], M = M_pl[i],id=(i+1))
    #sim.move_to_com()
    for p in sim.particles:
        print p
    return sim

def submit(abspath,subfile):
    #create the submission file on the cluster
    fout=open(subfile,mode='w')
    fout.write("#!/bin/csh\n")
    fout.write('#PBS -l nodes=1:ppn=1\n')
    fout.write('#PBS -q workq\n')
    fout.write('#PBS -r n \n')
    fout.write('#PBS -l walltime=06:00:00\n')
    fout.write('#PBS -N rebound_kepler\n')
    fout.write('source /home/edeibert/src/virtualenv-1.5.2/ve/bin/activate\n')
    fout.write('cd %s\n' % abspath)
    fout.write('/home/edeibert/src/virtualenv-1.5.2/ve/bin/python submit_rebound.py')
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
    for p in sim.particles:
        if p.id==0:
            continue
        line=orbit2str(p)
        fout.write("%f %d %s\n" % (sim.t,p.id,line))
    fout.close()
    return

def read_init(infile):
    #need to reload orbit elements from end result of a file.
    data=np.loadtxt(infile)
    t=data[:,0][0]
    mass_pl=data[:,2]
    r_pl=data[:,3]
    a_pl=data[:,4]
    e_pl=data[:,5]
    i_pl=data[:,6]
    Omega_pl=data[:,7]
    omega_pl=data[:,8]
    M_pl=data[:,9]
    return [t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl] 

def init_orbit():
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
    omega_pl=2.*np.pi*np.random.randn(N_pl)
    Omega_pl=2.*np.pi*np.random.randn(N_pl)
    M_pl=2.*np.pi*np.random.randn(N_pl)

    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl] 

def integrate(sim,times,outfile):
    #the main integration routine
    Ncurrent = sim.N
    finalstatus=np.zeros(Ncurrent-1)
    nstep=np.zeros(Ncurrent-1)
    end=np.zeros([Ncurrent-1,8])
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
    return [finalstatus,end,nstep]

def one_run(runnumber,infile=""):
    #flow for one run

    #set up the output files and directories (need further modify)
    np.random.seed()
    
    #initialize the run
    t=0
    if infile=="":
        mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=init_orbit()
    else:
        t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=read_init(infile)
    rundir="run%.4d" % runnumber
    os.system("mkdir %s" % rundir)
    os.chdir(rundir)

    outfile="rebound%.4d.txt" % runnumber
    fout=open(outfile,"w")
    fout.close()
    infofile="runrebound%.4d.pkl" % runnumber


    
    sim = callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=t)
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
    t_max=1.e5
    Noutputs=10000.

    sim.integrator="hybrid"
    sim.ri_hybarid.switch_ratio = 10  #units of Hill radii
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
    #sim.exit_max_distance = 100.
    #sim.exit_min_distance = 0.01
    #print sim.collisions[0]


    times = np.linspace(0,t_max,Noutputs)
    
    #call integration
    finalstatus,end,nstep=integrate(sim,times,outfile)


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

    datadump=[init,end,nstep,finalstatus,npcount,necount]
    
    def write_outcome(infofile,datadump):
        pickle.dump(datadump,open(infofile,"w"))
        return
    write_outcome(infofile, datadump)
    os.chdir('../')
    print 'exit run', runnumber
    return

def main():
    pool = mp.Pool(processes=num_proc)
    pool.map(one_run,range(1,max_runs+1))

    return
if __name__=='__main__':
    if len(sys.argv)==1:
        #for run in xrange(1, num_runs+1):
        #    p = mp.Process(target=main)
        #    p.start()
        
        main()
    elif sys.argv[1]=='submit':
        abspath=basename
        subfile="qsubrebound"
        submit(abspath,subfile)
    elif sys.argv[1]=='restart':
        abspath=basename
        infile="testinput.txt"
        one_run(1,infile)
