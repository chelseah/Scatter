#!/usr/bin/env python
from setting import *
import os
import numpy as np
import scipy as sp
from scipy.stats import rayleigh
def set_hill(mass_pl):
    i_hill=0
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



def writetobig(mass_pl,a_pl,d_pl,e_pl,i_pl,g_pl,n_pl,M_pl,outfile):
    fout=open(outfile,mode='w')
    #file header
    line=")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n"
    fout.write(line)
    line=") Lines beginning with `) are ignored.\n"
    fout.write(line)
    line=")---------------------------------------------------------------------\n"
    fout.write(line)
    line=" style (Cartesian, Asteroidal, Cometary) = Asteroidal\n"
    fout.write(line)
    line=" epoch (in days) = 0.\n"
    fout.write(line)
    line=")---------------------------------------------------------------------\n"
    fout.write(line)
    for i in xrange(N_pl):
        line='pl%d' % (i+1) + ' m = %15.5e' % (mass_pl[i]) + ' r=1.d0 d= %15.5e\n' % (d_pl[i])
        fout.write(line)
        line='%15.12e\t%15.12e\t%15.12e\n%15.12e\t%15.12e\t%15.12e\n' % (a_pl[i],e_pl[i],i_pl[i],g_pl[i],n_pl[i],M_pl[i])
        fout.write(line)
        line=' 0. 0. 0.\n'
        fout.write(line)
    fout.close()
    return

def submit(abspath,subfile):
    fout=open(subfile,mode='w')
    fout.write("#!/bin/csh\n")
    fout.write('#PBS -l nodes=1:ppn=1\n')
    fout.write('#PBS -q workq\n')
    fout.write('#PBS -r n \n')
    fout.write('#PBS -l walltime=48:00:00\n')
    fout.write('#PBS -N scatt_kepler\n')
    fout.write('cd %s\n' % abspath)
    fout.write('./mercury6 >& /dev/null \n')
    fout.write('./element6 \n')
    fout.close()
    os.system('qsub %s' % subfile)
    return

def main():
    for i in xrange(1,max_runs+1):
        rundir="run%d" % i
        os.system("mkdir %s" % rundir)
        for ex in excutable:
            os.system("cp %s %s/" % (ex,rundir))
        for infile in inlist:
            os.system("cp %s %s/" % (infile,rundir))

        #initial the basic param array for a system
        mass_pl=np.zeros(N_pl)
        a_pl=np.zeros(N_pl)
        d_pl=np.zeros(N_pl)

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
        d_pl[:2]=1.3
        d_pl[3:]=4.

        #set up other orbital elements
        e_pl=rayleigh.rvs(scale=sigma_e,size=6)
        i_pl=rayleigh.rvs(scale=sigma_i,size=6)
        g_pl=360.*np.random.randn(N_pl)
        n_pl=360.*np.random.randn(N_pl)
        M_pl=360.*np.random.randn(N_pl)
        t_max=t_orb*365.25*(a_inner)**1.5
        t_dump=t_max/1000.



        #write to big.in
        outfile = rundir+"/big.in"
        writetobig(mass_pl,a_pl,d_pl,e_pl,i_pl,g_pl,n_pl,M_pl,outfile)

        abspath=basename+rundir+"/"
        subfile=basename+rundir+"/submit_mercury"
        #subfile=rundir+"/submit_mercury"
        submit(abspath,subfile)
    return
if __name__=='__main__':
    main()
