#!/usr/bin/env python
import numpy as np
import scipy as sp
import matplotlib
from matplotlib import pyplot as plt
import os
import re
import commands
import pickle
def gather_run(rundir,nps,nele=8):
    initarr=np.zeros([int(nps),8])
    endarr=np.zeros([int(nps),8])
    nsteparr=np.zeros(nps)
    for i in xrange(nps):
        infile=rundir+"pl%d.aei" % (i+1)
        #print infile
        data=np.loadtxt(infile,skiprows=4)
        #get init condition
        #print data[0,:]
        try:
            initarr[i,:]=data[0,:]
        except IndexError:
            initarr[i,:]=data
            endarr[i,:]=data
            nsteparr[i]=1
            continue

        #get end result
        endarr[i,:]=data[-1,:]
        nsteparr[i]=data.shape[0]

    #get summary info
    #statuscode:
    #0-survive
    #1-collision
    #2-eject
    #3-star

    npstatus=np.zeros(nps)
    if (nsteparr==1).any():
        #bad run
        npstatus+=np.nan
        return [initarr,endarr,nsteparr,npstatus]



    infile=rundir+"info.out"
    status,outputs=commands.getstatusoutput("grep ^' pl' %s" % infile)
    #print outputs
    lensum=len(outputs.split('\n'))
    #print lensum
    #print infile,lensum,outputs
    if outputs=="":
        return [initarr,endarr,nsteparr,npstatus]
    for i in xrange(lensum):
        line = outputs.split('\n')[i]
        n,sc=classify(line)
        if sc==1:
            nc=int(line.split()[4][2])
            if nsteparr[n-1]<nsteparr[nc-1]:
                npstatus[n-1]=sc
            else:
                npstatus[nc-1]=sc
        else:
            npstatus[n-1]=sc



    return [initarr,endarr,nsteparr,npstatus]

def classify(line):
    n=int(line.split()[0][2])
    if 'ejected' in line:
        sc=2
        return [n,sc]
    if 'central' in line:
        sc=3
        return [n,sc]
    if 'hit' in line:
        if 'ast' in line:
            return [n,0]
        return [n,1]
        #nc=int(line.split()[4][2])
        #if nc<n:
        #    return [n,1]
        #else:
        #    return [nc,1]
    return

def output_status(statustotal):

    print "#Summary of the run:"
    print "Nexpriments: %d" % statustotal.shape[0]
    SE,SJ,EE,JJ,EJ,EEE,EEJ,EJJ,JJJ,EEEJ,EEJJ,EJJJ,EEEJJ,EEJJJ,Nsix,Nbad=[0]*16

    npcount=np.zeros(statustotal.shape[0])
    necount=np.zeros(statustotal.shape[0])
    for i in xrange(statustotal.shape[0]):
        if np.isnan(statustotal[i,:]).any():
            Nbad+=1
            npcount[i]=-1
            necount[i]=-1
        index=np.where(statustotal[i,:]==0)[0]
        if len(index)==6:
            Nsix+=1
            npcount[i]=6
            necount[i]=3
        if len(index)==5:
            if(index[2]<3):
                EEJJJ+=1
                necount[i]=2
            else:
                EEEJJ+=1
                necount[i]=3
            npcount[i]=5
        if len(index)==4:
            if(index[2]<3):
                EJJJ+=1
                necount[i]=1
            elif(index[1]<3):
                EEJJ+=1
                necount[i]=2
            else:
                EEEJ+=1
                necount[i]=3
            npcount[i]=4
        if len(index)==3:
            if (index<3).all():
                JJJ+=1
            elif (index<3).any():
                if index[1]<3:
                    EJJ+=1
                    necount[i]=1
                else:
                    EEJ+=1
                    necount[i]=2
            else:
                EEE+=1
                necount[i]=3
            npcount[i]=3
        if len(index)==2:
            if (index<3).all():
                JJ+=1
            elif (index<3).any():
                EJ+=1
                necount[i]=1
            else:
                EE+=1
                necount[i]=2
            npcount[i]=2
        if len(index)==1:
            if index[0]<3:
                SJ+=1
            else:
                SE+=1
                necount[i]=1
            npcount[i]=1
    print "Nsingle system:%d (total) %d (E) %d (J)" % (SE+SJ,SE,SJ)
    print "Ntwo planet system:%d (total) %d (EE) %d (JJ) %d (EJ)" % (EE+JJ+EJ,EE,JJ,EJ)
    print "Ntriple planet system: %d (total) %d (EEE) %d (EEJ) %d (EJJ) %d (JJJ)" % (EEE+EEJ+EJJ+JJJ,EEE,EEJ,EJJ,JJJ)
    print "Nfour planet system: %d (total) %d (EEEJ) %d (EEJJ) %d (EJJJ) " % (EEEJ+EEJJ+EJJJ,EEEJ,EEJJ,EJJJ)
    print "Nfive planet system: %d (total) %d (EEEJJ) %d (EEJJJ) " %(EEEJJ+EEJJJ,EEEJJ,EEJJJ)
    print "Nsix planet system: %d (total)" % (Nsix)
    print "Bad run: %d" % (Nbad)

    #output specific runs
    print npcount
    print 'five planets',np.where(npcount==5)[0]+1
    print 'four planets',np.where(npcount==4)[0]+1
    print 'three planets',np.where(npcount==3)[0]+1
    print 'two planets',np.where(npcount==6)[0]+1
    print 'three earth',np.where(necount==3)[0]+1
    print 'two earth',np.where(necount==2)[0]+1
    print 'one earth',np.where(necount==1)[0]+1
    return [npcount,necount]

def main():
    rmin=1
    rmax=151
    nps=6
    inpath="/home/xuhuang/Project/data/Kepler/scatt_kepler/"
    inittotal=np.zeros([rmax-rmin,nps,8])
    endtotal=np.zeros([rmax-rmin,nps,8])
    nsteptotal=np.zeros([rmax-rmin,nps])
    npstatustotal=np.zeros([rmax-rmin,nps])
    for i in xrange(rmin,rmax):
        rundir=inpath+"run%d/" % i
        initarr,endarr,nsteparr,npstatus=gather_run(rundir,nps)
        inittotal[i-rmin,:]=initarr
        endtotal[i-rmin,:]=endarr
        nsteptotal[i-rmin,:]=nsteparr
        npstatustotal[i-rmin,:]=npstatus
        #print result

    #analyze the result
    plot=False
    if (plot):
        plot_init(inittotal)

        plot_init_survive(inittotal,nsteptotal)

        plot_init_status(inittotal,npstatustotal)

        plot_init_end(inittotal,endtotal)

        plot_end_survive(endtotal,nsteptotal)

    npcount,necount=output_status(npstatustotal)
    datadump=[inittotal,endtotal,nsteptotal,npstatustotal,npcount,necount]
    pickle.dump(datadump,open('runsummary.pkl','w'))
    return

if __name__=='__main__':
    main()
