#!/usr/bin/env python
import numpy as np
import scipy as sp
import os
import re
import commands
import pickle
def gather_run(runfile,nps,nele=8):
    init=np.zeros([int(nps),8])
    end=np.zeros([int(nps),8])
    nstep=np.zeros(nps)
    data=np.loadtxt(runfile)
    maxstep=0
    pend=0
    for i in xrange(nps):
        orbit_i=data[data[:,2]==(i+1)]
        init[i,:]=orbit_i[0,3:]
        end[i,:]=orbit_i[-1,3:]
        nstep[i]=orbit_i.shape[0]
        if nstep[i]>maxstep:
            maxstep=nstep[i]
            pend=i+1
    bad_dt=np.zeros(maxstep)
    #dEs=np.zeros(maxstep)
    orbit_i=data[data[:,1]==pend]
    dE=orbit_i[:,1]
    HSR=np.nan
    dt=np.nan
    time=48.*3600.

    survived=np.where(nstep==maxstep)[0]
    npcount=len(survived)
    necount=len(survived[survived>2])
    finalstatus=np.zeros(nps)
    collision=False
    for i in xrange(nps):
        if i in survived:
            continue
        else:
          if i==0:
            finalstatus[i]=2
          else:
            finalstatus[i]=1

    #get summary info
    #statuscode:
    #0-survive
    #1-collision
    #2-eject
    #3-star



    return [init,end,nstep,finalstatus,npcount,necount,HSR,dt,dE,bad_dt,time]

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
    rmax=160
    nps=6
    inpath="untracked_testing/run_expEmily/run_exp1_emily/"
    for i in xrange(rmin,rmax):
        runfile=inpath+"rebound%.4d.txt" % i
        pklfile=inpath+"runrebound%.4d.pkl" % i
        if not os.path.exists(pklfile):
    #for line in open("unfinished.txt"):
    #    pklfile=inpath+line.rstrip()
    #    runfile=inpath+os.path.splitext(line.rstrip())[0][3:]+".txt"
    #    if(1):
            print pklfile
            init,end,nstep,finalstatus,npcount,necount,HSR,dt,dE,bad_dt,time=gather_run(runfile,nps)
    	    datadump=[init,end,nstep,finalstatus,npcount,necount,HSR,dt,dE,bad_dt,time]
            #print datadump
            pickle.dump(datadump,open(pklfile,mode="w"))
            #break

    return

if __name__=='__main__':
    main()
