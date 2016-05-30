#!/usr/bin/env python
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import os
import re
import commands
import pickle
import glob
#from mutual_inclination import *
from matplotlib import rc

def gather_run(infile):
    data = pickle.load(open(infile, 'rb'))
    init = data[0]
    end = data[1]
    endtime = data[2]
    npstatus = data[3]
    npcount = data[4]
    se = np.array([3,4,5])
    remain = np.where(npstatus==0)[0]
    necount = len(np.intersect1d(remain, se)) # temporary fix since necount isn't working for some
    #necount = data[5]
    return init, end, endtime, npstatus, npcount, necount

def gather_txtrun(infile):
    """To use when runs didn't finish and only .txt is available."""
    data = np.loadtxt(infile)
    init = data[:6, 3:]
    remaining = np.asarray(list(set(data[-6:, 2])))
    end = np.zeros([6, 8])
    #end = data[-(len(remaining)):, 3:]
    endt = np.zeros(6)
    for i in range(len(endt)):
        ind = np.where(data[:, 2]==(i+1))[0][-1]
        endt[i] = data[ind, 0]
        end[i, :] = data[ind, 3:]
    npstatus = np.zeros(6)
    for i in range(len(npstatus)):
        if (i+1) in remaining:
            npstatus[i] = 0
        else:
            npstatus[i] = 3
    npcount = len(remaining)
    SEs = [4, 5, 6]
    necount = len(list(set(remaining) & set(SEs)))
    return init, end, endt, npstatus, npcount, necount

def data_dump(infiles, nps, t='pkl'):
    runs = len(infiles)
    inittotal = np.zeros([runs, nps, 8])
    endtotal = np.zeros([runs, nps, 8])
    endttotal = np.zeros([runs, nps])
    npstatustotal = np.zeros([runs, nps])
    npcounttot = np.zeros(runs)
    necounttot = np.zeros(runs)
    if t=='pkl':
        for i in range(runs):
            init, end, endtime, npstatus, npcount, necount = gather_run(infiles[i])
            inittotal[i, :, :] = init
            endtotal[i, :, :] = end
            endttotal[i, :] = endtime
            npstatustotal[i, :] = npstatus
            npcounttot[i] = npcount
            necounttot[i] = necount
    if t=='txt':
        for i in range(runs):
            init, end, endtime, npstatus, npcount, necount = gather_txtrun(infiles[i])
            inittotal[i, :, :] = init
            endtotal[i, :, :] = end
            endttotal[i, :] = endtime
            npstatustotal[i, :] = npstatus
            npcounttot[i] = npcount
            necounttot[i] = necount
    return inittotal, endtotal, endttotal, npstatustotal, npcounttot, necounttot

def plot_scatter(dataarr,param1,param2,xlim=[],ylim=[],color=''):
    rc('text', usetex=True)
    rc('font', family='serif')
    rc('font', serif='cm')
    paramdic={"mass":0,"radius":1,"a":2,"e":3,"i":4,"peri":5,"node":6,"M":7}
    index1=paramdic[param1]
    index2=paramdic[param2]
    fig=plt.figure(figsize=[9,5])
    ax=fig.add_subplot(111)
    try:
        cax=ax.plot(dataarr[:,:,index1],dataarr[:,:,index2],'.')
    except IndexError:
        #we allow user to use one of the columns to color code the scatter.
        if not color=='':
            if type(color) is str:
                cax=ax.scatter(dataarr[:,index1],dataarr[:,index2],c=dataarr[:,paramdic[color]], cmap='hot')
            else:
                cax=ax.scatter(dataarr[:,index1],dataarr[:,index2],c=color, cmap='hot')
        else:
            cax=ax.plot(dataarr[:,index1],dataarr[:,index2],'.')
    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    if len(xlim)>0:
        ax.set_xlim(xlim)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    if not color=="":
        cbar = fig.colorbar(cax,shrink=0.5,aspect=10,pad=0.05)
    else:
        legend = fig.legend(cax,tuple(np.arange(1,7).astype(str)),bbox_to_anchor=(1,0.75),borderpad=1,numpoints=1)
    plt.show()
    return

def plot_hist(ncount, str):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    bins=np.arange(0,8)-0.5
    n,bins,patches=ax.hist(ncount,bins=bins,rwidth=0.5)
    plt.xlabel(str)
    print n, bins
    plt.show()
    return

def main():
    os.chdir('run_exp1')
    infiles = glob.glob('*.pkl')
    inittotal, endtotal, endttotal, npstatustotal, npcounttot, necounttot = data_dump(infiles, 6)
    os.chdir('../')
    plot_scatter(inittotal, "a", "e")
    plot_scatter(endtotal, "a", "e", xlim=[0, 10], ylim=[0, 1])
    survive=npstatustotal==0
    plot_scatter(endtotal[survive,:], "a", "e", xlim=[0,10], ylim=[0,1],color="mass")
    SEsys=necounttot>0
    survive=npstatustotal[SEsys,:]==0
    plot_scatter(endtotal[SEsys,:,:][survive,:], "a", "e", xlim=[0,100], ylim=[0, 1], color="mass")
    plot_hist(npcounttot, 'Number of Planets')
    plot_hist(necounttot, 'Number of Super Earths')
    necount_all = np.rollaxis(np.tile(necounttot, (8,1)),-1)
    plot_scatter(endtotal[SEsys,:,:][survive,:],"a", "i", xlim=[0, 10], ylim=[0, 0.8], color=necount_all[SEsys,:][survive])
    #datadump=[inittotal,endtotal,endttotal,npstatustotal,npcounttot,necounttot]
    #pickle.dump(datadump,open('runsummary_Rebound_exp1.pkl','w'))
    return

if __name__ == "__main__":
    main()
