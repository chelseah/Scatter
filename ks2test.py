#!/usr/bin/env python
import numpy as np
from scipy import stats
import pickle
import matplotlib.pylab as plt
from matplotlib import rc

def ks2(merc, merc_Exp2, rb, param, status='init', planets='all', torads=False, todegs=False):
    rc('text', usetex=True)
    rc('font', family='serif')
    rc('font', serif='cm')
    mercury = pickle.load(open(merc, 'rb'))
    mercury_Exp2 = pickle.load(open(merc_Exp2, 'rb'))
    rebound = pickle.load(open(rb, 'rb'))
    merc_paramdic = {"a":1,"e":2,"i":3,"peri":4,"node":5,"M":6,"mass":7}
    rb_paramdic = {"mass":0,"radius":1,"a":2,"e":3,"i":4,"peri":5,"node":6,"M":7}
    indexm = merc_paramdic[param]
    indexr = rb_paramdic[param]
    if status=='init':
        m = mercury[0]
        me = mercury_Exp2[0]
        r = rebound[0]
    elif status=='end':
        m = mercury[1]
        me = mercury_Exp2[1]
        npstatus = rebound[3]
        survive=npstatus==0
        r = rebound[1][survive,:]
    if planets=="all":
        m_tocomp = np.hstack(m[:,:,indexm])
        me_tocomp = np.hstack(me[:,:,indexm])
        if status=="end":
            r_tocomp = np.hstack(r[:,indexr])
        else:
            r_tocomp = np.hstack(r[:,:,indexr])
    elif planets=="inner":
        m_tocomp = np.hstack(m[:,3:,indexm])
        me_tocomp = np.hstack(me[:,3:,indexm])
        r_tocomp = np.hstack(r[:,3:,indexr])
    elif planets=="outer":
        m_tocomp = np.hstack(m[:,:3,indexm])
        me_tocomp = np.hstack(me[:,:3,indexm])
        r_tocomp = np.hstack(r[:,:3,indexr])
    else:
        m_tocomp = np.hstack(m[:,planets,indexm])
        me_tocomp = np.hstack(me[:,planets,indexm])
        r_tocomp = np.hstack(r[:,planets,indexr])
    if torads:
        m_tocomp = np.radians(m_tocomp)
        me_tocomp = np.radians(me_tocomp)
    elif todegs:
        r_tocomp = np.degrees(r_tocomp)
    dmme, pmme = stats.ks_2samp(m_tocomp, me_tocomp)
    dmr, pmr = stats.ks_2samp(m_tocomp, r_tocomp)
    dmmer, pmmer = stats.ks_2samp(me_tocomp, r_tocomp)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid('on')
    textstr='$\mathrm{p}_{\mathrm{m}/\mathrm{me}} = %.4g$\n$\mathrm{p}_{\mathrm{m}/\mathrm{r}} = %.4g$\n$\mathrm{p}_{\mathrm{me}/\mathrm{r}} = %.4g$'%(pmme, pmr, pmmer)
    nm, binsm, patchesm = ax.hist(m_tocomp, bins=20, histtype='step', normed=True, alpha=1.0, cumulative=False, label='Mercury')
    nme, binsme, patchesme = ax.hist(me_tocomp, bins=20, histtype='step', normed=True, alpha=1.0, cumulative=False, label='MercuryExp2')
    nr, binsr, patchesr = ax.hist(r_tocomp, bins=20, histtype='step', normed=True, alpha=1.0, cumulative=False, label='Rebound')
    nfake, binsfake, patchesfake = ax.hist([], histtype='step', alpha=0, edgecolor=None, label=textstr)
    if planets=="all":
        plt.xlabel('%s (%s)'%(param, status), fontsize=18)
    else:
        plt.xlabel('%s (%s) for %s planets'%(param, status, planets), fontsize=18)
    plt.ylabel('Normalized Count', fontsize=18)
    legend = ax.legend(numpoints=1, loc='best', fancybox=True)
    #plt.ylim(min(min(nm), min(nme), min(nr)), max(max(nm), max(nme), max(nr))*1.05)
    #plt.xlim(min(min(binsm), min(binsme), min(binsr)), max(max(binsm), max(binsme), max(binsr))*1.05)
    plt.show()
    return [dmme, pmme], [dmr, pmr], [dmmer, pmmer]
