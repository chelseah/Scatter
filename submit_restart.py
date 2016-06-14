#!/usr/bin/env python
import os
from setting import *
from submit_rebound import submit
def main():
    rmin=0
    rmax=max_runs*nnodes
    inpath=rundir
    for i in xrange(rmin,rmax):
        runfile=inpath+"rebound%.4d.txt" % i
        pklfile=inpath+"runrebound%.4d.pkl" % i
        if not os.path.exists(pklfile):
            subfile=basename+"qsubrebound_%d" % i
            submit(basename,subfileauto="qsubrebound_%d" % i, start=1+max_runs*i)
            os.system('qsub %s' % subfile)
        #break
    return


if __name__=='__main__':
    main()
