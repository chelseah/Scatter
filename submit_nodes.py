#!/usr/bin/env python
from setting import *
import os

def submit(abspath,start=1,subfileauto=""):
    #create the submission file on the cluster
    if subfileauto=="":
        subfileauto=subfile
    fout=open(subfileauto,mode='w')
    fout.write("#!/bin/bash\n")
    fout.write('#PBS -l nodes=1:ppn=%d\n' % max_runs)
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



def main():
    #nnodes=20
    for i in xrange(nnodes):
        print i
        subfile=basename+"qsubrebound_%d" % i
        submit(basename,subfileauto="qsubrebound_%d" % i,start=1+max_runs*i)
        os.system('qsub %s' % subfile)

if __name__=='__main__':
    main()
