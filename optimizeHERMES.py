import rebound
import numpy as np
from submit_rebound import one_run
import sys
import multiprocessing as mp
def submit(abspath,subfile,start=1):
    #create the submission file on the cluster
    fout=open(subfile,mode='w')
    fout.write("#!/bin/csh\n")
    fout.write('#PBS -l nodes=1:ppn=8\n')
    fout.write('#PBS -q workq\n')
    fout.write('#PBS -r n \n')
    fout.write('#PBS -l walltime=06:00:00\n')
    fout.write('#PBS -N rebound_kepler\n')
    fout.write('source %sbin/activate\n'% pythonpath)
    fout.write('cd %s\n' % abspath)
    fout.write('%sbin/python optimizeHERMES.py %d' % (pythonpath,start))
    fout.close()
    #os.system('qsub %s' % subfile)
    return



def sim_problem(params):
    one_run(params[0],HSR=params[2],dt=params[1])
    return



def main(start):
    n_trials = 8
    dt_arr = [-4,-1,n_trials]      #min/max limits in logspace, i.e. 10**min - 10**max.
    HSR_arr = [np.log10(2),2,n_trials]      #min/max limits in logspace, i.e. 10**min - 10**
    minpowdt,maxpowdt,numdt = dt_arr
    minpowHSR,maxpowHSR,numHSR = HSR_arr
    dt = np.logspace(minpowdt,maxpowdt,numdt)
    HSR = np.logspace(minpowHSR,maxpowHSR,numHSR)
    print dt
    print HSR
    pool = mp.Pool(processes=8)
    params=[]
    for i in xrange(8):
        param=[start+i,dt[int(start/8.)],HSR[i]]
        params.append(param)
        print param
        #sim_problem(param)
    TASK=[(params[i]) for i in xrange(len(params))] 
    pool.map(sim_problem,params)
    return

def submitnnode():
    nnodes=8
    for i in xrange(nnodes):
        print i
        submit(basename,"optimize_%d" % i,start=1+8*i)
        os.system('qsub %s' % subfile)

if __name__=='__main__':
    if  len(sys.argv)==1:
        submitnnode()
    elif sys.argv[1]=="test":
        sim_problem([1,1.e-4,0])
    else:
        start=eval(sys.argv[1])
        main(start)


