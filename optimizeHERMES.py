import rebound
import numpy as np
from submit_rebound import one_run


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
    one_run(params[0],HSR=prams[1],tmax=params[2])
    return



def main(start):
    n_trials = 8
    dt_arr = [-4,0,n_trials]      #min/max limits in logspace, i.e. 10**min - 10**max.
    HSR_arr = [0,2,n_trials]      #min/max limits in logspace, i.e. 10**min - 10**
    minpowdt,maxpowdt,numdt = dt_array
    minpowHSR,maxpowHSR,numHSR = HSR_array
    dt = np.logspace(minpowdt,maxpowdt,numdt)
    HSR = np.logspace(minpowHSR,maxpowHSR,numHSR)
    pool = mp.Pool(processes=8)
    params=[]
    for i in xrange(8):
        param=[i+1,dt[start],HSR[i]]
        params.append(param)
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
    else:
        start=eval(sys.argv[1])
        main(start)


