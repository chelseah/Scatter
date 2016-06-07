import numpy as np
from numpy import random
np.random.seed()
#orbit settings
sigma_e=0.01
sigma_i=0.01
k_Hill=12.
P_inner=5.
doubling_time=1.e4
acc_time=doubling_time/1000.

#parallel settings
max_runs=1
num_proc = 1
nnodes=8

#default configuration
N_pl = np.random.random_integers(6,10)

#integration settings
integrator="hybarid"
t_max=1.0e6*(P_inner/365) 
Noutputs=1000


#path settings
basename='/home/cxhuang/mnt/scatt_Kepler/'
pythonpath='/home/cxhuang/src/virtualenv-1.5.2/ve/'
rundir='./'
subfile="qsubrebound"


#restart settings
restartinput="testeject.txt"
#restart from binary file
frombin=False
binfile="checkpoint.txt"

#output settings
checkpoint=True

#other flags
debug=False
addmass=False
gridsearch=False

#
statuscode={"eject":2,"star":3,"collision":1,"survive":0}
