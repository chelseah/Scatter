import numpy as np

#orbit settings
sigma_e=0.01
sigma_i=0.01
a_min=2.
a_max=5.
a_max_min=3.
k_Hill=3.

#parallel settings
max_runs=8
num_proc = 8
nnodes=20

#default configuration,never change
t_orb=1.e7
N_pl=6
statuscode={"eject":2,"star":3,"collision":1,"survive":0}

#path settings
basename='/home/cxhuang/mnt/scatt_Kepler/'
pythonpath='/home/cxhuang/src/virtualenv-1.5.2/ve/'
rundir='./'
subfile="qsubrebound"

#integration settings
integrator="ias15"
t_max=1.e6*2.*np.pi
Noutputs=1000

#restart settings
restartinput="testeject.txt"
#restart from binary file
frombin=False
binfile="checkpoint.txt"

#output settings
checkpoint=True
