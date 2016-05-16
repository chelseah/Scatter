#!/usr/bin/env python
from setting import *
from submit_rebound import submit
import os
def main():
    #nnodes=20
    for i in xrange(nnodes):
        print i
	#subfile=basename+"optimize_%d" % i
	subfile=basename+"qsubrebound_%d" % i
        submit(basename,"qsubrebound_%d" % i,start=1+8*i)
        os.system('qsub %s' % subfile)

if __name__=='__main__':
    main()
