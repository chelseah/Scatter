#!/usr/bin/env python
from setting import *
from submit_rebound import submit

def main():
    for i in xrange(nnodes):
        print i
        submit(basename,"qsubrebound_%d" % i,start=1+8*i)
        os.system('qsub %s' % subfile)

if __name__=='__main__':
    main()
