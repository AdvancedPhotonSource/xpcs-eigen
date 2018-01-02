#!/bin/env python

import h5py
import numpy as np
import sys

f = h5py.File(sys.argv[1])

sg21 = f['/original/G2']
sip1 = f['/original/IP']
sif1 = f['/original/IF']


sg22 = f['/exchange/G2']
sip2 = f['/exchange/IP']
sif2 = f['/exchange/IF']


g21 = np.zeros(sg21.shape, dtype='float32')
ip1 = np.zeros(sip1.shape, dtype='float32')
if1 = np.zeros(sif1.shape, dtype='float32')

sg21.read_direct(g21)
sip1.read_direct(ip1)
sif1.read_direct(if1)

g22 = np.zeros(sg22.shape, dtype='float32')
ip2 = np.zeros(sip2.shape, dtype='float32')
if2 = np.zeros(sif2.shape, dtype='float32')

sg22.read_direct(g22)
sip2.read_direct(ip2)
sif2.read_direct(if2)
