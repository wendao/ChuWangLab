#!/usr/bin/env python
import sys
from math import fabs

# grab mono/+1/+2 peak

fn = sys.argv[1]

dC = 1.0033548
tol = 0.1 #Da
ignore = []
if len(sys.argv)>2:
  for c in sys.argv[2].split(','):
    ignore.append(int(c))

dd_lst = []
n_count = {}
for i in xrange(4):
  n_count[i] = 0

lines = open( fn, 'r' ).readlines()
for l in lines:
  es = l.strip().split()
  dm = float(es[2])
  n = round( dm / dC )
  if n in ignore: continue
  dd = dm - n*dC
  if fabs(dd)<tol:
    dd_lst.append(dd)
    n_count[int(n)] += 1

#print numbers
for dd in dd_lst:
  print dd

import numpy as np
#print n_count
print "# dm=", np.median(dd_lst), "std:", np.std(dd_lst)
for i in xrange(4):
  if n_count[i]>0:
    print "# +"+str(i), n_count[i]


