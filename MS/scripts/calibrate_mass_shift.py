#!/usr/bin/env python
import sys
from math import fabs

#grab protein and peptide info from diffmass.txt

fn = sys.argv[1]
pro = sys.argv[2]
pep = sys.argv[3]
ref = 0.0
if len(sys.argv)>4:
  ref = float(sys.argv[4])

dC = 1.0033548
tol = 0.1 #Da

dd_lst = []
n_count = {}
for i in xrange(4):
  n_count[i] = 0

lines = open( fn, 'r' ).readlines()
for l in lines:
  es = l.strip().split()
  if pro in es[0] and pep == es[1]:
    dm = float(es[2])
    n = round( dm / dC )
    dd = dm - n*dC
    if fabs(dd)<tol:
      dd_lst.append(dd)
      n_count[int(n)] += 1

import numpy as np
#print n_count
print "# dm=", ref, np.median(dd_lst), "std:", np.std(dd_lst)
for i in xrange(4):
  if n_count[i]>0:
    print "# +"+str(i), n_count[i]

if ref<-999:
  for dd in dd_lst:
    print dd

