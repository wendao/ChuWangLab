#!/usr/bin/env python
import sys
import numpy as np

diffm = []
lines = open( sys.argv[1], 'r' ).readlines()
for line in lines:
    elems = line.split()
    diffm.append(float(elems[2]))

sortm = sorted(diffm)
ndat = len(sortm)
nmid = 50
nwin = 101
dtol = 0.15

#for i in xrange(ndat-nwin):
#    dwin = sortm[i:i+nwin]
#    if dwin[-1]-dwin[0]<dtol:
#        print dwin[nmid]
#print sortm
ps = 0
pe = nwin
while ps < ndat and pe < ndat:
  if sortm[pe-1]-sortm[ps]>dtol:
    ps = ps + 1
    pe = pe + 1
    if pe == ndat: break
  else:
    dmid = np.median( sortm[ps:pe] )
    while sortm[pe-1]-dmid < dtol and pe < ndat:
      #extend window
      pe = pe + 1
      dmid = np.median( sortm[ps:pe] )
    print dmid, 0
    print dmid, pe - ps
    print
    ps = pe
    pe = ps + nwin
