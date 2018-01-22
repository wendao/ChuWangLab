#!/usr/bin/env python
import sys
import ms
from math import fabs, exp

dta = sys.argv[1]
star = sys.argv[2]
diffmass = float(sys.argv[3])
mods = sys.argv[4]
ms2_lst = sys.argv[5:]

map_spectra = {}
for ms2 in ms2_lst:
  map_spectra.update( ms.build_ms2_map(ms2) )

for hit in ms.load_dta(dta):
  #print hit.seq, hit.charge, hit.tol_Int, hit.mass, hit.calc_mass
  realseq = hit.seq.split('.')[1]
  if star in realseq:
    np = 0
    while np > -1:
      np = realseq.find( star, np )
      if np == -1: break

      print realseq, hit.scan
      rawseq = realseq[:np]+realseq[np+1:]
      sp = map_spectra[hit.scan]
      #check
      calc_mass = ms.theory_spectra( realseq, mods )[0]
      topseq = ms.scan_hyperscore( rawseq, diffmass, mods, star, sp )
      #print calc_mass, hit.calc_mass, "diff=", fabs(calc_mass-hit.calc_mass)
      #print topseq
      maxhs = topseq[0][1]
      for n, ( seq, hs ) in enumerate( topseq ):
        if hs < maxhs:
          sechs = hs
          break
      results = []
      sum_w = 0.0
      max_w = 0.0
      for n, ( seq, hs ) in enumerate( topseq ):
        delta_hs = hs - sechs
        weight = exp(delta_hs)
        if weight > max_w: max_w = weight
        results.append([seq, weight])
        sum_w += weight
      best = []
      for result in results:
        w = result[1]
        p = w / sum_w
        if p>0.01: print "--", result[0], p
        if w>max_w-0.01: best.append(result[0])
      if realseq in best: print "HIT: True", len(best)
      else: print "HIT: False"
      print
      #next
      np = np + 1

