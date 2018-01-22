#!/usr/bin/env python
import sys

scan = ""
n = 0

lines = open(sys.argv[1], 'r').readlines()
for l in lines:
  elems = l.split()

  if l[0] == "H":
    continue
  elif l[0] == "S":
    scan = int(elems[1])
    pepm = float(elems[3])
  elif l[0] == "I":
    if elems[1]=="RetTime":
      rt = float(elems[2])*60
  elif l[0] == "Z":
    charge = int(elems[1])
    #output header
    if n>0:
      print "END IONS"
      print
    print "BEGIN IONS"
    print "TITLE=/fake"
    n = n+1
    print "SCANS=%d" % (n)
    print "RTINSECONDS=%f"% (rt)
    print "CHARGE=%d+" % charge
    print "PEPMASS=%f" % pepm
  else:
    print elems[0], elems[1]

print "END IONS"
print
