#!/usr/bin/env python
import sys

for fn in sys.argv[1:]:
    elems = fn.split("/")
    if len(elems)==1: continue
    newfn = elems[-1]
    print "writing", newfn, "..."
    lines = open(fn, 'r').readlines()
    newfp = open(newfn, 'w')
    for l in lines:
        if l[:6] == "TITLE=":
            title = l.strip()
        elif l[:6] == "SCANS=":
            scan = l.split("=")[1]
            newfp.write(title+"|"+scan)
            newfp.write(l)
        else:
            newfp.write(l)
    newfp.close()

print "done!"
