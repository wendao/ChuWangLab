#!/usr/bin/env python
import xml.etree.ElementTree as ET
import sys
from math import log

ns = {'pep': 'http://regis-web.systemsbiology.net/pepXML'}

sc_list = []
perct = float(sys.argv[1])

for fn in sys.argv[2:]:
    try:
        tree = ET.parse(fn)
    except:
        print "Open pepXML file", fn, "failed!"
        sys.exit()

    root = tree.getroot()
    for summary in root.findall('pep:msms_run_summary', ns):
        for query in summary.findall('pep:spectrum_query', ns):
            for result in query.findall('pep:search_result', ns):
                for hit in result.findall('pep:search_hit', ns):
                    #process the hit
                    score_map = {}
                    for sc in hit.findall('pep:search_score', ns ):
                        score_map[sc.attrib['name']]= float(sc.attrib['value'])
                    score = -log(score_map["expect"])
                    pro = hit.attrib['protein'].strip()
                    if pro[:7] == "Reverse":
                        sc_list.append( ("R", score) )
                    else:
                        sc_list.append( ("S", score) )

N_all = len(sc_list)
N_R_tol = int(N_all*perct)
sorted_sc_list = sorted( sc_list, key=lambda x: x[1], reverse=True )

N_Rev = 0
for sc in sorted_sc_list:
    if sc[0] == "R": N_Rev = N_Rev + 1
    if N_Rev >= N_R_tol:
        print sc[1]
        break
