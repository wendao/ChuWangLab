#!/usr/bin/env python
import xml.etree.ElementTree as ET
from math import fabs, log
import sys

ns = {'pep': 'http://regis-web.systemsbiology.net/pepXML'}

sc_cut = float(sys.argv[1])
modstr = sys.argv[2]

mod_lst = []
for mod in modstr.split('|'):
    elems = mod.split(',')
    AA = elems[0]
    mass = float(elems[1])
    tag = elems[2]
    mod_lst.append((AA,mass,tag))

def process_scan( hit, scan_id, neutral_mass, rt, charge ):
  score_map = {}
  for sc in hit.findall('pep:search_score', ns ):
    score_map[sc.attrib['name']]= float(sc.attrib['value'])

  #score = score_map["hyperscore"]
  score = -log(score_map["expect"])
  if score < sc_cut: return

  #good scan
  pro = hit.attrib['protein'].strip()
  seq = hit.attrib['peptide']
  pep_mass = float( hit.attrib['calc_neutral_pep_mass'] )
  uAA = hit.attrib['peptide_prev_aa']
  dAA = hit.attrib['peptide_next_aa']
  num_tryp_ter = int(hit.attrib['num_tol_term'])
  #num_tot_proteins
  rank = int( hit.attrib['hit_rank'] )
  tot_num_ions = int( hit.attrib['tot_num_ions'] )
  num_matched_ions = int( hit.attrib['num_matched_ions'] )
  delta_mass = float(hit.attrib['massdiff'])
  num_miss_clv = int( hit.attrib['num_missed_cleavages'] )

  mod_pos = []
  for mod_info in hit.findall('pep:modification_info', ns):
    for mod in mod_info.findall('pep:mod_aminoacid_mass', ns):
      ndx = int(mod.attrib['position'])-1
      mass = float(mod.attrib['mass'])
      for m in mod_lst:
          if m[0] == seq[ndx] and fabs(m[1]-mass)<0.01:
              #save pos
              mod_pos.append(ndx)
  #output
  if len(mod_pos)>0:
    for n, pos in enumerate(mod_pos):
      mark_seq = seq[:pos+n+1] + "*" + seq[pos+n+1:]
    print pro.split()[0], uAA+"."+mark_seq+"."+dAA, delta_mass

#loop all pepXML file
for fn in sys.argv[3:]:
    raw_fn = fn.split('/')[-1]
    com_fn = raw_fn[:-10]
    fn_ndx = raw_fn[-9:-7]

    try:
      tree = ET.parse(fn)
    except:
      print "Open pepXML file", fn, "failed!"
      sys.exit()

    root = tree.getroot()
    for summary in root.findall('pep:msms_run_summary', ns):
      for query in summary.findall('pep:spectrum_query', ns):
        #print query.attrib['spectrum'], query.attrib['index'], query.attrib['assumed_charge'], query.attrib['precursor_neutral_mass']
        scan_id = 0
        mass = float(query.attrib['precursor_neutral_mass'])
        rt = float(query.attrib['retention_time_sec'])
        charge = int(query.attrib['assumed_charge'])
        for result in query.findall('pep:search_result', ns):
          for hit in result.findall('pep:search_hit', ns):
            process_scan(hit, scan_id, mass, rt, charge)

