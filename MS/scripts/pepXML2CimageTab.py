#!/usr/bin/env python
import xml.etree.ElementTree as ET
from math import fabs,log
import sys

ns = {'pep': 'http://regis-web.systemsbiology.net/pepXML'}

mass_tol = 0.01

mass_norm_AA = {
"G" : 57.02146,
"A" : 71.03711,
"S" : 87.03203,
"P" : 97.05276,
"V" : 99.06841,
"T" :101.04768,
"C" :103.00918,
"L" :113.08406,
"I" :113.08406,
"N" :114.04293,
"D" :115.02694,
"Q" :128.05858,
"K" :128.09496,
"E" :129.04259,
"M" :131.04048,
"H" :137.05891,
"F" :147.06841,
"R" :156.10111,
"Y" :163.06333,
"W" :186.07931
}

mass_nPter_mod = 42.01060
mass_C_blocker = 57.02416

sc_cut = float(sys.argv[1])
modstr = sys.argv[2]
style = sys.argv[3]

# parsing the modification description string
# elements split by |
# each element has 3/4 components split by ,
# AA, diffmass, [H/L], [marker]
mod_lst = []
std_lst = {}
for mod in modstr.split('|'):
    elems = mod.split(',')
    AA = elems[0]
    mass = float(elems[1])
    if AA=="C": mass = mass + mass_C_blocker
    tag = elems[2]
    if mass == 0:
      std_lst[AA] = tag
    else:
      if len(elems)>3:
        mod_lst.append((AA,mass,tag,elems[3]))
      else:
        mod_lst.append((AA,mass,tag))

ipi_lst = {}
scan_lst = []
cross_tab = {}

def get_label( labels, seq ):
  # label H/L based on modified pos
  # return all H or all L, otherwise NULL
  for n,a in enumerate(seq):
    # has aa
    if a in std_lst.keys():
      # not labeled
      if n not in labels.keys():
        labels[n]=std_lst[a]
  uniq_label = []
  for label in labels.values():
    if label not in uniq_label: uniq_label.append(label)
  if len(uniq_label) == 0:
    if style in ["pro", "pos"]:
      return "light"
    else:
      return "NULL"
  if len(uniq_label) == 1:
    return uniq_label[0]
  return "NULL"

def mark_seq(uAA, seq, dAA, markers):
  #print markers
  ans = ""
  if len(markers)==0:
    ans = uAA + "." + seq + "." + dAA
  else:
    #all_seq = uAA + "." + seq[:ndx+1] + "*" + seq[ndx+1:] + "." + dAA
    ans = uAA + "."
    prev = 0
    for marker, ndx in markers:
      ans = ans + seq[prev:ndx+1] + marker
      prev = ndx+1
    ans = ans + seq[prev:] + "." + dAA
  #print ans
  return ans

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

  labels = {}
  markers = []
  for mod_info in hit.findall('pep:modification_info', ns):
    for mod in mod_info.findall('pep:mod_aminoacid_mass', ns):
      ndx = int(mod.attrib['position'])-1
      mass = float(mod.attrib['mass'])
      ideal_mass = mass_norm_AA[seq[ndx]]
      diff_mass = mass - ideal_mass
      #print seq[ndx], diff_mass
      for m in mod_lst:
          if m[0] == "n" and ndx == 0: #nTerm
            #ideal_mass = mass_norm_AA[seq[0]]
            if fabs(diff_mass-m[1])<mass_tol or ( uAA=="-" and fabs(diff_mass-mass_nPter_mod-m[1])<mass_tol ):
              #print "nTerm"
              if len(m[2])>0: labels[0]=m[2]
              if len(m)>3: markers.append((m[3],0))
          elif m[0] == seq[ndx] and fabs(diff_mass-m[1])<mass_tol: #normal
            #print "K-heavy"
            if len(m[2])>0: labels[ndx]=m[2]
            if len(m)>3: markers.append((m[3],ndx))

  if style in ["pep", "pos"]:
    if len(markers)==0:
      return
  # check labels
  label = get_label( labels, seq )
  if label != "NULL":
    # add protein
    if pro not in ipi_lst.keys():
      ipi_lst[pro] = 1
    else:
      ipi_lst[pro] = ipi_lst[pro]+1
    gene = pro.split()[0].split('|')[1]
    if pro[:7] == "Reverse": gene = "Reverse_" + gene
    # add scan
    all_seq = mark_seq(uAA, seq, dAA, markers)
    #print all_seq
    scan_key = gene + ":" + all_seq + ":" + str(charge) + ":" + fn_ndx
    scan_lst.append(scan_key + " " + com_fn + " " + str(scan_id) + " " + label + "\n")
    # add cross
    if scan_key not in cross_tab.keys():
        cross_tab[scan_key] = []
    cross_tab[scan_key].append((neutral_mass, scan_id, score))

def get_symbol(disc):
  gene = ""
  elems = disc.split()
  if "GN=" in disc:
    #get from record
    for elem in elems:
      if elem[:3] == "GN=":
        gene = elem[3:]
  else:
    gene = elems[0]
  return gene

#loop all pepXML file
for fn in sys.argv[4:]:
    raw_fn = fn.split('/')[-1]
    elems = raw_fn.split('_')
    com_fn = ""
    for elem in elems[:-1]:
        com_fn = com_fn + elem + "_"
    com_fn = com_fn[:-1]
    fn_ndx = elems[-1].split(".")[0]

    try:
      tree = ET.parse(fn)
    except:
      print "Open pepXML file", fn, "failed!"
      sys.exit()

    root = tree.getroot()
    for summary in root.findall('pep:msms_run_summary', ns):
      for query in summary.findall('pep:spectrum_query', ns):
        #print query.attrib['spectrum'], query.attrib['index'], query.attrib['assumed_charge'], query.attrib['precursor_neutral_mass']
        scan_id = int(query.attrib['spectrum'].split('|')[-1])
        mass = float(query.attrib['precursor_neutral_mass'])
        rt = float(query.attrib['retention_time_sec'])
        charge = int(query.attrib['assumed_charge'])
        for result in query.findall('pep:search_result', ns):
          for hit in result.findall('pep:search_hit', ns):
            process_scan(hit, scan_id, mass, rt, charge)

#save ipi_name.table
ipiout = open("ipi_name.table", 'w')
ipiout.write("name\n")

if style == "pro":
  out_ipi_lst = [ p for p in ipi_lst.keys() if ipi_lst[p]>=2 ]
else:
  out_ipi_lst = [ p for p in ipi_lst.keys() if ipi_lst[p]>=1 ]

for pro in out_ipi_lst:
    #skip rev decoys ?
    elems = pro.split('|')
    ipi = elems[1] #uniprot actually
    disc = elems[2]
    tmp = disc.split()
    gene = tmp[0]
    disc = disc[len(gene)+1:]
    symbol = get_symbol(disc)
    if pro[:7] == "Reverse":
        ipi = "Reverse_" + ipi
    disc = disc.split("[REVIEWED]")[0]
    label = disc.find("[NOT_REVIEWED]")
    if label>0:
         disc = disc[:label]
    elems = disc.split()
    disc = ""
    for elem in elems:
      if "=" in elem: break
      disc = disc + " " + elem
    disc = disc[:min(60,len(disc))]
    disc = disc.replace("'", "")
    #ipiout.write(pro)
    #ipiout.write("\n")
    ipiout.write(ipi)
    ipiout.write("\t")
    ipiout.write(symbol)
    ipiout.write(" ")
    ipiout.write(disc)
    ipiout.write("\n")
ipiout.close()

#save all_scan.table
scanout = open("all_scan.table", 'w')
scanout.write("key run scan HL\n")
for scan in scan_lst:
    scanout.write(scan)
scanout.close()

#save cross_scan.table
crossout = open("cross_scan.table", 'w')
crossout.write( "key mass %s\n" % (com_fn) )
for scan_core in cross_tab.keys():
    rank = sorted( cross_tab[scan_core], key=lambda x: x[2] )
    neutral_mass = rank[-1][0]
    id_scan = rank[-1][1]
    crossout.write(scan_core+" "+str(neutral_mass)+" "+str(id_scan)+"\n")
crossout.close()
