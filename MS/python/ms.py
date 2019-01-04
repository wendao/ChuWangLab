import math
from scipy import signal
import numpy as np
import atomParams as atom
import peptideParams as pep

ppm = 1e-6

#sp.mass == sp.mz * sp.charge - (sp.charge-1)*ms.mass["Hplus"]
class spectra(object):
  def __init__(self):
    self.charge = 0
    self.scan = 0
    self.pscan = 0  #precursor scan
    self.pint = 0.0 #precursor int
    self.maxint = 0.0
    self.mz = 0.0
    self.mass = 0.0
    self.rt = 0.0
    self.peaks = []

class DTAhit(object):
  def __init__(self):
    self.seq = ""
    self.calc_mass = 0.0
    self.mass = 0.0
    self.scan = ""
    self.prots = []
    self.tol_Int = 0.0
    self.conf = 0.0
    self.deltCN = 0.0
    self.xcorr = 0.0
    self.redundancy = 0
    self.charge = 0

def load_ms(fn):
  with open(fn, 'r') as fp:
    sp = None
    mint = 0.0
    for l in fp:
      es = l.split()
      if l[0] == "S":
        if sp:
          sp.maxint = mint
          yield sp
        sp = spectra()
        mint = 0.0
        sp.scan = int(es[1])
        sp.mz = float(es[3])
      elif l[0] == "I":
        if es[1] == "RetTime":
          sp.rt = float(es[2])
        elif es[1] == "PrecursorInt":
          sp.pint = float(es[2])
        elif es[1] == "PrecursorScan":
          sp.pscan = int(es[2])
      elif l[0] == "Z":
        sp.charge = int(es[1])
        sp.mass = float(es[2])
      elif l[0] == "H":
        pass
      else:
        sp.peaks.append((float(es[0]),float(es[1])))
        if mint < float(es[1]):
          mint = float(es[1])
  #spit the last
  yield sp

#dict1.update(dict2)
def build_ms2_map( fn ):
  mapms = {}
  for sp in load_ms(fn):
    mapms[ fn.split('/')[-1][:-4]+"."+str(sp.scan)+"."+str(sp.scan)+"."+str(sp.charge) ] = sp
  return mapms

def load_dta(fn):
  hits = []
  lines = open(fn, 'r').readlines()
  pro_tags = {}
  pep_tags = {}
  for l in lines:
    elems = l.strip().split('\t')
    N = len(elems)
    if N == 9:
      if elems[0] == "Locus":
        for i, e in enumerate(elems):
          pro_tags[e] = i 
        pro_lst = []
        continue
      else:
        #load protein info
        pro = elems[pro_tags["Locus"]]
        pro_lst.append(pro)
    elif N == 13 or N == 12:
      shift = int(N==12)
      if len(pro_lst)>0:
        current_pro_lst = pro_lst
        pro_lst = []
      if elems[0] == "Unique":
        for i, e in enumerate(elems):
          pep_tags[e] = i 
      else:
        #load peptide info
        hit = DTAhit()
        hit.seq = elems[pep_tags["Sequence"]-shift]
        hit.calc_mass = float(elems[pep_tags["CalcM+H+"]-shift])
        hit.mass = float(elems[pep_tags["M+H+"]-shift])
        hit.conf = float(elems[pep_tags["Conf%"]-shift])
        hit.deltCN = float(elems[pep_tags["DeltCN"]-shift])
        hit.xcorr = float(elems[pep_tags["XCorr"]-shift])
        hit.redundancy = float(elems[pep_tags["Redundancy"]-shift])
        hit.prots = current_pro_lst
        hit.tol_Int = float(elems[pep_tags["TotalIntensity"]-shift])
        hit.scan = elems[pep_tags["FileName"]-shift]
        hit.charge = int(hit.scan.split('.')[-1])
        hits.append(hit)
    else:
      #else line
      continue
  #print pro_tags, pep_tags
  return hits

def nCr(n,r):
  f = math.factorial
  return f(n) / f(r) / f(n-r)

def gaussian(x, alpha, r):
  return np.exp(-alpha*np.power((x - r), 2.))

#counts of each elements
#ave_count(1000)
def ave_count(mass):
    ave_mass = 111.1254
    elements = [ "C", "H", "O", "N" ]
    comp = [ 4.9348, 7.7583, 1.4773, 1.3577 ]
    map_ave_count = {}
    for e, c in zip(elements, comp):
      map_ave_count[e] = int(round(c*mass/ave_mass))
    return map_ave_count

#isotope_dist(ave_count(1000))
def isotope_dist(count_map, N15_enrich=1.0):
    count_local = {}
    for e in atom.mass.keys():
        if e not in count_map.keys():
            count_map[e] = 0
        count_local[e] = count_map[e]
    #print count_local

    atom.heavy_map["N15"] = N15_enrich
    atom.light_map["Se"] = 0.2377 #Special for selenium

    single_prob = {}
    all_prob = []

    for e in atom.mass.keys():
        single_prob[e] = []
        count = count_local[e]
        if count == 0:
            continue
        v = range(count+1)
        l = atom.light_map[e]
        h = atom.heavy_map[e]
        for i in v:
            t1 = round(nCr(count,i)*(l**(count-i))*(h**i),4)
            single_prob[e].append(t1)
        new_prob = single_prob[e]
        if e=="O" or e=="S" or e=="Cl" or e=="Br" or e=="O18":
            new_prob = [0] * (2*count+1)
            for i in xrange(count+1):
                new_prob[i*2] = single_prob[e][i]
            single_prob[e] = new_prob
        if e=="Se":
            single_prob[e] = [0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873]
            new_prob = single_prob[e]
        if ( len(all_prob)>0 ):
            d = len(all_prob) - len(new_prob)
            if d > 0:
                for i in xrange(d): new_prob.append(0)
            elif d < 0:
                for i in xrange(-d): all_prob.append(0)
            all_prob = signal.convolve(np.array(all_prob),np.array(new_prob))
        else:
            all_prob = single_prob[e]
    return all_prob

# generate theory spectra
# for almost all cases "C:57.02146"
def theory_spectra( seq, mods ):
  mod_dict = {}
  dHplus = pep.MonoHplus

  for mod in mods.split("|"):
    es = mod.split(":")
    if len(es) != 2: break
    symbol = es[0]
    d_mass = float(es[1])
    mod_dict[symbol] = d_mass

  mass = 0.0
  an = []
  for c in seq:
    if c not in pep.AA_MonoMassMap.keys():
      if c in mod_dict.keys():
        mass = mass + mod_dict[c]
        if len(an)>0:
          an[-1] += mod_dict[c]
      else:
        print "unrecognized sequence", c
        assert 0
    else:
      mass = mass + pep.AA_MonoMassMap[c]
      if c in mod_dict.keys():
        mass += mod_dict[c]
      an.append(mass)

  mass += (pep.MonoH + pep.MonoOH)

  bp = []
  bpp = []
  bppp = []
  yp = []
  ypp = []
  yppp = []

  for i in range(1, len(an)):
    bp.append(an[i-1]+dHplus)
    yp.append(mass-an[i-1]+dHplus)
    bpp.append((an[i-1]+2*dHplus)/2)
    ypp.append((mass-an[i-1]+2*dHplus)/2)
    bppp.append((an[i-1]+3*dHplus)/3)
    yppp.append((mass-an[i-1]+3*dHplus)/3)

  #[ full-mass, b+, y+, b++, y++ ]
  return (mass+dHplus, bp, yp, bpp, ypp, bppp, yppp)

def logNfactorialSumInt( peaks, preds ):
  #calc hyperscore from sorted list
  #peaks = sorted( peaks, lambda x:x[0] )
  #preds = sorted( preds )

  Nhit = 0
  SumInt = 1.0 #if no any hits return 0
  for peak in peaks:
    tol = 20.0 * ppm * peak[0]
    for mz in preds:
      #no need: #or math.fabs(peak[0]-pep.deltaC12C13-mz)<tol
      if math.fabs(peak[0]-mz)<tol:
        Nhit += 1
        SumInt += peak[1]
        break
  return math.log( math.factorial(Nhit)*SumInt )

def hyperscore( sp, seq, mods ):
  mass, b1, y1, b2, y2, b3, y3 = theory_spectra( seq, mods )
  b_ions = b1 +b2 #+b3
  y_ions = y1 +y2 #+y3
  hs = logNfactorialSumInt( sp, b_ions ) 
  hs += logNfactorialSumInt( sp, y_ions ) 
  return hs, mass 

def scan_hyperscore( seq, diffmass, mods, star, spectra, ref_mass=0.0 ):
  #put diff mass on each position and calculate hyperscore
  #print seq, diffmass

  #filter spectra by intensity and sort by m/z
  #ntop = int(len(spectra.peaks))
  #topeaks = sorted( spectra.peaks, key=lambda x:-x[1] ) #[:ntop]
  topeaks = [ x for x in spectra.peaks if x[1]>spectra.maxint*0.01 ]

  npos = len(seq)
  try_lst = []
  max_hs = 0.0
  for i in xrange(npos):
    if seq[i] in pep.AA_MonoMassMap.keys():
      if i==npos-1 or seq[i+1] in pep.AA_MonoMassMap.keys():
        tryseq = seq[:i+1]+star+seq[i+1:]
        hs, mass = hyperscore( topeaks, tryseq, mods )
        try_lst.append( (tryseq, hs) )
        if ref_mass>0.0:
          print "try:", tryseq, mass, ref_mass
  #print try_lst
  return sorted( try_lst, key=lambda x:-x[1] )

