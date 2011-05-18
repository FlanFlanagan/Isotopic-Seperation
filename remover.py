import bright
from mass_stream import *



def remove_(ms, removed):
  ns={}
  isos = set(ms.comp.keys())
  int_m = list(isos.difference(removed))
  for t in range(len(int_m)):
	    i=int_m[t]
	    ns[i]=ms.comp[i]
  ms = MassStream(ns)
  return ps.comp