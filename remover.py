import bright
from mass_stream import *

a={922350: .01, 922360: .01, 942390: .02, 922380: .96}
b=[942390]

nms={}
def remove_(ms, removed):
     isos = set(ms.keys())
     int_ms = list(isos.difference(removed))
     for t in range(len(int_ms)):
	    i=int_ms[t]
	    nms[i]=[]
	    nms[i].append(0.0)
     	    nms[i]=ms[i]

     
	
     return nms	


remove_(a,b)
fms = MassStream(nms)
print fms.comp
	
