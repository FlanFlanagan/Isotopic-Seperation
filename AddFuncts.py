from metasci.nuke import origen as msno
from mass_stream import *
import bright
import subprocess
import matplotlib.pyplot as plt
import pickle
import numpy

def parse(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     instream=inputfile
     a = MassStream(instream, 1)
     a.comp
     msno.write_tape4(a.comp, name='TAPE4.INP')
  
      ## Running Origen ##
     msno.write_tape5_irradiation("IRP", time , 0, (204, 205, 206), name='TAPE5.INP', out_table_num=tables)
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=[]
     ty=0
     for ii in tables:
	  r=msno.parse_tape6(name='TAPE6.OUT')
          WattsT.append(r['table_'+str(ii)]['summary']['data'])
          totWT.append(0.0)

          isos = set(WattsT[0].keys())
          for i in range(len(WattsT)):
              for iso in WattsT[i].keys():
                  isos.add(iso)

          for iso in isos:
              isoT = []
              for i in range(len(WattsT)):
                  if iso in WattsT[i].keys():
                      isoT.append(WattsT[i][iso])
                      totWT[i] = totWT[i] + WattsT[i][iso]
                  else:
                      isoT.append(0.0)
    
      #plt.plot( range(1,20), isoT, label=str(iso))
#          for t in range(len(WattsT)):
#               for key in WattsT[t].keys():
#                   el=key/10000
#                   if el in el_dictT.keys():
#                       el_dictT[el][t]=el_dictT[el][t]+WattsT[t][key]
#                   else:
#                       el_dictT[el]=[]
#                       for tt in range(len(WattsT)):
#                          el_dictT[el].append(0.0)
#                       el_dictT[el][t]=el_dictT[el][t]+WattsT[t][key]
#          l=[(max(el_dictT[el]), el) for el in el_dictT.keys()]
#          l.sort(reverse=True)
#          tlt=[totWT, el_dictT[l[0][1]], l[0][1]]

     #Creating the Table
	  writer2 = open('data.txt','a')
#             print >> 	ii, str(totWT)
	  print >> writer2, '{0} {1}'.format(ii, totWT[ty])
	  writer2.close
	  ty=ty+1
	  #     Storage.calc(instream, time)

def remove(ms, removed):
  ns={}
  isos = set(ms.comp.keys())
  int_m = list(isos.difference(removed))
  for t in range(len(int_m)):
	    i=int_m[t]
	    ns[i]=ms.comp[i]
  ms = MassStream(ns)
  return ps.comp


def heat(inputfile, time):
     ## Setting up the Mass Stream ##
     a = MS.MassStream(inputfile, 1)
     a.comp
     msno.write_tape4(a.comp, name='TAPE4.INP')
     WattsT =[]
      ## Running Origen ##
     t = 0.1
     ty = 0
     el_dictT={}
     totWT=[]
     while t < time:
	msno.write_tape5_irradiation("IRP", t , 0, (204, 205, 206), name='TAPE5.INP', out_table_num=9)
	subprocess.call('o2_therm_linux.exe',shell=True)

	#Watts#
	r=msno.parse_tape6(name='TAPE6.OUT')
	WattsT.append(r['table_'+str(9)]['summary']['data'])
	totWT.append(0.0)

	isos = set(WattsT[0].keys())
	for i in range(len(WattsT)):
	    for iso in WattsT[i].keys():
		isos.add(iso)

	for iso in isos:
	    isoT = []
	    for i in range(len(WattsT)):
		if iso in WattsT[i].keys():
		    isoT.append(WattsT[i][iso])
		    totWT[i] = totWT[i] + WattsT[i][iso]
		else:
		    isoT.append(0.0)

	#Creating the Table
	if t <.1:
	  writer2 = open('heat.txt','a')
	  mt = ((t+0)/2)*exp((log1p(totWT[ty])+0)/2)
	  print >> writer2, '{0} {1}'.format(t,mt)
	  writer2.close
	else:
	  writer2 = open('heat.txt','a')
	  ms = (((t+t/2)/2)-(t/2-t/4)/2)*exp((log1p(totWT[ty])+log1p(totWT[ty-1])/2))
	  print >> writer2, '{0} {1}'.format(t,ms)
	  writer2.close
	
	ty=ty+1
	t = t*2
	#     Storage.calc(instream, time)


