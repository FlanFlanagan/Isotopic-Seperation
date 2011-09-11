from pyne.origen22 import *
from pyne.material import *
from math import *
import bright
import subprocess
import matplotlib.pyplot as plt
import pickle
import numpy
import sys
from mass_stream import *
from BUd import *
import os


def parse(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     a = Material()
     a.load_from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
  
      ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=tables)
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=[]
     ty=0
     r=parse_tape6('TAPE6.OUT')
     for ii in tables:
	  WattsT.append(r['table_'+str(ii)]['summary']['data'])
          totWT.append(0.0)
          
	  
          isos = set(WattsT[0].keys())
          isos.remove(0)
          for i in range(len(WattsT)):
              for iso in WattsT[i].keys():
                  isos.add(iso)
	  isos.remove(0)
          for iso in isos:
              isoT = []
              for i in range(len(WattsT)):
                  if iso in WattsT[i].keys():
                      isoT.append(WattsT[i][iso])
                      totWT[i] = totWT[i] + WattsT[i][iso][0][1]
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
	  writer2 = open('data_'+str(time)+'.txt','a')
#             print >> 	ii, str(totWT)
	  print >> writer2, '{0} {1}'.format(ii, totWT[ty]/(LWR_BUd+FR_BUd))
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
  ps = MassStream(ns)
  return ps 


def add_mass_stream(ms, added):
    writer3 = open('adding.txt', 'a')
    for iso in ms.comp.keys():
	writer3.write('{0}{1}', format(iso, ms.comp[iso]))
    for iso in added.comp.keys():
	writer3.write('{0}{1}', format(iso, added.comp[iso]))
    writer3.close
    ps = MassStream('added.txt', 1)
    os.remove('adding.txt')
    return ps

def heat(inputfile, time):
     ## Setting up the Mass Stream ##
     a = Material()
     a.load_from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     WattsT =[]
      ## Running Origen ##
     n=0
     t = 0.1
     ty = 0
     el_dictT={}
     totWT=[]
     while t < time*2:
	if t < time:
	    write_tape5_irradiation("IRP", t, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
	    subprocess.call('o2_therm_linux.exe',shell=True)

	    #Watts#
	    r=parse_tape6('TAPE6.OUT')
	    WattsT.append(r['table_'+str(9)]['summary']['data'])
	    totWT.append(0.0)

	    isos = set(WattsT[0].keys())
	    isos.remove(0)
	    for i in range(len(WattsT)):
		for iso in WattsT[i].keys():
		    isos.add(iso)
	    isos.remove(0)
	    for iso in isos:
		isoT = []
		for i in range(len(WattsT)):
		    if iso in WattsT[i].keys():
			isoT.append(WattsT[i][iso])
			totWT[i] = totWT[i] + WattsT[i][iso][0][1]
		    else:
			isoT.append(0.0)

	    #Creating the Table
	    if t <.2:
		writer2 = open('heat_'+str(time)+'.txt','a')
		#mt = (t-0)*(expm1((log1p(totWT[ty]))/2))
		mt = (t-0)*(sqrt(totWT[0]))
		n = n + mt
		print >> writer2, '{0} {1}'.format(t,n)
		writer2.close
	    else:
		writer2 = open('heat_'+str(time)+'.txt','a')
		#ms = (t-(t/2))*(expm1(log1p(totWT[ty] * totWT[ty-1])/2))
		mt = (t-(t/2))*(sqrt(totWT[0] * Prev_tot))
		n = n + mt
		print >> writer2, '{0} {1}'.format(t,n)
		writer2.close
	    os.remove('TAPE6.OUT')
	elif t > time*2:
	    exit
	elif t > time:
	    write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
	    subprocess.call('o2_therm_linux.exe',shell=True)

	    #Watts#
	    r=parse_tape6('TAPE6.OUT')
	    WattsT.append(r['table_'+str(9)]['summary']['data'])
	    totWT.append(0.0)

	    isos = set(WattsT[0].keys())
	    isos.remove(0)
	    for i in range(len(WattsT)):
		for iso in WattsT[i].keys():
		    isos.add(iso)
	    isos.remove(0)
	    for iso in isos:
		isoT = []
		for i in range(len(WattsT)):
		    if iso in WattsT[i].keys():
			isoT.append(WattsT[i][iso])
			totWT[i] = totWT[i] + WattsT[i][iso][0][1]
		    else:
			isoT.append(0.0)
	    writer2 = open('heat_'+str(time)+'.txt','a')
	    #mp = (time-(t/2))*(expm1((log1p(totWT[ty] * totWT[ty-1]))/2))
	    mt = (time-(t/2))*(sqrt(totWT[0] * Prev_tot))
	    n = n + mt
	    print >> writer2, '{0} {1}'.format(time,n/(LWR_BUd+FR_BUd))
	    writer2.close
	Prev_tot= totWT[0]
	ty=ty+1
	t=t*2
	el_dictT={}
	isoT=[]
	isos=[]
	WattsT=[]
	totWT=[]
		#     Storage.calc(instream, time)


