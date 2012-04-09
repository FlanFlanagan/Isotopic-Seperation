from pyne.origen22 import *
from pyne.material import *
from math import *
import bright
import subprocess
import matplotlib.pyplot as plt
import pickle
import numpy
import sys
import os
from BUd import *


totals =['H','HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR',
	'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD',
	'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER',
	'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TI', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',
	'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'HA', 'SG', 'NS', 'HS', 'MT']

def fuelout(inputfile):
    diff = [922350, 922380, 962420, 962430, 962440, 962450, 962460, 962470,962480, 952410, 952430, 942380, 942390, 942400, 942410, 942420, 400930, 430990, 461070, 340790,531290, 551350]
    a = from_text(inputfile)
    writer1 = open('fuelcomps_'+str(inputfile), 'a')
    b = {}
    for iso in diff:
	if iso in a.keys():
	    writer1.write(str(iso)+' '+str(a[iso])+'\n')
    writer1.close()
	
    
	
def ansource(inputfile, time):
     ## Setting up the Mass Stream ##
     a = from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     writer1 = open('ans_'+str(time)+'.txt', 'a')
     writer1.write('\n'+'Standard_')
     ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[1])
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=0
     r=parse_tape6('TAPE6.OUT')
     WattsT.append(r['alpha_neutron_source'])
     isos = set(WattsT[0].keys())
     c = ['units', 'title']
     isos = isos.difference(c)
     
     for iso in isos:
	 totWT = totWT + WattsT[0][iso][1]
     writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
     writer1.close()

def sfsource(inputfile, time):
     ## Setting up the Mass Stream ##
     a = from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     disos = os.getcwd()
     writer1 = open('sfs_'+str(time)+'.txt', 'a')
     ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[1])
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=0
     r=parse_tape6('TAPE6.OUT')
     WattsT.append(r['spont_fiss_neutron_source'])
     isos = set(WattsT[0].keys())
     c = ['units', 'title']
     isos = isos.difference(c)
     
     for iso in isos:
	 totWT = totWT + WattsT[0][iso][1]
     writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
     writer1.close()
	     
def sfsmax(inputfile, number, time):
    a = from_text(inputfile)
    write_tape4(a, 'TAPE4.OUT')
    write_tape5_irradiation('IRP', time, 0, "TAPE5.INP", (1,2,3), (204,205,206), out_table_num=[1])
    subprocess.call('o2_therm_linux.exe', shell = True)
    r=parse_tape6('TAPE6.OUT')
    watts=[]
    dict={}

    watts.append(r['spont_fiss_neutron_source'])
    isos = set(watts[0].keys())
    isos = isos.difference(totals)
    extra = ['units','title']
    isos = isos.difference(extra)
    for i in range(len(watts)):
	for iso in isos:
	    if iso in watts[i].keys():
		dict[iso] = watts[i][iso][1]
    
    b={}
    c = dict
    d = dict
    n = 0
    x = 0
    writer5 = open('max.txt', 'a')
    for it in c:
	x = x + c[it]
    writer5.write('sfs_'+str(time)+' '+str(x)+'\n')
    while n < number:
	em={}
	for itt in d.keys():
	    em[d[itt]]=itt
	b[em[max(em)]] = max(em)
	f=[em[max(em)]]
	isos2 = set(d.keys())
	isos2 = isos2.difference(f)
	d = {}
	for iiso in isos2:
	    d[iiso] = c[iiso]
	writer5.write(str(em[max(em)])+'_'+str(max(em))+'\n')
	n=n+1
	isos2=[]
    return b
	     
def parse_fp(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     a = Material()
     a.from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     writer1 = open('fp_data_'+str(time)+'.txt', 'a')
     if len(tables) > 1:
	writer1.write('\n'+str(RmIsos)+'_')

      ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=tables)
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=0
     ty=0
     r=parse_tape6('TAPE6.OUT')
     for ii in tables:
	  WattsT.append(r['table_'+str(ii)]['summary']['activation_products'])
	                     	  
          isos = set(WattsT[0].keys())
	  isos = isos.difference(totals)
	  
	  for iso in isos:
	      if iso in WattsT[0].keys():
		  totWT = totWT + WattsT[0][iso][1]
	      

	  if time < 10:
	      if len(tables)==1:
		  writer1.write(str(ii)+' '+str(totWT)+'\n')
	      if len(tables)>1:
		  writer1.write(str(totWT)+'__')
	  if time > 11:
	      if len(tables)==1:
		  writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
	      if len(tables)>1:
		  writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'__')
	  totWT=0
	  ty=ty+1
     
def parse(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     a = Material()
     a.from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     writer1 = open('data_'+str(time)+'.txt', 'a')
     if len(tables) > 1:
	writer1.write('\n'+'Standard_')
     ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=tables)
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=0
     ty=0
     r=parse_tape6('TAPE6.OUT')
     for ii in tables:
	  WattsT.append(r['table_'+str(ii)]['summary']['actinides'])
	  WattsT.append(r['table_'+str(ii)]['summary']['fission_products'])
	  WattsT.append(r['table_'+str(ii)]['summary']['activation_products'])
                    	  
          isos = set(WattsT[0].keys())
	  isos.update(set(WattsT[1].keys()), set(WattsT[2].keys()))
	  isos = isos.difference(totals)
	  
	  for iii in range(len(WattsT)):
	      for iso in isos:
		  if iso in WattsT[iii].keys():
		      totWT = totWT + WattsT[iii][iso][1]
	      
	  
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
	  if len(tables)==1:
	      writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
	      writer1.write(str(totWT)+'\n')
	  if len(tables)>1:
	      writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'__')
	      writer1.write(str(totWT)+'__')
	  totWT=0
	  ty=ty+1
	  
def parse_act(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     a = Material()
     a.from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     writer1 = open('act_data_'+str(time)+'.txt', 'a')
     if len(tables) > 1:
	writer1.write('\n'+str(RmIsos)+'_')
     ## Running Origen ##
     write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=tables)
     subprocess.call('o2_therm_linux.exe',shell=True)

      #Watts#
     WattsT=[]
     el_dictT={}
     totWT=0
     ty=0
     r=parse_tape6('TAPE6.OUT')
     for ii in tables:
	  WattsT.append(r['table_'+str(ii)]['summary']['actinides'])
                    	  
          isos = set(WattsT[0].keys())
	  isos=isos.difference(totals)
	  
	  for iso in isos:
	    if iso in WattsT[0].keys():
		totWT = totWT + WattsT[0][iso][1]
     if len(tables)==1:
	  writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
	  writer1.write(str(totWT)+'\n')
     if len(tables)>1:
	  writer1.write(str(totWT/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'__')
	  writer1.write(str(totWT)+'__')
     totWT=0
     ty=ty+1
	  
	  #     Storage.calc(instream, time)
def dictmax_fp(inputfile, num, time, tables):
    a = from_text(inputfile)
    write_tape4(a, 'TAPE4.OUT')
    write_tape5_irradiation('IRP', time, 0, "TAPE5.INP", (1,2,3), (204,205,206), out_table_num=tables)
    subprocess.call('o2_therm_linux.exe', shell = True)
    r=parse_tape6('TAPE6.OUT')
    watts=[]
    dict={}
    for iii in tables:
	watts.append(r['table_'+str(iii)]['summary']['activation_products'])
	isos = set(watts[0].keys())
	isos = isos.difference(totals)
	for iso in isos:
	    if iso in watts[0].keys():
		dict[iso] = watts[0][iso][1]
	b={}
	c = dict
	d = dict
	n = 0
	x = 0
	writer5 = open('max.txt', 'a')
	for it in c:
	    x = x+c[it]
	writer5.write('table_'+str(iii)+'_'+str(time)+' '+str(x)+'\n')
	while n < num:
	    em={}
	    for itt in d.keys():
		em[d[itt]]=itt
	    b[em[max(em)]] = max(em)
	    f=[em[max(em)]]
	    isos2 = set(d.keys())
	    isos2 = isos2.difference(f)
	    d = {}
	    for iiso in isos2:
		d[iiso] = c[iiso]
	    writer5.write(str(em[max(em)])+'_'+str(max(em))+'\n')
	    n=n+1
	    isos2=[]
	watts=[]
	  
def dictmax(inputfile, num, time, tables):
    a = from_text(inputfile)
    write_tape4(a, 'TAPE4.OUT')
    write_tape5_irradiation('IRP', time, 0, "TAPE5.INP", (1,2,3), (204,205,206), out_table_num=tables)
    subprocess.call('o2_therm_linux.exe', shell = True)
    r=parse_tape6('TAPE6.OUT')
    watts=[]
    dict={}
    for iii in tables:
	watts.append(r['table_'+str(iii)]['summary']['actinides'])
	watts.append(r['table_'+str(iii)]['summary']['fission_products'])
	watts.append(r['table_'+str(iii)]['summary']['activation_products'])
	isos = set(watts[0].keys())
	isos.update(set(watts[1].keys()), set(watts[2].keys()))
	isos = isos.difference(totals)
	for i in range(len(watts)):
	    for iso in isos:
		if iso in watts[i].keys():
		    dict[iso] = watts[i][iso][1]
	b={}
	c = dict
	d = dict
	n = 0
	x = 0
	writer5 = open('max.txt', 'a')
	for it in c:
	    x = x+c[it]
	writer5.write('table_'+str(iii)+'_'+str(time)+' '+str(x)+'\n')
	while n < num:
	    em={}
	    for itt in d.keys():
		em[d[itt]]=itt
	    b[em[max(em)]] = max(em)
	    f=[em[max(em)]]
	    isos2 = set(d.keys())
	    isos2 = isos2.difference(f)
	    d = {}
	    for iiso in isos2:
		d[iiso] = c[iiso]
	    writer5.write(str(em[max(em)])+'_'+str(max(em))+'\n')
	    n=n+1
	    isos2=[]
	watts=[]
	  
	  
def remove(ms, removed):
  ns={}
  isos = set(ms.comp.keys())
  int_m = list(isos.difference(removed))
  for t in range(len(int_m)):
	    i=int_m[t]
	    ns[i]=ms.comp[i]
  ps = Material(ns)
  return ps 



def heat(inputfile, time):
     ## Setting up the Mass Stream ##
     a = Material()
     a.from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     WattsT =[]
      ## Running Origen ##
     n=0
     t = 0.1
     el_dictT={}
     totWT=0
     writer1 = open('heat_total_'+str(time)+'.txt', 'a')
     while t < time*3:
	if t < time:
	    write_tape5_irradiation("IRP", t, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
	    subprocess.call('o2_therm_linux.exe',shell=True)

	    #Watts#
	    r=parse_tape6('TAPE6.OUT')
	    WattsT.append(r['table_'+str(9)]['summary']['actinides'])
	    WattsT.append(r['table_'+str(9)]['summary']['fission_products'])
	    WattsT.append(r['table_'+str(9)]['summary']['activation_products'])
			    
	    isos = set(WattsT[0].keys())
	    isos.update(set(WattsT[1].keys()), set(WattsT[2].keys()))
	    isos = isos.difference(totals)
	    
	    for iii in range(len(WattsT)):
		for iso in isos:
		    if iso in WattsT[iii].keys():
			totWT = totWT + WattsT[iii][iso][1]
	    #Creating the Table
	    if t <.2:
		writer2 = open('heat_'+str(time)+'.txt','a')
		#mt = (t-0)*(expm1((log1p(totWT[ty]))/2))
		mt = (t-0)*(sqrt(totWT))
		n = n + mt
		print >> writer2, '{0} {1}'.format(t,n)
		writer2.close()
	    else:
		writer2 = open('heat_'+str(time)+'.txt','a')
		#ms = (t-(t/2))*(expm1(log1p(totWT[ty] * totWT[ty-1])/2))
		mt = (t-(t/2))*(sqrt(totWT * Prev_tot))
		n = n + mt
		print >> writer2, '{0} {1}'.format(t,n)
		writer2.close()
	    os.remove('TAPE6.OUT')
	elif t > time*2:
	    break
	elif t > time:
	    write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
	    subprocess.call('o2_therm_linux.exe',shell=True)

	    #Watts#
	    r=parse_tape6('TAPE6.OUT')
	    WattsT.append(r['table_'+str(9)]['summary']['actinides'])
	    WattsT.append(r['table_'+str(9)]['summary']['fission_products'])
	    WattsT.append(r['table_'+str(9)]['summary']['activation_products'])
			    
	    isos = set(WattsT[0].keys())
	    isos.update(set(WattsT[1].keys()), set(WattsT[2].keys()))
	    isos = isos.difference(totals)
	    
	    for iii in range(len(WattsT)):
		for iso in isos:
		    if iso in WattsT[iii].keys():
			totWT = totWT + WattsT[iii][iso][1]
			
	    mt = (time-(t/2))*(sqrt(totWT * Prev_tot))
	    n = n + mt
			
	    writer1.write(str(time)+' '+str(n/(LWR_BUd*(.33*1000/24)+FR_BUd*(1000*.42/24)))+'\n')
	    writer1.close()
	Prev_tot= totWT
	t=t*2
	el_dictT={}
	isos=[]
	WattsT=[]
	totWT=0
		#     Storage.calc(instream, time)


def heatmax(inputfile, time, num):
     ## Setting up the Mass Stream ##
     a = Material()
     a.from_text(inputfile)
     write_tape4(a, outfile='TAPE4.INP')
     WattsT =[]
      ## Running Origen ##
     n=0
     dict={}
     dict2={}
     t = 0.1
     el_dictT={}
     totWT=0
     writer1 = open('heat_total_max'+str(time)+'.txt', 'a')
     write_tape5_irradiation("IRP", t, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
     subprocess.call('o2_therm_linux.exe',shell=True)

     r=parse_tape6('TAPE6.OUT')
     WattsT.append(r['table_'+str(9)]['summary']['actinides'])
     WattsT.append(r['table_'+str(9)]['summary']['fission_products'])
     WattsT.append(r['table_'+str(9)]['summary']['activation_products'])
		    
     isos = set(WattsT[0].keys())
     isos.update(set(WattsT[1].keys()), set(WattsT[2].keys()))
     isos = isos.difference(totals)
     for i in range(len(WattsT)):
	 for iso in isos:
	     if iso in WattsT[i].keys():
	  	 dict2[iso] = WattsT[i][iso][1]
     for iso2 in dict2.keys():
	while t < time*3:
	    if t < time:
		write_tape5_irradiation("IRP", t, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
		subprocess.call('o2_therm_linux.exe',shell=True)

		#Watts#
		r=parse_tape6('TAPE6.OUT')
		WattsT.append(r['table_'+str(9)]['summary']['actinides'])
		WattsT.append(r['table_'+str(9)]['summary']['fission_products'])
		WattsT.append(r['table_'+str(9)]['summary']['activation_products'])
				
		for iv in range(len(WattsT)):
		    if iso2 in WattsT[iv].keys():
			totWT = WattsT[iv][iso2][1]
		#Creating the Table
		if t <.2:
		    writer2 = open('heat_'+str(time)+'.txt','a')
		    #mt = (t-0)*(expm1((log1p(totWT[ty]))/2))
		    mt = (t-0)*(sqrt(totWT))
		    n = n + mt
		else:
		    writer2 = open('heat_'+str(time)+'.txt','a')
		    #ms = (t-(t/2))*(expm1(log1p(totWT[ty] * totWT[ty-1])/2))
		    mt = (t-(t/2))*(sqrt(totWT * Prev_tot))
		    n = n + mt
		    
		os.remove('TAPE6.OUT')
	    elif t > time:
		write_tape5_irradiation("IRP", time, 0, "TAPE5.INP", (1,2,3), (204, 205, 206), out_table_num=[9])
		subprocess.call('o2_therm_linux.exe',shell=True)

		#Watts#
		r=parse_tape6('TAPE6.OUT')
		WattsT.append(r['table_'+str(9)]['summary']['actinides'])
		WattsT.append(r['table_'+str(9)]['summary']['fission_products'])
		WattsT.append(r['table_'+str(9)]['summary']['activation_products'])

		for iv in range(len(WattsT)):
		    if iso2 in WattsT[iv].keys():
			totWT = WattsT[iv][iso2][1]
			    
		mt = (time-(t/2))*(sqrt(totWT * Prev_tot))
		n = n + mt
		dict[iso2] = n

	    Prev_tot= totWT
	    t=t*2
	    isos=[]
	    WattsT=[]
	    totWT=0
	    n=0
     b={}
     c = dict
     d = dict
     n = 0
     x = 0
     writer5 = open('max.txt', 'a')
     for it in c:
	 x = x+c[it]
     writer1.write(str(time)+' '+str(x)+'\n')
     while n < num:
	 em={}
	 for itt in d.keys():
	     em[d[itt]]=itt
	 b[em[max(em)]] = max(em)
	 f=[em[max(em)]]
	 isos2 = set(d.keys())
	 isos2 = isos2.difference(f)
	 d = {}
	 for iiso in isos2:
	     d[iiso] = c[iiso]
	 writer1.write(str(em[max(em)])+'_'+str(max(em))+'\n')
	 n=n+1
	 isos2=[]
     writer1.close()
     return b	