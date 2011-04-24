from metasci.nuke import origen as msno
import BriPy
import MassStream as MS
import subprocess
import matplotlib.pyplot as plt
import pickle
import numpy


## Setting up the Mass Stream ##
a = MS.MassStream ("HLW_CooledIsos.txt", 1000)
a.comp
msno.write_tape4(a.comp, name='Tape4.INP')

## Running Origen ##
msno.write_tape5_irradiation("IRP", 365.2425, 0, (204, 205, 206), name='Tape5.INP', out_table_num=(None))
subprocess.call('O2_Therm',shell=True)

#Watts#
WattsT=[]
totWT=[]
el_dictT={}
r=msno.parse_tape6(name='TAPE6.OUT')
WattsT.append(r['table_9']['summary']['data'])
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

for t in range(len(WattsT)):
    for key in WattsT[t].keys():
        el=key/10000
        if el in el_dictT.keys():
            el_dictT[el][t]=el_dictT[el][t]+WattsT[t][key]
        else:
            el_dictT[el]=[]
            for tt in range(len(WattsT)):
                el_dictT[el].append(0.0)
            el_dictT[el][t]=el_dictT[el][t]+WattsT[t][key]


l=[(max(el_dictT[el]), el) for el in el_dictT.keys()]
l.sort(reverse=True)
tlt=[totWT, el_dictT[l[0][1]], l[0][1]]
#writer =open ('heat.txt','w')
#pickle.dump(tlt, writer )     
#plt.plot(range(1,20), totWT, label='total')

#Toxicity#
ToxIng=[]
totTI=[]
el_dictTI={}
r=msno.parse_tape6(name='TAPE6.OUT')
ToxIng.append(r['table_17']['summary']['data'])
totTI.append(0.0)

isos = set(ToxIng[0].keys())
for i in range(len(ToxIng)):
    for iso in ToxIng[i].keys():
        isos.add(iso)

for iso in isos:
    isoT = []

    for i in range(len(ToxIng)):
        if iso in ToxIng[i].keys():
            isoT.append(ToxIng[i][iso])
            totTI[i] = totTI[i] + ToxIng[i][iso]
        else:
            isoT.append(0.0)
    
    #plt.plot( range(1,20), isoT, label=str(iso))

for t in range(len(ToxIng)):
    for key in ToxIng[t].keys():
        el=key/10000
        if el in el_dictTI.keys():
            el_dictTI[el][t]=el_dictTI[el][t]+ToxIng[t][key]
        else:
            el_dictTI[el]=[]	
            for tt in range(len(ToxIng)):
                el_dictTI[el].append(0.0)
            el_dictTI[el][t]=el_dictTI[el][t]+ToxIng[t][key]
l=[(max(el_dictTI[el]), el) for el in el_dictTI.keys()]
l.sort(reverse=True)
tlt=[totTI, el_dictTI[l[0][1]], l[0][1]]
#writer =open ('heat.txt','w')
#pickle.dump(tlt, writer )     
#plt.plot(range(1,20), TotTI, label='total')

#Radioactivity#
Rad=[]
TotR=[]
el_dictR={}
r=msno.parse_tape6(name='TAPE6.OUT')
Rad.append(r['table_7']['summary']['data'])
TotR.append(0.0)

isos = set(Rad[0].keys())
for i in range(len(Rad)):
    for iso in Rad[i].keys():
        isos.add(iso)

for iso in isos:
    isoT = []

    for i in range(len(Rad)):
        if iso in Rad[i].keys():
            isoT.append(Rad[i][iso])
            TotR[i] = TotR[i] + Rad[i][iso]
        else:
            isoT.append(0.0)
    
    #plt.plot( range(1,20), isoT, label=str(iso))

for t in range(len(Rad)):
    for key in Rad[t].keys():
        el=key/10000
        if el in el_dictR.keys():
            el_dictR[el][t]=el_dictR[el][t]+Rad[t][key]
        else:
            el_dictR[el]=[]
            for tt in range(len(Rad)):
                el_dictR[el].append(0.0)
            el_dictR[el][t]=el_dictR[el][t]+Rad[t][key]
l=[(max(el_dictR[el]), el) for el in el_dictR.keys()]
l.sort(reverse=True)
tlt=[TotR, el_dictR[l[0][1]], l[0][1]]
#writer =open ('heat.txt','w')
#pickle.dump(tlt, writer )     
#plt.plot(range(1,20), TotR, label='total')

#Creating the Table
writer2 = open('data.txt','w')
for metric in ['Watts'+str(totWT), 'Toxicity'+str(totTI), 'Rads'+str(TotR)]:
    print >> writer2, metric
writer2.close

#from tables import *
#h5file = openFile("data.txt", mode = "w", title = "Fuel Metrics")
#group = h5file.createGroup("/", 'Metrics', 'Metrics Output')
#table = h5file.creatTable(group, 'Metrics', Metrics, "Fuel Metrics")
#row = table.row
#row['name'] = 'IsoRemoved'
#row['Heat'] = totWT
#row['Toxicity'] = totTI
#row['Rad'] = TotR


