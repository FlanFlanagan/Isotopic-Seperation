from metasci.nuke import origen as msno
import bright
import mass_stream as MS
import subprocess
import matplotlib.pyplot as plt
import pickle
import numpy

def parsers(inputfile, time, tables):
     ## Setting up the Mass Stream ##
     instream=inputfile
     a = MS.MassStream(instream, 1)
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
    

#from tables import *
#h5file = openFile("data.txt", mode = "w", title = "Fuel Metrics")
#group = h5file.createGroup("/", 'Metrics', 'Metrics Output')
#table = h5file.creatTable(group, 'Metrics', Metrics, "Fuel Metrics")
#row = table.row
#row['name'] = 'IsoRemoved'
#row['Heat'] = totWT
#row['Toxicity'] = totTI
#row['Rad'] = TotR


