import os
import mass_stream as MS

h=os.getcwd()
a = MS.MassStream('CooledIsos.txt', 1)
for keys in a.comp.keys():    
   if keys/10000 > 88:
      try:
	os.mkdir(str(keys))
	os.chdir(str(keys))
	print(str(keys))
	writer = open('rmIsos.py','w')
	writer.write('RmIsos =' + str(keys))
	writer.close
	writer2 = open('BUd.py','w')
	writer2.close
      except OSError:
	pass
   
   os.chdir(h)