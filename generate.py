import os
import mass_stream as MS

h=os.getcwd()
a = MS.MassStream('CooledIsos.txt', 1)
for keys in a.comp.keys():    
   try:
	os.mkdir(str(keys))
   except OSError:
	pass
   os.chdir(str(keys))
   print(str(keys))
   writer = open('rmIsos.py','w')
   writer.write('RmIsos =' + str(keys))
   writer.close
   os.chdir(h)
   