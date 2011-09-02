from AddFuncts import *
import os
import mass_stream as MS

parse('LWR_CooledIsos.txt', 1, [7,9,13,15,17])
parse('HLW_CooledIsos.txt', 36525, [7,9,13,15,17])
parse('HLW_CooledIsos.txt', 3652500, [7,9,13,15,17])
heat('HLW_CooledIsos.txt', 36525)
heat('HLW_CooledIsos.txt', 3652500)    
    
    
