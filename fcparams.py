from bright.reactor_parameters import lwr_defaults
LWR_Params = lwr_defaults()

from bright.reactor_parameters import fr_defaults
FR_Params = fr_defaults()

#General Specifications
Quiet = True

#LWR Specifications
LWR_Params.BUt =              51.0 
LWR_Params.batches   =        3 
LWR_Params.pnl      =         0.98 

LWR_Fuel2Mod  =  0.30146786330
                            
#LWR Storage                
LWR_SNF_Storage_Time        =  6
                            
#LWR Reprocessing                         
LWR_SE_U      =              0.999 
LWR_SE_NP     =              0.99 
LWR_SE_PU     =              0.99 
LWR_SE_AM     =              0.99 
LWR_SE_CM     =              0.99 
LWR_SE_CS     =              0.0 
LWR_SE_SR     =              0.0
                            
#FR Specifications                        
#FR_Params.BUt = 176.0						#FR Burnup
FR_Params.BUt = 180.0						#FR Burnup
FR_Params.batches = 3						#Number of FR batches 
FR_Params.pnl =  0.65		 				#FR Non-Leakage Probability

FR_TRU_CR     =  0.500000
FR_LAN_FF_Cap =  0.0 * (10**-6)
                            
#FR Storage                 
FR_SNF_Storage_Time  = 3
                            
#FR Reprocessing                          
FR_SE_U       =              0.999 
FR_SE_NP      =              0.99 
FR_SE_PU      =              0.99 
FR_SE_AM      =              0.99 
FR_SE_CM      =              0.99 
FR_SE_CS      =              0.0 
FR_SE_SR      =              0.0 


#Interim Storage
INT_SNF_Storage_Time        =  30
