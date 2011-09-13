#! /usr/bin/env python
from __future__ import print_function

import time
import math
from rmIsos import *
import tables
from AddFuncts import *
import metasci
from fcparams import *
from bright import *
from mass_stream import *
from isoname import *
from pyne.data import *
##################
### Prep Work! ###
##################

bright_config.write_hdf5 = True
bright_config.write_text = False


#Various Variables
snf_need = []
if not Quiet:
    verbosity(100)

#General Functions
def MakeSep(s):
    "Makes a dictionary of separation efficiencies based on values in fcparams."
    "s is the reactor string."
    prefix = s + "_SE_"
    prelen = len(prefix)
    seps = {}
    for g in globals().keys():
	if g.startswith(prefix):
	    seps[g[prelen:]] = globals()[g]
    return seps

#redefine iUsotrak
trackfile = tables.openFile("FR.h5", 'r')
itrack = trackfile.root.ToIso_zz.read()
trackfile.close()
bright_config.track_isos = set([int(i) for i in itrack])

#Converts storage times to seconds
for key in vars().keys():
    if key[-17:] == "_SNF_Storage_Time":
        val = float( vars()[key] )
	vars()[key] = metasci.time2sec(val, 'y')

#Calculates LWR Pin Size based on Fuel-to-Moderator Ratio
if 'LWR_Fuel2Mod' in vars().keys():
    LWR_Params.fuel_radius = LWR_Params.unit_cell_pitch * math.sqrt(LWR_Fuel2Mod / math.pi)

#Seperation Dictionaries
sepeffLWR = {"92": 0.9, "93": 0.9, "94": 1, "95": 0, "96": 0, "55": 0, "38": 0}
sepeffFR  = {"92": 0.9, "93": 0.9, "94": 1, "95": 0, "96": 0, "55": 0, "38": 0}


#Fuel Cycle Components
LWR      = LightWaterReactor1G(libfile="LWR.h5",reactor_parameters=lwr_defaults(),name= "LWR")
FR       = FastReactor1G(libfile= "FR.h5", reactor_parameters=FR_Params, name= "FR")
LWR_Rep  = Reprocess(sepeff=sepeffLWR)
FR_Rep   = Reprocess(sepeff=sepeffFR)
LWR_Stor = Storage(name = "LWR_Storage")
FR_Stor  = Storage(name = "FR_Storage")
INT_Stor = Storage(name = "INT_Storage")
Rmv_Stor = Storage(name = "Rmv_Storage")

#######################
### LWR Computation ###
#######################

def LWR_delR_BU_(ms):
    "Calculates the delta Reaction Rates at the target burnup."
    LWR.ms_feed = ms
    LWR.fold_mass_weights()
    dR = LWR.batch_average(LWR_Params.BUt, "p") - LWR.batch_average(LWR_Params.BUt,"D")
    return dR

U235 = MassStream({922350: 1.0}, 0.04, "U235")
U238 = MassStream({922380: 1.0}, 0.96, "U238")

delR_U235 = LWR_delR_BU_(U235)
delR_U238 = LWR_delR_BU_(U238)

#Calculate delta R for the Guess
LWR_CoreInput = U238 + U235
LWR_CoreInput.name = "LWR_CoreInput"
LWR_CoreInput.normalize()
LWR_delR_Guess = LWR_delR_BU_(LWR_CoreInput)


k = LWR.batch_average_k(LWR_Params.BUt)
n = 0
if not Quiet:
    print("{0}) {1}".format(1, k), end=" ") 

while 0.001 < abs(1.0 - k) and n < 10:
    #Adjust Masses based on pertubation guess.
    LWR_DeltaM_U238 = - LWR_delR_Guess / (delR_U238 - delR_U235)
    U238.mass = U238.mass + LWR_DeltaM_U238
    U235.mass = U235.mass - LWR_DeltaM_U238

    #Recalculate core parameters for new masses guess
    LWR_CoreInput = U238 + U235
    LWR_CoreInput.name ="LWR_CoreInput"
    LWR_CoreInput.normalize()
    LWR_delR_Guess = LWR_delR_BU_(LWR_CoreInput)
    k = LWR.batch_average_k(LWR_Params.BUt)
    n = n+1
    if not Quiet:
	print(k, end=" ") 
if not Quiet:
    print()
    print()

#Calculate and write output
LWR.BUd_bisection_method()
LWR.calc_ms_prod()
LWR.write()
LWR_SNF = 1.0 * LWR.ms_prod
LWR_SNF.name = "LWR_SNF"

LWR_Cooled = LWR_Stor.calc(LWR_SNF, LWR_SNF_Storage_Time)
LWR_Cooled.name = "LWR_Cooled"

RmvIsos = [RmIsos]
RmStor = MassStream({RmIsos: LWR_Cooled.comp[RmIsos]})
stor_t = 0 * half_life(RmIsos)
Rm_Stor = Rmv_Stor.calc(RmStor, stor_t)
LWR_Cooled = remove(LWR_Cooled, RmvIsos)

LWR_stream = LWR_Cooled.mult_by_mass()
with open('LWR_CooledIsos.txt', 'w') as f:
    for iso in LWR_stream.keys():
	f.write("{0:10}{1:.5E}\n".format(isoname.zzaaam_2_LLAAAM(iso), LWR_stream[iso]))
	
	
LWR_RepOut = LWR_Rep.calc(LWR_Cooled)
LWR_RepOut.name = "LWR_Reprocessing_Product"

#LWR_Rep.write_ms_pass()
#LWR_Stor.write_ms_pass()


######################
### FR Computation ###
######################

#Mass Streams
UTopUp       = MassStream("U-TopUp.txt", 0.50, "UTopUp")
TRUTopUp     = LWR_RepOut.get_tru("TRUTopUp")
FR_RepUout   = MassStream({922350: 1.0}, 0.0, "RepUout")
FR_RepTRUout = MassStream({942380: 1.0}, 0.0, "RepTRUout")
FR_RepLANout = MassStream({591440: 1.0}, 0.0, "RepLANout")

TRU_per_kgLWR_FF = TRUTopUp.mass
TRUTopUp.mass = 0.5

def FR_delR_BU_(ms):
    "Calculates the delta Reaction Rates at the target burnup."
    FR.ms_feed = ms
    FR.fold_mass_weights()
    dR = FR.batch_average(FR_Params.BUt, "p") - FR.batch_average(FR_Params.BUt, "D")
    return dR

def FR_Mass_Ratio_Calc():
    delR_UTop   = FR_delR_BU_(UTopUp)
    delR_TRUTop = FR_delR_BU_(TRUTopUp)

    #First Guess for UTopUp and TRUTopUp masses; each get half of the remaining mass space.
    TopUpMassSpace = 1.0 - FR_RepUout.mass - FR_RepTRUout.mass - FR_RepLANout.mass

    #Find bound for All U
    UTopUp.mass   = TopUpMassSpace * 1.0
    TRUTopUp.mass = TopUpMassSpace * 0.0
    CoreInput = UTopUp + TRUTopUp + FR_RepUout + FR_RepTRUout + FR_RepLANout
    CoreInput.name = "CoreInput"
    CoreInput.normalize()
    delR_Guess = FR_delR_BU_(CoreInput)
    k_AllU = FR.batch_average_k(FR_Params.BUt)
    sign_U = (1.0 - k_AllU) / abs(1.0 - k_AllU)

    #Find bound for All TRU
    UTopUp.mass   = TopUpMassSpace * 0.0
    TRUTopUp.mass = TopUpMassSpace * 1.0
    CoreInput = UTopUp + TRUTopUp + FR_RepUout + FR_RepTRUout + FR_RepLANout
    CoreInput.name = "CoreInput"
    CoreInput.normalize()
    delR_Guess = FR_delR_BU_(CoreInput)
    k_AllTRU = FR.batch_average_k(FR_Params.BUt)
    sign_TRU = (1.0 - k_AllTRU) / abs(1.0 - k_AllTRU)


    if sign_U == sign_TRU:
	raise RuntimeError("BadFuelForm: Multiplication Factor Opperates on Range {0}".format([k_AllU, k_AllTRU]))
    else:        
	#Continue nomrally
	UTopUp.mass   = TopUpMassSpace * 0.5
	TRUTopUp.mass = TopUpMassSpace - UTopUp.mass 

	#Calculate delta R for the Guess
	CoreInput = UTopUp + TRUTopUp + FR_RepUout + FR_RepTRUout + FR_RepLANout
	CoreInput.name = "CoreInput"
	CoreInput.normalize()
	delR_Guess = FR_delR_BU_(CoreInput)

	k = FR.batch_average_k(FR_Params.BUt)
	n = 0
	if not Quiet:
	    print("{0}) {1}".format(cyc+1, k), end=" ") 

	while 0.001 < abs(1.0 - k) and n < 10:
	    #Adjust Masses based on pertubation guess.
	    DeltaM_U = - delR_Guess / (delR_UTop - delR_TRUTop)
	    UTopUp.mass   = UTopUp.mass   + DeltaM_U
	    TRUTopUp.mass = TRUTopUp.mass - DeltaM_U

	    #Recalculate core parameters for new masses guess
	    CoreInput = UTopUp + TRUTopUp + FR_RepUout + FR_RepTRUout + FR_RepLANout
	    CoreInput.name = "CoreInput"
	    CoreInput.normalize()
	    delR_Guess = FR_delR_BU_(CoreInput)
	    k = FR.batch_average_k(FR_Params.BUt)
	    n = n+1
	    if not Quiet:
		print(k, end=" ") 
	if not Quiet:
	    print()
	    print()

    #Calculate and write output
    FR.BUd_bisection_method()
    FR.calc_ms_prod()
    FR.calcSubStreams()
    FR.calc_tru_cr()

    return


def FR_Calibrate_PNL_2_TRUCR():
#      delta = 0.1
#      delta = 0.05
    delta = 0.01

    #Determine Lower Bound
#   pnl_a  = 0.1
    pnl_a  = 0.30
    FoundA = False
    while not FoundA:
	try: 
	    FR.P_NL = pnl_a
	    FR_Mass_Ratio_Calc()
	    trucr_a = FR.tru_cr
	    sign_a = (trucr_a - FR_TRU_CR) / abs(trucr_a - FR_TRU_CR)
	    FoundA = True
	except RuntimeError as e:
	    if ("BadFuelForm" not in str(e)) and ("FUEL COMPOSITION NOT COMPUTABLE!" not in str(e)): 
		raise e
	    pnl_a = pnl_a + delta

    print("Yay!")

    #Determine Upper Bound
    #pnl_b  = 1.2
    pnl_b  = 0.8
    FoundB = False
    while not FoundB:
	try: 
	    FR.P_NL = pnl_b
	    FR_Mass_Ratio_Calc()
	    trucr_b = FR.tru_cr
	    sign_b = (trucr_b - FR_TRU_CR) / abs(trucr_b - FR_TRU_CR)
	    FoundB = True
	except RuntimeError as e:
	    if ("BadFuelForm" not in str(e)) and ("FUEL COMPOSITION NOT COMPUTABLE!" not in str(e)): 
		raise e
	    pnl_b = pnl_b - delta

    print("Yay!")

    DoA = 10.0**(-5)        #Degree of accuracy to carry out calculations to.
    q = 0
    while (DoA < abs(pnl_a - pnl_b)) and (DoA < abs(trucr_a - trucr_b)) and q < 30:
	#WARNING! This next block is a quick hack that sometimes fails.
	GoodBoundary = False
	n = 1
	while not GoodBoundary:
	    try:
		if not Quiet:
		    print("P_NL_c calculation try number {0}.".format(n))
		pnl_c = (pnl_a + pnl_b) / 2.0
		FR.P_NL = pnl_c
		FR_Mass_Ratio_Calc()
		GoodBoundary = True
	    except RuntimeError as e:
		if ("BadFuelForm" not in str(e)) and ("FUEL COMPOSITION NOT COMPUTABLE!" not in str(e)): 
		    raise e
		pnl_a = pnl_a + 0.1*pnl_a
		pnl_b = pnl_b - 0.1*pnl_b
		n = n + 1

	trucr_c = FR.tru_cr
	sign_c = (trucr_c - FR_TRU_CR) / abs(trucr_c - FR_TRU_CR)
	q = q + 1

	if (sign_a == sign_c) and not (sign_b == sign_c):
	    pnl_a   = pnl_c
	    trucr_a = trucr_c
	    sign_a  = sign_c
	elif (sign_b == sign_c) and not (sign_a == sign_c):
	    pnl_b   = pnl_c
	    trucr_b = trucr_c
	    sign_b  = sign_c
	else:
	    if not Quiet:
		print() 
		print("SOMETHING WENT WRONG WHILE FINDING THE TRU CONVERSION RATIO!!!")
		print("Here is some information that might help you debug ^_^")
		print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'a', 'pnl': pnl_a, 'trucr': trucr_a, 'sign': sign_a})
		print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'b', 'pnl': pnl_b, 'trucr': trucr_b, 'sign': sign_c})
		print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'c', 'pnl': pnl_c, 'trucr': trucr_c, 'sign': sign_c})
		print()

    if not Quiet:
	print()
	print("Final Result P_NL Calibration to TRU_CR via Bisection Method Calculation:")
	print("q = %i"%q)
	print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'a', 'pnl': pnl_a, 'trucr': trucr_a, 'sign': sign_a})
	print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'b', 'pnl': pnl_b, 'trucr': trucr_b, 'sign': sign_c})
	print("pnl_%(ltr)s = %(pnl).16f\ttrucr_%(ltr)s = %(trucr)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'c', 'pnl': pnl_c, 'trucr': trucr_c, 'sign': sign_c})
	print()

for cyc in range(10):
    if cyc in [0]:
	FR_Calibrate_PNL_2_TRUCR()
    else:
	delta = 0.001
	BadKRange = True
	while BadKRange:
	    try:
		FR_Mass_Ratio_Calc()
		BadKRange = False
	    except RuntimeError as e:
		if "FUEL COMPOSITION NOT COMPUTABLE!" not in str(e): 
		    raise e

		pnl_regime = float(str(e).split()[-1][:-1])
		if pnl_regime < 1.0:
		    FR.P_NL = FR.P_NL + delta
		elif 1.0  <= pnl_regime:
		    FR.P_NL = FR.P_NL - delta
    FR.write()

    #Calculate the LWR SNF Top up needed
    snf_need.append( TRUTopUp.mass / TRU_per_kgLWR_FF )

    StorOut = FR_Stor.calc(FR.ms_prod, FR_SNF_Storage_Time)
    StorOut.name = "StorOut"
    FR_Stor.write()

    FR_RepOut = FR_Rep.calc(StorOut)
    FR_RepOut.name = "RepOut"
    FR_Rep.write()

    FR_RepUout   = FR_RepOut.get_u(FR_RepUout.name)
    FR_RepTRUout = FR_RepOut.get_tru(FR_RepTRUout.name)
    FR_RepLANout = FR_RepOut.get_lan(FR_RepLANout.name)
    if FR_LAN_FF_Cap < FR_RepLANout.mass:
	FR_RepLANout.mass = FR_LAN_FF_Cap

#Write the SNF Needed line to output file
params = open(FR.name + "Params.txt", 'a')
params.write("LWR_SNF\t")
for el in snf_need:
    params.write( "%.6E\t%.6E\t"%(el, 0.0) )
params.write("\n")
params.close()


FR_stream = StorOut.mult_by_mass()
with open('FR_CooledIsos.txt', 'w') as f:
    for iso in FR_stream.keys():
	f.write("{0:10}{1:.5E}\n".format(isoname.zzaaam_2_LLAAAM(iso), FR_stream[iso]))


#################################
### Construct HLW Mass Stream ###
#################################
#Define other FP
other_FP = []
for i in FP:
    if i in LAN:
	continue
    elif i in ["CS", "SR"]:
	continue
    else:
	other_FP.append(i)

#Fist get LWR HLW
LWR_SNF_U   = LWR_Stor.ms_prod.get_u()
LWR_SNF_NP  = LWR_Stor.ms_prod.get_sub_stream(["NP"])
LWR_SNF_PU  = LWR_Stor.ms_prod.get_pu()
LWR_SNF_AM  = LWR_Stor.ms_prod.get_sub_stream(["AM"])
LWR_SNF_CM  = LWR_Stor.ms_prod.get_sub_stream(["CM"])
LWR_SNF_CS  = LWR_Stor.ms_prod.get_sub_stream(["CS"])
LWR_SNF_SR  = LWR_Stor.ms_prod.get_sub_stream(["SR"])
LWR_SNF_LAN = LWR_Stor.ms_prod.get_lan()
LWR_SNF_oFP = LWR_Stor.ms_prod.get_sub_stream(other_FP)

LWR_SNF_oFP = LWR_SNF_oFP + MassStream({10010: 1.0 - LWR_SNF_U.mass - \
    LWR_SNF_NP.mass - LWR_SNF_PU.mass - LWR_SNF_AM.mass - LWR_SNF_CM.mass - \
    LWR_SNF_CS.mass - LWR_SNF_SR.mass - LWR_SNF_LAN.mass - LWR_SNF_oFP.mass})

LWR_HLW = LWR_SNF_oFP + LWR_SNF_LAN + \
    ((1.0 - LWR_SE_U)  * LWR_SNF_U) + \
    ((1.0 - LWR_SE_NP) * LWR_SNF_NP) + \
    ((1.0 - LWR_SE_PU) * LWR_SNF_PU) + \
    ((1.0 - LWR_SE_AM) * LWR_SNF_AM) + \
    ((1.0 - LWR_SE_CM) * LWR_SNF_CM) + \
    ((1.0 - LWR_SE_CS) * LWR_SNF_CS) + \
    ((1.0 - LWR_SE_SR) * LWR_SNF_SR)

#mass of LWR_HLW = (kgLWR_HLW / kgLWR_SNF) * (kgLWR_SNF / kgFR_FF)
LWR_HLW = LWR_HLW * snf_need[-1]  


#Then get FR HLW
FR_SNF_U   = FR_Stor.ms_prod.get_u()
FR_SNF_NP  = FR_Stor.ms_prod.get_sub_stream(["NP"])
FR_SNF_PU  = FR_Stor.ms_prod.get_pu()
FR_SNF_AM  = FR_Stor.ms_prod.get_sub_stream(["AM"])
FR_SNF_CM  = FR_Stor.ms_prod.get_sub_stream(["CM"])
FR_SNF_CS  = FR_Stor.ms_prod.get_sub_stream(["CS"])
FR_SNF_SR  = FR_Stor.ms_prod.get_sub_stream(["SR"])
FR_SNF_LAN = FR_Stor.ms_prod.get_lan()
FR_SNF_oFP = FR_Stor.ms_prod.get_sub_stream(other_FP)

FR_SNF_oFP = FR_SNF_oFP + MassStream({10010: 1.0 - FR_SNF_U.mass - \
    FR_SNF_NP.mass - FR_SNF_PU.mass - FR_SNF_AM.mass - FR_SNF_CM.mass - \
    FR_SNF_CS.mass - FR_SNF_SR.mass - FR_SNF_LAN.mass - FR_SNF_oFP.mass})

FR_HLW = FR_SNF_oFP + FR_SNF_LAN + \
    ((1.0 - FR_SE_U)  * FR_SNF_U) + \
    ((1.0 - FR_SE_NP) * FR_SNF_NP) + \
    ((1.0 - FR_SE_PU) * FR_SNF_PU) + \
    ((1.0 - FR_SE_AM) * FR_SNF_AM) + \
    ((1.0 - FR_SE_CM) * FR_SNF_CM) + \
    ((1.0 - FR_SE_CS) * FR_SNF_CS) + \
    ((1.0 - FR_SE_SR) * FR_SNF_SR)

#Finally
HLW = FR_HLW + LWR_HLW + Rm_Stor
HLW.normalize()

######################################
### Do Interim storage calculation ###
######################################
HLW_Cooled = INT_Stor.calc(HLW, INT_SNF_Storage_Time)
INT_Stor.calc_params()
INT_Stor.write()

HLW_stream = HLW_Cooled.mult_by_mass()
with open('HLW_CooledIsos.txt', 'w') as f:
    for iso in HLW_stream.keys():
	f.write("{0:10}{1:.5E}\n".format(isoname.zzaaam_2_LLAAAM(iso), HLW_stream[iso]))
       
writer2 = open('BUd.py','a')
n = LWR.BUd*snf_need[-1]
m = FR.BUd
writer2.write('LWR_BUd =' + str(n)+'\n')
writer2.write('FR_BUd =' + str(m)+'\n')
writer2.close