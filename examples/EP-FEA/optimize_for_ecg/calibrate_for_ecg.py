#Author: Oguz  Ziya Tikenogullari
################################################################################

import numpy as np
import csv
import os
import math as m
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import subprocess

#for sherlock 
perigee_build_dir = "/home/users/oguzziya/build-perigee/su36_calibrate/" 
PSv2_ecg_path = "/home/users/oguzziya/PERIGEE/examples/EP-FEA/optimize_for_ecg/SU36-V2-1stbeat.csv"
PSv6_ecg_path = "/home/users/oguzziya/PERIGEE/examples/EP-FEA/optimize_for_ecg/SU36-V6-1stbeat.csv"
Sim_ecg_path = "ecg_recording.csv"

##for mika
#perigee_build_dir = "/Users/oguz/build-PERIGEE/ep_main/"
#PS_ecg_path = "/Users/oguz/build-PERIGEE/ep_main/ecg_PS.csv"

def getrmse(x):

    global outfile
    global count
    count = count +1 

    #print ("running the analysis with %f, %f" %(x[0],x[1]) )
    Process = subprocess.call([#"mpirun", "-np","6",
                               "srun",
                               (perigee_build_dir+"analysis_3d_EP"),
                               "-myo_cond_scaler",   str(x[0]),
                               "-LVpur_cond_scaler", str(x[1]),
                               "-RVpur_cond_scaler", str(x[2]),
                               "-LV_pur_delay",      str(x[3]), 
                               "-RV_pur_delay",      str(x[4]) ],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.STDOUT)
                              #    stdin=subprocess.DEVNULL,
    
    #print ("running compute ecg" )    
    Process = subprocess.call([#"mpirun", "-np","6",
                               "srun",
                               (perigee_build_dir+"compute_ECG")],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.STDOUT)

    #print ("finished running with %f, %f, %f, %f" %(x[0],x[1],x[2],x[3]) )

    # get the patient specific ecg that is the aimed ecg to match and normalize it.
    PSv2_time_org = np.genfromtxt(PSv2_ecg_path, delimiter = ',', usecols=0)
    PSv6_time_org = np.genfromtxt(PSv6_ecg_path, delimiter = ',', usecols=0)
    
    #normalize amp=1,     #bring the rest state=0
    PSv2_ecg_org  = np.genfromtxt(PSv2_ecg_path, delimiter = ',', usecols=1)
    PSv2_ecg_org  = PSv2_ecg_org/(np.amax(PSv2_ecg_org)-np.amin(PSv2_ecg_org))
    PSv2_ecg_org  = PSv2_ecg_org - ((PSv2_ecg_org[0]+PSv2_ecg_org[-1])/2.0)

    PSv6_ecg_org  = np.genfromtxt(PSv6_ecg_path, delimiter = ',', usecols=1)
    PSv6_ecg_org  = PSv6_ecg_org/(np.amax(PSv6_ecg_org)-np.amin(PSv6_ecg_org)) 
    PSv6_ecg_org  = PSv6_ecg_org - ((PSv6_ecg_org[0]+PSv6_ecg_org[-1])/2.0) 
    
    # get the simulated ecg that is being modified to match PS ecg     
    Time        =np.genfromtxt(Sim_ecg_path,skip_header=1,delimiter =',',
                               usecols=0)  
    
    Simv2_ecg     =np.genfromtxt(Sim_ecg_path,skip_header=1,delimiter =',',
                                 usecols=1) 
    Simv2_ecg     = Simv2_ecg/(np.amax(Simv2_ecg)-np.amin(Simv2_ecg))

    Simv6_ecg     =np.genfromtxt(Sim_ecg_path,skip_header=1,delimiter =',',
                                 usecols=2) 
    Simv6_ecg     = Simv6_ecg/(np.amax(Simv6_ecg)-np.amin(Simv6_ecg))    

    #interpolate the PS ecg so that time points are the same as the simulated one.
    PSv2_ecg      = np.interp(Time, PSv2_time_org, PSv2_ecg_org)
    PSv6_ecg      = np.interp(Time, PSv6_time_org, PSv6_ecg_org)
    
    plt.plot(PSv2_time_org, PSv2_ecg_org, 'r-.', Time, Simv2_ecg, 'r',
             PSv6_time_org, PSv6_ecg_org, 'k-.', Time, Simv6_ecg, 'k')
    plt.savefig("ecg_%d.png"%count)
    plt.clf()
    #plt.show()
    
    rmse = ( np.sqrt(np.mean((Simv2_ecg-PSv2_ecg)**2))
            + np.sqrt(np.mean((Simv6_ecg-PSv6_ecg)**2)) )
    
    print ("iter%d, rmse with %1.2f, %1.2f, %1.2f, %1.2f, %1.2f is: "
           %(count,x[0],x[1],x[2],x[3],x[4]), rmse)
    
    f = open(outfile,"a")
    f.write("%d, \t %1.2f, %1.2f, %1.2f, %1.2f, %1.2f \t %f \n"
            %(count,x[0],x[1],x[2],x[3],x[4],rmse))
    f.close()

    return rmse

#===================================================================================

os.system("rm -rf *.png")

global Aim
outfile = "ecg_logfile.dat"

f = open(outfile,"a")
f.write("iter, \t X values, \t RMSE \n")
f.close()

count=0
#x values are : (1) myocardium conduction scaler,
#               (2) LV purkinje conduction scaler,
#               (3) RV purkinje conduction scaler, 
#               (4) LV purkinje activation delay,
#               (5) RV purkinje activation delay

x0 = np.array([1.0, 1.2, 0.3, 100.0, 140.0])

dx = np.array([0.1, 0.1, 0.1, 5.0, 5.0])
bnds = ((0.5, 1.5), ( 0.1, 3.0), ( 0.1, 3.0), (0.0, 400.0), ( 0.0, 400.0))

print ("initial error: ", getrmse(x0))

count=0

#res = minimize(getrmse, x0, args=(), method='SLSQP', jac=None, bounds=bnds,
#               options={'eps': 0.02,'disp': True})

res = minimize(getrmse, x0, method='nelder-mead',
               options={'xtol': 1e-8, 'disp': True})
#options={'disp': False, 'minfev': 0, 'scale': None, 'rescale': -1, 'offset': None, 'gtol': -1, 'eps': 1e-08, 'eta': -1, 'maxiter': None, 'maxCGit': -1, 'mesg_num': None, 'ftol': -1, 'xtol': -1, 'stepmx': 0, 'accuracy': 0}

print("res.x" , res.x)
print ("final error: ", getrmse(res.x))
