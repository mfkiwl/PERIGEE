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
#for mika
#perigee_build_dir = "/Users/oguz/build-PERIGEE/ep_main/"

#patient specific ecg's path. this is the aimed ecg
#PS_ecg_path = "/Users/oguz/PERIGEE/examples/EP-FEA/optimize for ecg/ecg_digitized.csv"
PS_ecg_path = "/home/users/oguzziya/SU36-V2-1stbeat.csv"

def getrmse(x):

    global outfile
    global count
    count = count +1 

    #print ("running the analysis with %f, %f" %(x[0],x[1]) )
    Process = subprocess.call([#"mpirun", "-np","6",
                               "srun",
                               (perigee_build_dir+"analysis_3d_EP"),
                               "-myo_cond_scaler", str(x[0]),
                               "-pur_cond_scaler", str(x[1]),
                               "-LV_pur_delay",    str(x[2]), 
                               "-RV_pur_delay",    str(x[3]) ])
    #stdout=subprocess.DEVNULL,
    #stderr=subprocess.STDOUT)
    #stdin=subprocess.DEVNULL,
    
    #print ("running compute ecg" )    
    Process = subprocess.call([#"mpirun", "-np","6",
                               "srun",
                               (perigee_build_dir+"compute_ECG")],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.STDOUT)

    #print ("finished running with %f, %f, %f, %f" %(x[0],x[1],x[2],x[3]) )

    # get the patient specific ecg that is the aimed ecg to match and normalize it.
    PS_time_org =np.genfromtxt(PS_ecg_path, delimiter = ',', usecols=0)
    
    PS_ecg_org  =np.genfromtxt(PS_ecg_path, delimiter = ',', usecols=1)
    PS_ecg_org  = PS_ecg_org/(np.amax(PS_ecg_org)-np.amin(PS_ecg_org)) #normalize amp=1
    PS_ecg_org  = PS_ecg_org - ((PS_ecg_org[0]+PS_ecg_org[-1])/2.0) #bring the rest state=0 
    
    # get the simulated ecg that is being modified to match PS ecg     
    Time        =np.genfromtxt('ecg_recording.csv',skip_header=1,delimiter =',',
                               usecols=0)  
    
    Sim_ecg     =np.genfromtxt('ecg_recording.csv',skip_header=1,delimiter =',',
                               usecols=1) 
    Sim_ecg     = Sim_ecg/(np.amax(Sim_ecg)-np.amin(Sim_ecg))    

    #interpolate the PS ecg so that time points are the same as the simulated one.
    PS_ecg      = np.interp(Time, PS_time_org, PS_ecg_org)
    
    plt.plot(PS_time_org, PS_ecg_org, 'r', Time, Sim_ecg, 'b')
    plt.savefig("ecg_%d.png"%count)
    plt.clf()
    #plt.show()
    
    rmse = np.sqrt(np.mean((Sim_ecg-PS_ecg)**2))
    
    print ("iter%d, rmse with %1.2f, %1.2f, %1.2f, %1.2f is: " %(count,x[0],x[1],x[2],x[3]), rmse)
    
    f = open(outfile,"a")
    f.write("%d, \t %1.2f, %1.2f, %1.2f, %1.2f \t %f \n" %(count,x[0],x[1],x[2],x[3],rmse))
    f.close()

    return rmse

#===================================================================================

os.system("trash *.png")
os.system("trash *logfile*")

global Aim
outfile = "ecg_logfile.dat"

f = open(outfile,"a")
f.write("iter, \t X values, \t RMSE \n")
f.close()

count=0
#x values are : (1) myocardium conduction scaler,
#               (2) purkinje conduction scaler, 
#               (3) LV purkinje activation delay,
#               (4) RV purkinje activation delay

x0 = np.array([1.0, 1.0, 140.0, 140.0])

dx = np.array([0.1, 0.1, 5.0, 5.0])
bnds = ((0.5, 1.5), ( 0.1, 3.0), (0.0, 400.0), ( 0.0, 400.0))

print ("initial error: ", getrmse(x0))

count=0

#res = minimize(getrmse, x0, args=(), method='SLSQP', jac=None, bounds=bnds,
#               options={'eps': 0.02,'disp': True})

res = minimize(getrmse, x0, method='nelder-mead',
               options={'xtol': 1e-8, 'disp': True})
#options={'disp': False, 'minfev': 0, 'scale': None, 'rescale': -1, 'offset': None, 'gtol': -1, 'eps': 1e-08, 'eta': -1, 'maxiter': None, 'maxCGit': -1, 'mesg_num': None, 'ftol': -1, 'xtol': -1, 'stepmx': 0, 'accuracy': 0}

print("res.x" , res.x)
print ("final error: ", getrmse(res.x))
