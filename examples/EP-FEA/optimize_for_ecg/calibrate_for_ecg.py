#Author: Oguz  Ziya Tikenogullari
################################################################################

import numpy as np
import csv
import os
import math as m
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import subprocess


def getrmse(x):

    global outfile
    global count
    count = count +1 

    #print ("running the analysis with %f, %f" %(x[0],x[1]) )
    Process = subprocess.call(["/Users/oguz/build-PERIGEE/ep_main/analysis_3d_EP",
                               "-myo_cond_scaler", str(x[0]),
                               "-pur_cond_scaler", str(x[1]),
                               "-LV_pur_delay",    str(x[2]), 
                               "-RV_pur_delay",    str(x[3]) ] )
                              #stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    #print ("running compute ecg" )    
    Process = subprocess.call(["/Users/oguz/build-PERIGEE/ep_main/compute_ECG"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    #print ("finished running with %f, %f, %f, %f" %(x[0],x[1],x[2],x[3]) )

    PS_time_org =np.genfromtxt('/Users/oguz/build-PERIGEE/ep_main/ecg_digitized.csv'
                               ,delimiter = ',',usecols=0)
    
    PS_ecg_org  =np.genfromtxt('/Users/oguz/build-PERIGEE/ep_main/ecg_digitized.csv'
                               ,delimiter = ',',usecols=1)
    
    Time        =np.genfromtxt('/Users/oguz/build-PERIGEE/ep_main/ecg_recording.csv'
                  ,skip_header=1,delimiter = ',',usecols=0)
    
    Sim_ecg     =np.genfromtxt('/Users/oguz/build-PERIGEE/ep_main/ecg_recording.csv'
                  ,skip_header=1,delimiter = ',',usecols=1)
    
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

x0 = np.array([0.8, 2.0, 0.0, 0.0])

dx = np.array([0.1, 0.1, 1.0, 1.0])
bnds = ((0.5, 1.5), ( 0.1, 3.0), (0.0, 100.0), ( 0.0, 100.0))

print ("initial error: ", getrmse(x0))

count=0

#res = minimize(getrmse, x0, args=(), method='SLSQP', jac=None, bounds=bnds,
#               options={'eps': 0.02,'disp': True})

res = minimize(getrmse, x0, method='nelder-mead',
               options={'xtol': 1e-8, 'disp': True})
#options={'disp': False, 'minfev': 0, 'scale': None, 'rescale': -1, 'offset': None, 'gtol': -1, 'eps': 1e-08, 'eta': -1, 'maxiter': None, 'maxCGit': -1, 'mesg_num': None, 'ftol': -1, 'xtol': -1, 'stepmx': 0, 'accuracy': 0}


print("res.x" , res.x)
print ("final error: ", getrmse(res.x))
