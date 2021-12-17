import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.patches as patches

PSv2_ecg_path = "SU36-V2-1stbeat.csv"
PSv6_ecg_path = "SU36-V6-1stbeat.csv"
Sim_ecg_path = "ecg_recording.csv"


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
plt.savefig("ecg.png")
plt.show()
