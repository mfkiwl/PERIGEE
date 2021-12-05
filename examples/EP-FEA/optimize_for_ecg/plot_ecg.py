import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.patches as patches 


PS_time_org = np.genfromtxt('./ecg_PS.csv',skip_header=1, delimiter =
                            ',',usecols=0) 
PSv2_ecg_org  = np.genfromtxt('./ecg_PS.csv',skip_header=1, delimiter =
                              ',',usecols=1)
PSv6_ecg_org  = np.genfromtxt('./ecg_PS.csv',skip_header=1, delimiter =
                              ',',usecols=2) 

PSv2_ecg_org  = PSv2_ecg_org/(np.amax(PSv2_ecg_org)-np.amin(PSv2_ecg_org))
PSv6_ecg_org  = PSv6_ecg_org/(np.amax(PSv6_ecg_org)-np.amin(PSv6_ecg_org))

#bring the rest state=0
PSv2_ecg_org  = PSv2_ecg_org-((PSv2_ecg_org[0]+PSv2_ecg_org[-1])/2.0) 
PSv6_ecg_org  = PSv6_ecg_org-((PSv6_ecg_org[0]+PSv6_ecg_org[-1])/2.0) 

    
Time        = np.genfromtxt('./ecg_recording.csv',skip_header=1,delimiter =
                            ',',usecols=0) 
Simv2_ecg   = np.genfromtxt('./ecg_recording.csv',skip_header=1,delimiter =
                            ',',usecols=1)
Simv6_ecg   = np.genfromtxt('./ecg_recording.csv',skip_header=1,delimiter =
                            ',',usecols=2)

Simv2_ecg     = Simv2_ecg/(np.amax(Simv2_ecg)-np.amin(Simv2_ecg))
Simv6_ecg     = Simv6_ecg/(np.amax(Simv6_ecg)-np.amin(Simv6_ecg))

#PS_ecg      = np.interp(Time, PS_time_org, PS_ecg_org)
    
plt.plot(PS_time_org, PSv2_ecg_org, 'r-.', Time, Simv2_ecg, 'r',
         PS_time_org, PSv6_ecg_org, 'k-.', Time, Simv6_ecg, 'k')
plt.savefig("ECG.png")
plt.show()
