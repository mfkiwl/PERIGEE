import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.patches as patches 


PS_time_org = np.genfromtxt('./ecg_digitized.csv',delimiter = ',',usecols=0)
PS_ecg_org  = np.genfromtxt('./ecg_digitized.csv',delimiter = ',',usecols=1)

Time        = np.genfromtxt('./ecg_recording.csv',skip_header=1,delimiter =
                            ',',usecols=0) 
Sim_ecg     = np.genfromtxt('./ecg_recording.csv',skip_header=1,delimiter =
                            ',',usecols=1) 

#PS_ecg      = np.interp(Time, PS_time_org, PS_ecg_org)
    
plt.plot(PS_time_org, PS_ecg_org, 'r', Time, Sim_ecg, 'b')
#plt.savefig("ECG.png")
plt.show()
