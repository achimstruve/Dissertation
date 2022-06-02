from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft


fn_fst = "Dlc_13T_+00deg_"
wind_speeds_names=np.array(['V04','V06','V08','V10','V11.4','V12','V14','V16','V18','V20','V22','V24'])
wind_speeds_names=np.array(['V04','V06','V08','V10','V12','V14'])
SD_tower_names = np.array(['0.0','1.0','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0','10.0'])
SD_tower_names = np.array(['7'])
file_number = len(wind_speeds_names) * len(SD_tower_names)

T_start = 0
T_add_zero_pad = 600 * 6

## FUNKTIONEN ##

def getnames(datanames):
    """This function returns the names of all data columns"""
    # count the output variables of file
    z=0
    for j in np.arange(len(datanames)-1):
        if '\t' in (datanames[j]):
            z=z+1
    names=["" for x in range(z+1)]
    # save the output variable names in "names" list
    z=0
    begin=0 # set begin 
    end=len(datanames)
    for j in np.arange(len(datanames)-1):
        if '\t' in (datanames[j]):
            end=j-1
            names[z]=datanames[begin:end]
            begin=end
            names[z]=names[z].replace(" ","")
            names[z]=names[z].replace("\t","")
            names[z]=names[z].replace("\n","")
            z=z+1
    # save the last name element in names
    names[len(names)-1]=datanames[end:len(datanames)-1]
    names[len(names)-1]=names[len(names)-1].replace(" ","")
    names[len(names)-1]=names[len(names)-1].replace("\t","")
    names[len(names)-1]=names[len(names)-1].replace("\n","")
    return names

def getdata(data_raw,names):
    # initialize data array
    data=np.zeros((len(data_raw),len(names)),dtype=float)
    # speichere data_raw - Daten in entsprechende Arrayelemente von data
    characters=int(0)
    for i in range(len(data_raw)):
        for k in range(len(names)):
            data[i][k]=float(data_raw[i][characters:characters+11])
            characters=characters+11
        characters=0
    return data

#print 'Schreibe vollstaendigen Dateinamen zum Plotten von Zeitreihen'
color_a = np.array(['b', 'g', 'r'])
#f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
f, ax1 = plt.subplots(1)
filename = "Dlc_63_+30deg.out"    
# load data from file
f=open(filename,'r') # open file
datanames=list(f)[6] # laod the list of output variables
f=open(filename,'r') # open file
data_raw=list(f)[8:1000] # load the data list of output variables

# call function getnames to get names of the data columns of the file
names=getnames(datanames)

# call function getdata to get data of the data columns of the file
data=getdata(data_raw,names)

column1 = 0
column2 = 8
# plot der entsprechenden Datenspalten
dt = np.transpose(data)[column1][1] - np.transpose(data)[column1][0]
start_index = int(T_start / dt)
N = len(np.transpose(data)[column1][start_index:]) + int(T_add_zero_pad/dt)
T = np.transpose(data)[column1][1] - np.transpose(data)[column1][0]
fft_data = np.transpose(data)[column2][start_index:] - np.mean(np.transpose(data)[column2][start_index:])
fft_zp_data = np.zeros([N])
fft_zp_data[0:len(fft_data)] = np.copy(fft_data)
fft_1 = fft(fft_zp_data)
xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
indexes_xfg04 = np.min(np.where(xf > 0.4))
fft_y = 2.0/N * np.abs(fft_1[0:N//2])**2
ax1.plot(xf, fft_y, color = color_a[0])
#ax1.set_title("$\mathrm{\overline{v}_{w}}$: " +  ws_n[1:] + " $\mathrm{m/s,\ \ \ \overline{f}_{3p}}$: " + str(round(rot_spd_mean/60 * 3,4)) + " $\mathrm{Hz}$")

params = {'legend.fontsize': 24}
plt.rcParams.update(params)
font = {'family' : 'normal',
        'size'   : 24}
plt.rc('font', **font)

ax1.legend(loc = 0)
ax1.set_xlabel('frequency / Hz')
#ax1.set_xlim([0,0.8])

ax1.grid('on')
 