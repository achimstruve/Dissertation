# -*- coding: utf-8 -*-

#import sys
#sys.modules[__name__].__dict__.clear()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

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


root = "Out_Files/"
filename_out1 = "Dlc_61_"
#filename_out2 = "deg_S1_V" # DLC 13
filename_out2 = "deg"
filename_end = ".out"
file_spec = np.array(["03","05","07","09","11","13","15","17","19","21","23","25"])
#file_spec2 = np.array(["-10","+00","+10"]) # DLC 13
file_spec2 = np.array(["-15","-10","-05","+00","+05","+10","+15"]) # DLC 61
#file_spec2 = np.array(["-30","-25","-20","-15","-10","-05","+00","+05","+10","+15","+20","+25","+30"]) # DLC 61
filename_array = []
for j in range(len(file_spec2)):
    """
    # DLC 13
    for i in range(len(file_spec)):
        filename_array.append(root+filename_out1+file_spec2[j]+filename_out2+file_spec[i]+filename_end)
    """
    # DLC 61 or DLC 63
    filename_array.append(root+filename_out1+file_spec2[j]+filename_out2+filename_end)

# vx = 1
# Mx = 82
# My = 83
# get length of dataset
f=open(filename_array[0],'r') # open file
data_raw=list(f)[8:] # load the data list of output variables
# print all avalible sensors
f=open(filename_array[0],'r') # open file
datanames=list(f)[6] # laod the list of output variables
names=getnames(datanames)
for n in range(len(names)):
    print str(n)+': '+names[n]

# initialize all dataset maticies
data_set_Vx=np.empty([len(data_raw),len(filename_array)],dtype=list)
data_set_Mx=np.empty([len(data_raw),len(filename_array)],dtype=list)
data_set_My=np.empty([len(data_raw),len(filename_array)],dtype=list)
meanVx=np.empty([len(filename_array)],dtype=list)
meanMx=np.empty([len(filename_array)],dtype=list)
meanMy=np.empty([len(filename_array)],dtype=list)
maxMx=np.empty([len(filename_array)],dtype=list)
maxMy=np.empty([len(filename_array)],dtype=list)
minMx=np.empty([len(filename_array)],dtype=list)
minMy=np.empty([len(filename_array)],dtype=list)
Mres=np.empty([len(filename_array)],dtype=list)
z=0
for filename in filename_array:
    
    # load data from file
    f=open(filename,'r') # open file
    datanames=list(f)[6] # laod the list of output variables
    f=open(filename,'r') # open file
    data_raw=list(f)[8:] # load the data list of output variables
    
    # call function getnames to get names of the data columns of the file
    names=getnames(datanames)
    
    # call function getdata to get data of the data columns of the file
    data=getdata(data_raw,names)
    
    meanVx[z]=np.mean(np.transpose(data)[1][0:])
    meanMx[z]=np.mean(np.transpose(data)[82][0:])
    meanMy[z]=np.mean(np.transpose(data)[83][0:])
    maxMx[z]=np.max(np.transpose(data)[82][0:])
    maxMy[z]=np.max(np.transpose(data)[83][0:])
    minMx[z]=np.min(np.transpose(data)[82][0:])
    minMy[z]=np.min(np.transpose(data)[83][0:])
    Mres[z]=np.sqrt(np.max(np.abs(np.transpose(data)[82][0:]))**2 + np.max(np.abs(np.transpose(data)[83][0:]))**2)
    for i in range(len(np.transpose(data)[1][0:])):
        data_set_Vx[i,z]=np.transpose(data)[1][i]
        data_set_Mx[i,z]=np.transpose(data)[82][i]
        data_set_My[i,z]=np.transpose(data)[83][i]
    z=z+1


# sort datasets
meanMx_dataset = np.empty([len(meanVx)],dtype = list)
meanMy_dataset = np.empty([len(meanVx)],dtype = list)
maxMx_dataset = np.empty([len(meanVx)],dtype = list)
maxMy_dataset = np.empty([len(meanVx)],dtype = list)
minMx_dataset = np.empty([len(meanVx)],dtype = list)
minMy_dataset = np.empty([len(meanVx)],dtype = list)
Mres_dataset = np.empty([len(meanVx)],dtype = list)
meanVxs = np.empty([len(meanVx)],dtype = list)
meanMxs = np.empty([len(meanVx)],dtype = list)
meanMys = np.empty([len(meanVx)],dtype = list)
maxMxs = np.empty([len(meanVx)],dtype = list)
maxMys = np.empty([len(meanVx)],dtype = list)
minMxs = np.empty([len(meanVx)],dtype = list)
minMys = np.empty([len(meanVx)],dtype = list)
Mress = np.empty([len(meanVx)],dtype = list)

for i in range(len(meanVx)):
    meanMx_dataset[i] = (meanVx[i],meanMx[i])
    maxMx_dataset[i] = (meanVx[i],maxMx[i])
    minMx_dataset[i] = (meanVx[i],minMx[i])
    meanMy_dataset[i] = (meanVx[i],meanMy[i])
    maxMy_dataset[i] = (meanVx[i],maxMy[i])
    minMy_dataset[i] = (meanVx[i],minMy[i])
    Mres_dataset[i] = (meanVx[i],Mres[i])

meanMx_dataset_sorted = sorted(meanMx_dataset,key = lambda vw: vw[0])
maxMx_dataset_sorted = sorted(maxMx_dataset,key = lambda vw: vw[0])
minMx_dataset_sorted = sorted(minMx_dataset,key = lambda vw: vw[0])
meanMy_dataset_sorted = sorted(meanMy_dataset,key = lambda vw: vw[0])
maxMy_dataset_sorted = sorted(maxMy_dataset,key = lambda vw: vw[0])
minMy_dataset_sorted = sorted(minMy_dataset,key = lambda vw: vw[0])
Mres_dataset_sorted = sorted(Mres_dataset,key = lambda vw: vw[0])

for i in range(len(meanVx)):
    meanMx_dataset_sorted_tmp = meanMx_dataset_sorted[i]
    meanVxs[i] = meanMx_dataset_sorted_tmp[0]
    meanMxs[i] = meanMx_dataset_sorted_tmp[1]
    maxMx_dataset_sorted_tmp = maxMx_dataset_sorted[i]
    maxMxs[i] = maxMx_dataset_sorted_tmp[1]
    minMx_dataset_sorted_tmp = minMx_dataset_sorted[i]
    minMxs[i] = minMx_dataset_sorted_tmp[1]
    
    meanMy_dataset_sorted_tmp = meanMy_dataset_sorted[i]
    meanMys[i] = meanMy_dataset_sorted_tmp[1]
    maxMy_dataset_sorted_tmp = maxMy_dataset_sorted[i]
    maxMys[i] = maxMy_dataset_sorted_tmp[1]
    minMy_dataset_sorted_tmp = minMy_dataset_sorted[i]
    minMys[i] = minMy_dataset_sorted_tmp[1]
    Mres_dataset_sorted_tmp = Mres_dataset_sorted[i]
    Mress[i] = Mres_dataset_sorted_tmp[1]

#for i in range(len(filename_array)):
    #plt.scatter(data_set_Vx[0:,i],np.abs(data_set_My[0:,i]),color='r')
    #plt.scatter(data_set_Vx[0:,i],np.abs(data_set_Mx[0:,i]),color='b')

fig = plt.figure()
ax = plt.subplot(111)

font = {'family' : 'normal',
        'size'   : 26}
plt.rc('font', **font)
ax.grid('on')
scaling=1000
wind_direction = file_spec2.astype(np.float)
# DLC 61 or 63
ax.plot(wind_direction,Mres/scaling,color='k')
ax.plot(wind_direction,maxMx/scaling,color='b')
ax.plot(wind_direction,maxMy/scaling,color='r')
ax.plot(wind_direction,minMx/scaling,'--',color='b')
ax.plot(wind_direction,minMy/scaling,'--',color='r')
ax.plot(wind_direction,meanMx/scaling,'-.',color='b')
ax.plot(wind_direction,meanMy/scaling,'-.',color='r')
"""
# DLC 13
ax.plot(meanVx[24:36],Mres[24:36]/scaling,color='k')
ax.plot(meanVx[24:36],maxMx[24:36]/scaling,color='b')
ax.plot(meanVx[24:36],maxMy[24:36]/scaling,color='r')
ax.plot(meanVx[24:36],minMx[24:36]/scaling,'--',color='b')
ax.plot(meanVx[24:36],minMy[24:36]/scaling,'--',color='r')
ax.plot(meanVx[24:36],meanMx[24:36]/scaling,'-.',color='b')
ax.plot(meanVx[24:36],meanMy[24:36]/scaling,'-.',color='r')
"""
fontP = FontProperties()
#fontP.set_size('small')
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
ax.legend(['$M_{res,max}$','$M_{x,max}$','$M_{y,max}$','$M_{x,min}$','$M_{y,min}$','$M_{x,mean}$','$M_{y,mean}$'], prop = fontP,loc='upper center', bbox_to_anchor=(0.5, -0.09),fancybox=True, shadow=True, ncol=10)
ax.set_xlim([-15,15])
ax.set_xlabel('yaw misalignment / '+'$deg$')
ax.set_ylabel('tower base bending moments / $MNm$')
