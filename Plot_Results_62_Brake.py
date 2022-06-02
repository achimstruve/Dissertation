# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:02:48 2015

@author: onrok
"""

#import sys
#sys.modules[__name__].__dict__.clear()

import numpy as np
import matplotlib.pyplot as plt

## FUNKTIONEN ##

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

answer='j'
answer2='n'
#print 'Schreibe vollstaendigen Dateinamen zum Plotten von Zeitreihen'
#filename=raw_input()
filename='Out_Files/Dlc_13T_+00deg_V20_S6.out'
#while answer=='j':
for ii in range(3):

    if ii == 1:
        filename='Out_Files/Dlc_62T_-20deg_S1.out'
    if ii == 0:
        filename='Dlc_62T_-20deg_S1_idle.out'
    if ii == 2:
        filename='Dlc_62T_-20deg_S1.out'
        #filename='Out_Files/Dlc_62T_+20deg_S3.out'
    column2 = 76
    column3 = 76
    if answer2 <> 'j':
        # load data from file
        f=open(filename,'r') # open file
        names=list(f)[6].split('\t') # laod the list of output variables
        names[-1] = names[-1][0:-1]
        f=open(filename,'r') # open file
        units=list(f)[7].split('\t') # laod the list of output variables
        units[-1] = units[-1][0:-1]
        f=open(filename,'r') # open file
        data_raw=list(f)[8:] # load the data list of output variables

    
        # call function getdata to get data of the data columns of the file
        data = np.zeros([len(data_raw), len(names)], dtype = float)
        for di in range(len(data_raw)):
            data_tmp = data_raw[di]
            data[di,:] = data_tmp.split('\t')
        
        units[0] = 's'
        for ui in range(len(units)):
            u_tmp = units[ui]
            if (len(u_tmp)>=4) and (u_tmp[3] == '\xb7'):
                u_tmp = u_tmp[1:3] + u_tmp[4]
            while u_tmp[-1] == ' ':
                u_tmp = u_tmp[0:-1]
            if u_tmp[0] == '(':
                u_tmp = u_tmp[1:]
            if u_tmp[-1] == ')':
                u_tmp = u_tmp[0:-1]
            if u_tmp[-1] == 'c':
                u_tmp = u_tmp[0:-2]
            if u_tmp[-1] == '2':
                u_tmp = u_tmp[0:-4] + '^2'
            units[ui] = u_tmp
    
    """   
    if ii==0:
        for n in range(len(names)):
            print str(n)+': '+names[n]
    """
    """
    print 'Waehle Datenreihe (Zahl) für X-Achse.'
    column1=int(raw_input()) # Datenspalte 1, welche entlang der x-Achse geplottet werden soll
    
    column1=0
    for n in range(len(names)):
        print str(n)+': '+names[n]
    print 'Waehle Datenreihe (Zahl) für 1. Y-Achse.'
    column2=int(raw_input()) # Datenspalte 2, welche entlang der 1. y-Achse geplottet werden soll
    for n in range(len(names)):
        print str(n)+': '+names[n]
    print 'Waehle Datenreihe (Zahl) für 2. Y-Achse.'
    column3=int(raw_input()) # Datenspalte 3, welche entlang der 2. y-Achse geplottet werden soll
    """
    column1 = 0
    # plot der entsprechenden Datenspalten
    winddir=(180/np.pi)*np.arctan((-np.transpose(data)[2][0:])/np.transpose(data)[1][0:]) # calculate winddirection
    show_winddir = False
    
    params = {'legend.fontsize': 20}
    plt.rcParams.update(params)
    font = {'family' : 'normal',
            'size'   : 28}
    plt.rc('font', **font)
    if ii == 0:
        fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = True)
        lns1 = ax1.plot(np.transpose(data)[column1][0:],np.transpose(data)[column2][0:]/1e3, 'k-' ,linewidth = 2, label='idling / brake not applied')
        # Make the y-axis label and tick labels match the line color.
        ax1.set_ylabel(r'$\mathrm{'+names[column2] + '}$ / $\mathrm{MNm}$')
        if show_winddir:
            ax1.plot(np.transpose(data)[column1][0:],winddir, 'b-')
            ax1.set_ylabel('WindDir', color='b')
        ax1.grid('on')
        print 'Mean of ' + names[column2] + ' : ' + str(round(np.mean(np.transpose(data)[column2][0:]),1)) + ' ' + units[column2]
        print 'standard deviation of sensor: ', np.std(np.transpose(data)[column2][0:])
    if ii == 1:
        lns2 = ax2.plot(np.transpose(data)[column1][0:],np.transpose(data)[column3][0:]/1e3, 'k-',linewidth = 2, label='idling / brake applied')
        ax2.set_ylabel(r'$\mathrm{'+names[column2] + '}$ / $\mathrm{MNm}$')
        print 'Mean of ' + names[column3] + ' : ' + str(round(np.mean(np.transpose(data)[column3][0:]),1)) + ' ' + units[column3]
        print 'standard deviation of sensor: ', np.std(np.transpose(data)[column3][0:])
    if ii == 2:
        lns3 = ax3.plot(np.transpose(data)[column1][0:],np.transpose(data)[column3][0:]/1e3, 'k-',linewidth = 2, label='parked')
        print 'Mean of ' + names[column3] + ' : ' + str(round(np.mean(np.transpose(data)[column3][0:]),1)) + ' ' + units[column3]
        print 'standard deviation of sensor: ', np.std(np.transpose(data)[column3][0:])
        ax3.set_ylabel(r'$\mathrm{'+names[column2] + '}$ / $\mathrm{MNm}$')
        ax3.set_xlabel(r'$\mathrm{'+names[column1] + '}$ / $\mathrm{' + units[column1] +'}$')
        
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.set_ylim([-50,250])
    ax2.set_ylim([-50,250])
    ax3.set_ylim([-50,250])
    ax1.grid('on')
    ax2.grid('on')
    ax3.grid('on')
    #plt.title(filename)
    

    """
    # locus curve
    plt.plot(np.transpose(data)[column2][0:], np.transpose(data)[column3][0:], 'b-')
    
    plt.grid('on')
    print 'Weiterer Plott? (j/n)'
    answer=raw_input()
    plt.show()
    if answer=='j':
        print 'Gleiche Datei? (j/n)'
        answer2=raw_input()
        if answer2 == 'n':
            print 'Schreibe vollstaendigen Dateinamen zum Plotten von Zeitreihen'
            filename=raw_input()  
    """