# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:02:48 2015

@author: onrok
"""

#import sys
#sys.modules[__name__].__dict__.clear()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

## FUNKTIONEN ##
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

answer='j'
answer2='n'
#print 'Schreibe vollstaendigen Dateinamen zum Plotten von Zeitreihen'
#filename=raw_input()
filename='Dlc_63_+30deg.out'
#while answer=='j':
for ii in range(2):
    if ii == 0:
        filename='Dlc_63_+30deg.out'
        column2 = 8
        column3 = 9
    if ii == 1:
        filename='Dlc_63_+30deg1.out'
        column2 = 8
        column3 = 9

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
            
    if ii==0:
        for n in range(len(names)):
            print str(n)+': '+names[n]
    
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
    
    params = {'legend.fontsize': 28}
    plt.rcParams.update(params)
    font = {'family' : 'normal',
            'size'   : 36}
    plt.rc('font', **font)
    C2max, C2min = peakdet(np.transpose(data)[column2][0:],.3)
    C3max, C3min = peakdet(np.transpose(data)[column3][0:],.3)
    C2maxt = np.zeros(np.shape(C2max)[0])
    C2mint = np.zeros(np.shape(C2min)[0])
    C3maxt = np.zeros(np.shape(C3max)[0])
    C3mint = np.zeros(np.shape(C3min)[0])
    for ci in range(len(C2maxt)):
        C2maxt[ci] = np.transpose(data)[column1][C2max[ci,0]]
    for ci in range(len(C2mint)):
        C2mint[ci] = np.transpose(data)[column1][C2min[ci,0]]
    for ci in range(len(C3maxt)):
        C3maxt[ci] = np.transpose(data)[column1][C3max[ci,0]]
    for ci in range(len(C3mint)):
        C3mint[ci] = np.transpose(data)[column1][C3min[ci,0]]
    if ii == 0:
        fig, ax2 = plt.subplots()
        ax2.set_xlabel(r'$\mathrm{'+names[column1] + '}$ / $\mathrm{' + units[column1] +'}$')
        ax2.plot(C3maxt,C3max[:,1], 'k-', linewidth=2, label='DLC 6.2 '+r'$V_\mathrm{e50}$, $+30$ $\mathrm{deg}$')
        ax2.plot(C3mint,C3min[:,1], 'k-', linewidth=2)
        ax2.set_ylabel(r'$\mathrm{'+names[column3] + '}$ / $\mathrm{' + units[column3]+'}$')
        print 'Mean of ' + names[column3] + ' : ' + str(round(np.mean(np.transpose(data)[column3][0:]),1)) + ' ' + units[column3]
        print 'standard deviation of sensor: ', np.std(np.transpose(data)[column3][0:])
    if ii == 1:
        ax2.plot(C3maxt,C3max[:,1], 'k--', linewidth=2, label='DLC 6.3 '+r'$V_\mathrm{e01}$, $+30$ $\mathrm{deg}$')
        ax2.plot(C3mint,C3min[:,1], 'k--', linewidth=2)
        ax2.legend(loc=5)
        ax2.grid('on')
        ax2.set_xlim([110,600])
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