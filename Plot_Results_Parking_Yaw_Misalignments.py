import os
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

root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_Assessment\Out_Files\\"

All_files_folder = os.listdir(root)

# get all extreme event relevant files
All_files = []
for i in range(len(All_files_folder)):
    fn_tmp = All_files_folder[i]
    if (fn_tmp[0:7] == 'Dlc_62T') and not (fn_tmp[-1] == 'm'):
        All_files.append(fn_tmp)

directions_62=np.array([-180,-160,-140,-120,-100,-80,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,80,100,120,140,160,180])
res1 = np.zeros([len(directions_62),6,4])
res1b = np.zeros([len(directions_62),4])
res2 = np.zeros([len(directions_62),6,4])
res2b = np.zeros([len(directions_62),4])
for ii in range(len(All_files)):
    print round(ii*1.0/(len(All_files)*1.0) * 100,1)
    filename = All_files[ii]
    column2 = 94
    column3 = 117

    # get current wind direction from filename
    if len(filename) == 22:
        winddir = float(filename[8:12])
    if len(filename) == 21:
        winddir = float(filename[8:11])

    # load data from file
    f=open(root+filename,'r') # open file
    names=list(f)[6].split('\t') # laod the list of output variables
    names[-1] = names[-1][0:-1]
    f=open(root+filename,'r') # open file
    units=list(f)[7].split('\t') # laod the list of output variables
    units[-1] = units[-1][0:-1]
    f=open(root+filename,'r') # open file
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


    column1 = 0
    # plot der entsprechenden Datenspalten
    params = {'legend.fontsize': 28}
    plt.rcParams.update(params)
    font = {'family' : 'normal',
            'size'   : 30}
    plt.rc('font', **font)
    
    #ax1.scatter(np.transpose(data)[column1][0:], np.max(data[column2][0:]), color = 'r', label = r'$\mathrm{'+names[column2] + '}$ / $\mathrm{' + units[column2] +'}$')
    iwd = np.where(winddir == directions_62[:])
    si = int(filename[-5])-1
    data2 = np.transpose(data)[column2][0:] / 1e3
    data3 = np.transpose(data)[column3][0:] / 1e3
    res1[iwd,si,0] = np.max(data2)
    res1[iwd,si,1] = np.mean(data2)
    res1[iwd,si,2] = np.min(data2)
    res1[iwd,si,3] = np.std(data2)
    res2[iwd,si,0] = np.max(data3)
    res2[iwd,si,1] = np.mean(data3)
    res2[iwd,si,2] = np.min(data3)
    res2[iwd,si,3] = np.std(data3)

    if ii == 0:
        fig, (ax1,ax2) = plt.subplots(2, sharex=True)

    ax1.scatter(winddir, np.max(data2), color = 'r')
    ax1.scatter(winddir, np.mean(data2), color = 'g')
    ax1.scatter(winddir, np.min(data2), color = 'b')
    ax2.scatter(winddir, np.max(data3), color = 'r')
    ax2.scatter(winddir, np.mean(data3), color = 'g')
    ax2.scatter(winddir, np.min(data3), color = 'b')
    
    if winddir == 180:
        ax1.scatter(-winddir, np.max(data2), color = 'r')
        ax1.scatter(-winddir, np.mean(data2), color = 'g')
        ax1.scatter(-winddir, np.min(data2), color = 'b')
        ax2.scatter(-winddir, np.max(data3), color = 'r')
        ax2.scatter(-winddir, np.mean(data3), color = 'g')
        ax2.scatter(-winddir, np.min(data3), color = 'b')
        res1[0,si,0] = np.max(data2)
        res1[0,si,1] = np.mean(data2)
        res1[0,si,2] = np.min(data2)
        res1[0,si,3] = np.std(data2)
        res2[0,si,0] = np.max(data3)
        res2[0,si,1] = np.mean(data3)
        res2[0,si,2] = np.min(data3)
        res2[0,si,3] = np.std(data3)
    
    if ii == (len(All_files)-1):
        for wdi in range(len(directions_62)):
            res1b[wdi,0] = np.max(res1[wdi,:,0])
            res1b[wdi,1] = np.mean(res1[wdi,:,1])
            res1b[wdi,2] = np.min(res1[wdi,:,2])
            res1b[wdi,3] = np.mean(res1[wdi,:,3])
            res2b[wdi,0] = np.max(res2[wdi,:,0])
            res2b[wdi,1] = np.mean(res2[wdi,:,1])
            res2b[wdi,2] = np.min(res2[wdi,:,2])
            res2b[wdi,3] = np.mean(res2[wdi,:,3])
        lns1 = ax1.plot(directions_62, res1b[:,0], color = 'r', linewidth=2, label = r'$\mathrm{Max.}$')
        lns2 = ax1.plot(directions_62, res1b[:,1], color = 'g', linewidth=2, label = r'$\mathrm{Mean}$')
        lns3 = ax1.plot(directions_62, res1b[:,2], color = 'b', linewidth=2, label = r'$\mathrm{Min.}$')
        ax11 = ax1.twinx()
        lns4 = ax11.plot(directions_62, res1b[:,3], color = 'k', linewidth=2, label=r'$\mathrm{\bar{\sigma}}$')
        ax11.set_ylabel(r'$\mathrm{\bar{\sigma}}$ / $\mathrm{MN}$', color='k')
        lns12 = ax2.plot(directions_62, res2b[:,0], color = 'r', linewidth=2, label = r'$\mathrm{Max.}$')
        lns22 = ax2.plot(directions_62, res2b[:,1], color = 'g', linewidth=2, label = r'$\mathrm{Mean}$')
        lns32 = ax2.plot(directions_62, res2b[:,2], color = 'b', linewidth=2, label = r'$\mathrm{Min.}$')
        ax21 = ax2.twinx()
        lns42 = ax21.plot(directions_62, res2b[:,3], color = 'k', linewidth=2, label=r'$\mathrm{\bar{\sigma}}$')
        ax21.set_ylabel(r'$\mathrm{\bar{\sigma}}$ / $\mathrm{MNm}$', color='k')
        ax1.set_xlim([-181,181])
        ax11.set_xlim([-181,181])
        ax2.set_xlim([-181,181])
        ax21.set_xlim([-181,181])
        ax1.set_ylim([-1,1])
        ax11.set_ylim([0,np.max(res1b[:,3])*1.02])
        ax1.grid('on')
        ax2.grid('on')
        lns1 = lns1+lns2+lns3+lns4
        labs1 = [l.get_label() for l in lns1]
        ax1.legend(lns1, labs1, loc=1, ncol=2)
        lns2 = lns12+lns22+lns32+lns42
        labs2 = [l.get_label() for l in lns2]
        ax2.legend(lns2, labs2, loc=1, ncol=2)
        ax1.set_ylabel(r'$\mathrm{'+names[column2] + '}$ / $\mathrm{MN}$')
        ax2.set_xlabel(r'$\mathrm{Yaw}$ $\mathrm{Error}$ / $\mathrm{deg}$')
        ax2.set_ylabel(r'$\mathrm{'+names[column3] + '}$ / $\mathrm{MNm}$')
        ax2.set_xticks(np.arange(-180,200,20))
    

