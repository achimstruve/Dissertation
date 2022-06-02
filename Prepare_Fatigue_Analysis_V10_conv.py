from __future__ import division
import numpy as np
import shutil
import os
import sys
import time
import math
import matplotlib.pyplot as plt
from random import random
from operator import itemgetter
import configparser
from CS_coords_V02 import Leg_props, Hat_props
from Get_Coords_RectaR_V05 import creat_tower_nodes4

# Neues sechseckiges Eckstielprofiel
# Lese Querschnittseigenschaften aus den schon berechneten Querschnitten in Matlab
# parallelisierte Version

"""
# create batch files
ML_root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_cupled_to_SubDyn\\"
calcs = 48
N_members = 408
while np.mod(N_members,calcs) != 0:
    calcs -= 1

all_batch = ['' for xi in range(300)]
all_batch[0] = 'D:\n'
all_batch[1] = 'cd ' + ML_root[0:-1] + '\n'
for ci in range(calcs):
    batch = ['' for xi in range(300)]
    batch[0] = 'D:\n'
    batch[1] = 'cd ' + ML_root[0:-1] + '\n'
    batch[2] = 'python Prepare_Fatigue_Analysis_V10.py '+str(ci)+'\n'
    all_batch[2 + ci] = 'start /B Call _zpython_fatigue_instance'+str(ci+1)+'.bat\n'

    with open(ML_root + '_zpython_fatigue_instance'+str(ci+1)+'.bat', 'w') as file:
        file.writelines(batch)
        file.close()

with open(ML_root + '__python_fatigue_all.bat', 'w') as file:
    file.writelines(all_batch)
    file.close()
"""
"""
# search for remaining members
ML_root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_cupled_to_SubDyn\\"
All_files = os.listdir(ML_root + "Load_Results\Member_Fatigue\Member_Stresses_V23")
wind_speeds_names=np.array(['V04','V06','V08','V10','V11.4','V12','V14','V16','V18','V20','V22','V24'])
N_members = 408
remaining_files = []
for i in range(len(wind_speeds_names)):
    vwn = wind_speeds_names[i]
    for j in range(N_members):
        if (not os.path.exists(ML_root + "Load_Results\Member_Fatigue\Member_Stresses_V23\SM"+str(j+1)+"Dlc_11T_+00deg_"+vwn+"_S1.out")):
            remaining_files.append(str(j+1))

Nseg = 20
Nleg = 4
calcs = 48
while np.mod(N_members,calcs) != 0:
    calcs -= 1
N_members = (Nseg * Nleg * 5 + Nleg * 2)
zc = 0
zm = 0
Requested_Members_a = np.zeros([calcs,int(N_members/calcs)], dtype = int)
for mi in range(1, N_members + 1):
    Requested_Members_a[zc, zm] = mi
    if zm == 0:
        for i in range(len(remaining_files)):
            if mi == int(remaining_files[i]):
                print zc+1
                break
    zm += 1
    if zm >= int(N_members/calcs):
        zm = 0
        zc += 1
"""

def prepare_fatigue(ML_root, start_m_idx, tower_name, vers, calcs):
    
    root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_cupled_to_SubDyn\Out_Files\\"
    CS_root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_cupled_to_SubDyn\Cross_Sectional_Properties\\"
    SD_params_root = "D:\Struve\Promotion\Calculation\FAST_8.16\FAST_cupled_to_SubDyn\SubDyn\Rectangular_Tower_Single_V0\Tower_Results\\"
    
    PSF = 1.0
    
    All_files_folder = os.listdir(root)
    
    # create new directory for member stresses, if they not exist
    if not os.path.isdir(ML_root + 'Fatigue/Member_Stresses_V'+vers):
        os.makedirs(ML_root + 'Fatigue/Member_Stresses_V'+vers)
        print "\n"+ML_root + 'LFatigue/Member_Stresses_V'+vers + " has been created."
    
    ## all other parameters ##
    # read all parameters from parameter file
    conf = configparser.ConfigParser()
    conf.read(SD_params_root + tower_name + "_inputs.txt")

    Ht = 87.6
    rho = 8500
    E = 0.210000E+12
    G = 8.076923e+10
    
    Requested_Members = np.arange(1,12)

    
    #Requested_Members_a = np.array([np.arange(77,89),np.arange(77,89)])
    #Requested_Members = np.array([1])
    # get all fatigue relevant files
    """
    All_files = []
    for i in range(len(All_files_folder)):
        fn_tmp = All_files_folder[i]
        if (fn_tmp[0:15] == 'Dlc_11T_+00deg_') and (fn_tmp[-6:-4] != 'n2') and (fn_tmp[-5] != 'n') and (fn_tmp[-5] == '1'):
            All_files.append(fn_tmp)
    """
    All_files = ['Dlc_11T_+00deg_V04_S1.out', 'Dlc_11T_+00deg_V06_S1.out', 'Dlc_11T_+00deg_V08_S1.out', 'Dlc_11T_+00deg_V10_S1.out', 'Dlc_11T_+00deg_V11.4_S1.out', 'Dlc_11T_+00deg_V12_S1.out', 'Dlc_11T_+00deg_V14_S1.out', 'Dlc_11T_+00deg_V16_S1.out', 'Dlc_11T_+00deg_V18_S1.out', 'Dlc_11T_+00deg_V20_S1.out', 'Dlc_11T_+00deg_V22_S1.out', 'Dlc_11T_+00deg_V24_S1.out']       
    
    seed_names_11=np.array(['S1','S2','S3','S4','S5','S6'])
    #seed_names_11=np.array(['S1'])
    wind_speeds_11=np.array(['04','06','08','10','11.4','12','14','16','18','20','22','24'])
    All_files = ['' for xi in range(len(seed_names_11) * len(wind_speeds_11))]
    zs = 0
    for si in range(len(seed_names_11)):
        for wi in range(len(wind_speeds_11)):
          All_files[zs] = 'Dlc_11T_+00deg_V'+wind_speeds_11[wi]+'_'+seed_names_11[si]+'.out'
          zs += 1
    
    
    def mag(x):
        if x>0:
            return int(math.log10(x))
        if x<0:
            x = -x
            return int(math.log10(x))
        if x==0:
            return 0
    
    def fs(value, didgits):
        try:
            if value < 0:
                didgits -= 1
            if mag(value) < 0:
                didgits -= 1
            string_ = str.format('{0:.'+str(didgits)+'f}',(value)/(10**mag(value)))+"e"+str(mag(value))
        except:
            string_ = str(value)
        return string_
    
    
    # load cross-sectional parameters
    nodes_a, nodes_b, support_nodes, interface_nodes, leg_angles, nID_b, nID_ba, Cmass_a, NCmass = creat_tower_nodes4(At, Am, Ab, Bt, Bm, Bb, Ht, Htt, L1_max, L1_min, Nseg, Nsegt, Nleg, alpha_y, root, D_tower_top, precone, Overhang, Twr2Shft, NacCMxn, NacCMzn, mass_Nacelle, mass_hub, mass_blade, CM_blade_root, inertia_nacelle_yaw, inertia_hub_CM, CM_blade_downwind, fy, L2, rotatable, Damp_mass)
    CS_params = np.zeros([len(Requested_Members), 11], dtype = float) # cross sectional parameters
    m_CS_coords = np.zeros([len(Requested_Members), 12, 2], dtype = float) # cross sectional coordinates at which the stresses should be calculated
    zL = 0
    zS = 0
    for mi in range(len(Requested_Members)):
        CS_params[mi,10] = zS
        if mi == (Nseg+2) * Nleg:
            zL = 0
            zS = 0
        if mi < (Nseg+2) * Nleg:
            CS_params[mi,8] = 0 # shows, that it is a leg cross section
            if mi < ((Nseg-1) * Nleg):
                if (zL == 0) or (zL == 3):
                    CS_params[mi,9] = 0 # shows, that it is a leg cross section at the front of the tower
                if (zL == 1) or (zL == 2):
                    CS_params[mi,9] = 1 # shows, that it is a leg cross section at the back of the tower
                zL += 1
                if zL == 4:
                    zL = 0
            if mi >= ((Nseg-1) * Nleg):
                if mi == ((Nseg-1) * Nleg):
                    zL = 0
                if (zL == 0) or (zL == 1) or (zL == 2) or (zL == 9) or (zL == 10) or (zL == 11):
                    CS_params[mi,9] = 0 # shows, that it is a leg cross section at the front of the tower
                if (zL == 3) or (zL == 4) or (zL == 5) or (zL == 6) or (zL == 7) or (zL == 8):
                    CS_params[mi,9] = 1 # shows, that it is a leg cross section at the back of the tower
                zL += 1
            if zL == 0:
                zS += 1
        
        if mi >= (Nseg+2) * Nleg:        
            CS_params[mi,8] = 1 # shows, that it is a bracing cross section
            zL += 1
            if zL == 16:
                zL = 0
                zS += 1
    
    
    for mi in Requested_Members_a[start_m_idx,:]-1:
        zS = int(CS_params[mi,10])
        segment_hight = nodes_b[(zS) * Nleg, 3]
        base_hight = nodes_b[(Nseg-1) * Nleg, 3]
        H_F = round(c1_bot + (c1_top-c1_bot) / base_hight * segment_hight,0) / 1e3
        H_L = round(c2_bot + (c2_top-c2_bot) / base_hight * segment_hight,0) / 1e3
        #tU = round(tU_bot + (tU_top-tU_bot) / base_hight * segment_hight,0) / 1e3
        H_L = H_L * SFB_a[zS]
        tU = np.max([np.ceil((H_F * 1e3)/(ct_min_1support)*1.05), np.ceil((H_L * 1e3)/(2+ct_min_2support)*1.05)]) / 1e3
        if mi < (Nseg+2) * Nleg:        
            # rescale the leg dimension L_L of the leg profiles for the top segments
            if (CS_params[mi,9] == 0):
                SFL = SFLf_a[zS]
            if (CS_params[mi,9] == 1):
                SFL = SFLb_a[zS]
            L1 = leg_angles[zS,0,1]
            L_L = L1 * SFL
            A, Ix, Iy, Ixy, CM, coords = Leg_props(False, Nseg, E, G, rho, L_L, L2, tp, wp, dnL, dnt, H_L, H_F, tU, fy)
        
        if mi >= (Nseg+2) * Nleg:        
            A, Ix, Iy, Ixy, CM, coords, elements = Hat_props(False, E, G, rho, H_F, H_L, tU, wp, tp, dnL, dnt)    
        
        CS_params[mi,0] = A
        CS_params[mi,1] = CM[0] #xc
        CS_params[mi,2] = CM[1] #yc
        CS_params[mi,3] = Ix
        CS_params[mi,4] = Iy
        CS_params[mi,5] = Ixy
        if CS_params[mi,8] == 0:
            CS_params[mi,6] = L_L
            CS_params[mi,7] = L2
        if CS_params[mi,8] == 1:
            CS_params[mi,6] = H_L
            CS_params[mi,7] = H_F
        
        m_CS_coords[mi, 0:len(coords[:,0]),:] = coords
    
    """
    for i in range(len(m_CS_coords[:,0])):
        print i, m_CS_coords[i,:]
    """
    
    """
    # print chosen channel names
    for i in range(len(indexes)):
        print i, channels[indexes[i]]
    """
    
    # prepare Mlife input files
    # create matlab input file for all members
    mtlb_inp = ['\n' for x in range(int(1e4))]
    mtlb_inp[0] = 'close all\n'
    mtlb_inp[1] = 'clear all\n'
    mtlb_inp[2] = 'clc\n'
    mtlb_inp[3] = "addpath('D:\Struve\Programme\ASE\MLife\Source','D:\Struve\Programme\ASE\MLife\Source/rainflow','D:\Struve\Programme\ASE\MLife\Source\datatablepackage');\n"
    
    # read raw input file
    with open(ML_root + "Input_MLife_Lattice_raw_only_Seed1.mlif", 'r') as file:
        ML_ini = file.readlines()
        file.close()
    
    spots = 8
    for mi in range(len(Requested_Members)):
        ML_mi = ['\n' for x in range(int(1e4))]
        ML_mi[0:len(ML_ini)] = np.copy(ML_ini)
        ML_mi[7] = '"Results/Dlc12_V'+vers+'_M'+str(Requested_Members[mi])+'_only_Seed1"      	  RootName          Root name for output files.\n'
    
        counter = 0
        for Ji in range(2):      
            for spti in range(spots):
                ML_mi[56 + counter] = '$M'+str(Requested_Members[mi])+'szJ'+str(Ji+1)+str(spti+1)+'$       1    4              BN		20      	WM      3007682577.780582\n'
                counter += 1
    
        ML_mi[32] = str(int(counter)) + '               nFatigueChannels    The number of rainflow channels.  Next six lines ignored if zero.\n'
        ML_mi[56 + counter] = '1                 nGroups             Number of fatigue groups\n'
        ML_mi[57 + counter] = 'Group Name  NChannels      ChannelList\n'
        string = ''
        for i in range(1,counter+1):
            string += str(i) + ' '
        ML_mi[57 + counter + 1] = '"M'+str(Requested_Members[mi])+'"            '+str(int(counter))+'              '+string+'\n'
        ML_mi[57 + counter + 2] = '-----  Input Files  ------------------------------------------------------------\n'
        ML_mi[57 + counter + 3] = '1                 FileFormat         Format of input files.  1 = FAST ascii, 2 = FAST binary\n'
        ML_mi[57 + counter + 4] = str(len(All_files[:])) + '  1.1   1.3   1.5   1.7    (Weibull-Weighted Normal Operation: NumNormFiles, PSF1, PSF2, PSF3, PSF4)\n'
        
        counter2 = 0
        for i in range(len(All_files)):
            ML_mi[57 + counter + 5 + counter2] = '"Member_Stresses_V'+vers+'/SM'+str(Requested_Members[mi]) + All_files[i] + '"\n'
            counter2 += 1
        
        ML_mi[57 + counter + 5 + counter2] = '0  1.1   1.3   1.5   1.7    (Weibull-Weighted Idling: NumIdleFiles, PSF1, PSF2, PSF3, PSF4)\n'
        ML_mi[57 + counter + 5 + counter2 + 1] = '0  1.2   1.3   1.4   1.6    (Discrete Events: NumDiscFiles, PSF1, PSF2, PSF3, PSF4)\n'
        ML_mi[57 + counter + 5 + counter2 + 2] = '==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.\n'
        
        with open(ML_root + "Load_Results/Member_Fatigue/Input_MLife_Lattice_M"+str(Requested_Members[mi])+"_only_Seed1.mlif", 'w') as file:
            file.writelines(ML_mi)
            file.close()
        
        mtlb_inp[4 + mi * 2] = "settingsFile"+str(Requested_Members[mi])+" = 'Input_MLife_Lattice_M"+str(Requested_Members[mi])+"_only_Seed1.mlif';\n"
        mtlb_inp[4 + mi * 2 + 1] = "[Fatigue"+str(Requested_Members[mi])+", Statistics"+str(Requested_Members[mi])+"]= mlife( settingsFile"+str(Requested_Members[mi])+" );\n"
    
    """
    mtlb_inp[4 + len(Requested_Members) * 2 + 2] = 'n = 407\n'
    mtlb_inp[4 + len(Requested_Members) * 2 + 3] = 'fatigue_DMG = zeros(n,16)\n'
    mtlb_inp[4 + len(Requested_Members) * 2 + 4] = 'for i=1:n\n'
    for mi in range(len(Requested_Members)):
        mtlb_inp[4 + len(Requested_Members) * 2 + 5 + mi] = "fatigue_DMG("+str(mi)+",1:length(Fatigue"+str(Requested_Members[mi])+".Channel.lifetimeDamage    \n"
    """
    with open(ML_root + "Load_Results/Member_Fatigue/Fatigue_MLife_V01.m", 'w') as file:
        file.writelines(mtlb_inp)
        file.close()
    
    
    # save loads of all requested members to new array
    tstps = int(600/0.1 + 1)
    idx_start = int(0/0.1)
    
    # Get column indexes of the requested member channels
    channels_string = 'YawBrFxp  	YawBrFyp  	YawBrFzp  	YawBrMxp  	YawBrMyp  	YawBrMzp  	TwHt9FLxt 	TwHt9FLyt 	TwHt9FLzt 	TwHt9MLxt 	TwHt9MLyt 	TwHt9MLzt 	TwHt8FLxt 	TwHt8FLyt 	TwHt8FLzt 	TwHt8MLxt 	TwHt8MLyt 	TwHt8MLzt 	TwHt7FLxt 	TwHt7FLyt 	TwHt7FLzt 	TwHt7MLxt 	TwHt7MLyt 	TwHt7MLzt 	TwHt6FLxt 	TwHt6FLyt 	TwHt6FLzt 	TwHt6MLxt 	TwHt6MLyt 	TwHt6MLzt 	TwHt5FLxt 	TwHt5FLyt 	TwHt5FLzt 	TwHt5MLxt 	TwHt5MLyt 	TwHt5MLzt 	TwHt4FLxt 	TwHt4FLyt 	TwHt4FLzt 	TwHt4MLxt 	TwHt4MLyt 	TwHt4MLzt 	TwHt3FLxt 	TwHt3FLyt 	TwHt3FLzt 	TwHt3MLxt 	TwHt3MLyt 	TwHt3MLzt 	TwHt2FLxt 	TwHt2FLyt 	TwHt2FLzt 	TwHt2MLxt 	TwHt2MLyt 	TwHt2MLzt 	TwHt1FLxt 	TwHt1FLyt 	TwHt1FLzt 	TwHt1MLxt 	TwHt1MLyt 	TwHt1MLzt 	TwrBsFxt  	TwrBsFyt  	TwrBsFzt  	TwrBsMxt  	TwrBsMyt  	TwrBsMzt'
    channels = channels_string.split
    
    for i in range(len(All_files)):
        print "Progress: ", round(i/(len(All_files))*100,2), " %"
        with open(root + All_files[i], 'r') as file:
            data_i = file.readlines()
            file.close()
        
        indexes = np.zeros([2*6], dtype = int)
        for mi in Requested_Members_a[start_m_idx,:]-1:
            midx = Requested_Members[mi]
            cidx = 0
            for ci in range(len(channels)):
                c_tmp = channels[ci]
                if (c_tmp[0:mag(midx)+6] == 'M'+str(midx)+'J1FK') or (c_tmp[0:mag(midx)+6] == 'M'+str(midx)+'J2FK') or (c_tmp[0:mag(midx)+6] == 'M'+str(midx)+'J1MK') or (c_tmp[0:mag(midx)+6] == 'M'+str(midx)+'J2MK'):
                    indexes[cidx] = ci
                    cidx += 1
            stresses = np.zeros([tstps, 2, 8])
            new_data = np.zeros([tstps, len(indexes)+4])
            #print mi/len(Requested_Members)
            ## Choose the correct cross sectional parameter set for this member
            
    
            A = CS_params[mi,0]
            xc = CS_params[mi,1]
            yc = CS_params[mi,2]
            Ix = CS_params[mi,3]
            Iy = CS_params[mi,4]
            Ixy = CS_params[mi,5]
            
            if CS_params[mi,8] == 0:
                L_L = CS_params[mi,6]
                L_F = CS_params[mi,7]
                
                p1 = m_CS_coords[mi,0]
                p3 = m_CS_coords[mi,2]
                p4 = m_CS_coords[mi,3]
                p6 = m_CS_coords[mi,5]
                p7 = m_CS_coords[mi,6]
                p8 = m_CS_coords[mi,7]
                p10 = m_CS_coords[mi,9]
                p11 = m_CS_coords[mi,10]
                
            if CS_params[mi,8] == 1:
                H_L = CS_params[mi,6]
                H_F = CS_params[mi,7]
    
                p1 = m_CS_coords[mi,0]
                p3 = m_CS_coords[mi,2]
                p4 = m_CS_coords[mi,3]
                p6 = m_CS_coords[mi,5]
                p7 = m_CS_coords[mi,6]
                p9 = m_CS_coords[mi,8]
                p10 = m_CS_coords[mi,9]
                p12 = m_CS_coords[mi,11]
    
            
            for j in range(len(data_i[8+idx_start:])):
                new_data[j,0:3] = data_i[j+8+idx_start].split('\t')[0:3] # time and wind speed
                new_data[j,3] = np.sqrt(new_data[j,1]**2 + new_data[j,2]**2) # resultand wind speed
                new_data[j,4:] = itemgetter(*indexes.tolist())(data_i[j+8+idx_start].split('\t')) # loads
    
                # get loads of this member at both joints
                FxJ1 = new_data[j, 4]
                FyJ1 = new_data[j, 5]
                FzJ1 = new_data[j, 6]
                MxJ1 = new_data[j, 7]
                MyJ1 = new_data[j, 8]
                MzJ1 = new_data[j, 9]
                
                FxJ2 = new_data[j, 10]
                FyJ2 = new_data[j, 11]
                FzJ2 = new_data[j, 12]
                MxJ2 = new_data[j, 13]
                MyJ2 = new_data[j, 14]
                MzJ2 = new_data[j, 15]
                    
                
                ### calculate the stresses ###
                
                """
                # account for bracing eccentricities
                if ms > 0:
                    # calculate eccentricity for this tower segment with respect to the leg center of mass
                    if ms == 1:
                        iseg = 0
                    if ms == 2:
                        iseg = 9
                    if ms == 3:
                        iseg = 19
                    e_Leg_CM = e_Leg_CM_Bot + (e_Leg_CM_Top - e_Leg_CM_Bot) / (Nseg-1) * iseg
                    
                    #print Requested_Members[mi],MyJ1,FzJ1,e_Leg_CM, CSL_params[1], MyJ1/(FzJ1 * (e_Leg_CM - (-CSL_params[1])))
                    
                    if ecc_counter < 2:
                        MyJ1 = MyJ1 + FzJ1 * (e_Leg_CM - (-CSL_params[1]))
                        MyJ2 = MyJ2 + FzJ2 * (e_Leg_CM - (-CSL_params[1]))
                    if ecc_counter >= 2:
                        MyJ1 = MyJ1 + FzJ1 * (e_Leg_CM - (-CSL_params[1]))
                        MyJ2 = MyJ2 + FzJ2 * (e_Leg_CM - (-CSL_params[1]))
                    
                    ecc_counter += 1
                    if ecc_counter == 4:
                        ecc_counter = 0
                    
                    #print FzJ1 / CSL_params[0], (MxJ1 / CSL_params[3] * (CSL_params[9]/2.0))/(FzJ1 / CSL_params[0]), (MyJ1 / CSL_params[4] * (-CSL_params[1] - CSL_params[9]))/(FzJ1 / CSL_params[0])
                    #print mi,"\t", round(MyJ1 / CSL_params[4] * (-CSL_params[1] - CSL_params[9])/1e6,1),"\t", round(FzJ1 / CSL_params[0]/1e6,1)
                """
                
                
                if CS_params[mi,8] == 0:
                    ## normal stresses ##
                    # J1 spot 1
                    stresses[j, 0, 0] = FzJ1 / A + MxJ1 / Ix * (p1[1] - yc) + MyJ1 / Iy * (p1[0] - xc)
                    # J1 spot 2
                    stresses[j, 0, 1] = FzJ1 / A + MxJ1 / Ix * (p3[1] - yc) + MyJ1 / Iy * (p3[0] - xc)
                    # J1 spot 3
                    stresses[j, 0, 2] = FzJ1 / A + MxJ1 / Ix * (p4[1] - yc) + MyJ1 / Iy * (p4[0] - xc)
                    # J1 spot 4
                    stresses[j, 0, 3] = FzJ1 / A + MxJ1 / Ix * (p6[1] - yc) + MyJ1 / Iy * (p6[0] - xc)
                    # J1 spot 5
                    stresses[j, 0, 4] = FzJ1 / A + MxJ1 / Ix * (p7[1] - yc) + MyJ1 / Iy * (p7[0] - xc)
                    # J1 spot 6
                    stresses[j, 0, 5] = FzJ1 / A + MxJ1 / Ix * (p8[1] - yc) + MyJ1 / Iy * (p8[0] - xc)
                    # J1 spot 7
                    stresses[j, 0, 6] = FzJ1 / A + MxJ1 / Ix * (p10[1] - yc) + MyJ1 / Iy * (p10[0] - xc)
                    # J1 spot 8
                    stresses[j, 0, 7] = FzJ1 / A + MxJ1 / Ix * (p11[1] - yc) + MyJ1 / Iy * (p11[0] - xc)
                    # J2 spot 1
                    stresses[j, 1, 0] = FzJ2 / A + MxJ2 / Ix * (p1[1] - yc) + MyJ2 / Iy * (p1[0] - xc)
                    # J2 spot 2
                    stresses[j, 1, 1] = FzJ2 / A + MxJ2 / Ix * (p3[1] - yc) + MyJ2 / Iy * (p3[0] - xc)
                    # J2 spot 3
                    stresses[j, 1, 2] = FzJ2 / A + MxJ2 / Ix * (p4[1] - yc) + MyJ2 / Iy * (p4[0] - xc)
                    # J2 spot 4
                    stresses[j, 1, 3] = FzJ2 / A + MxJ2 / Ix * (p6[1] - yc) + MyJ2 / Iy * (p6[0] - xc)
                    # J2 spot 5
                    stresses[j, 1, 4] = FzJ2 / A + MxJ2 / Ix * (p7[1] - yc) + MyJ2 / Iy * (p7[0] - xc)
                    # J2 spot 6
                    stresses[j, 1, 5] = FzJ2 / A + MxJ2 / Ix * (p8[1] - yc) + MyJ2 / Iy * (p8[0] - xc)
                    # J2 spot 7
                    stresses[j, 1, 6] = FzJ2 / A + MxJ2 / Ix * (p10[1] - yc) + MyJ2 / Iy * (p10[0] - xc)
                    # J2 spot 8
                    stresses[j, 1, 7] = FzJ2 / A + MxJ2 / Ix * (p11[1] - yc) + MyJ2 / Iy * (p11[0] - xc)
              
                if CS_params[mi,8] == 1:
                    ## normal stresses ##
                    # J1 spot 1
                    stresses[j, 0, 0] = FzJ1 / A + MxJ1 / Ix * (p1[1] - yc) + MyJ1 / Iy * (p1[0] - xc)
                    # J1 spot 2
                    stresses[j, 0, 1] = FzJ1 / A + MxJ1 / Ix * (p3[1] - yc) + MyJ1 / Iy * (p3[0] - xc)
                    # J1 spot 3
                    stresses[j, 0, 2] = FzJ1 / A + MxJ1 / Ix * (p4[1] - yc) + MyJ1 / Iy * (p4[0] - xc)
                    # J1 spot 4
                    stresses[j, 0, 3] = FzJ1 / A + MxJ1 / Ix * (p6[1] - yc) + MyJ1 / Iy * (p6[0] - xc)
                    # J1 spot 5
                    stresses[j, 0, 4] = FzJ1 / A + MxJ1 / Ix * (p7[1] - yc) + MyJ1 / Iy * (p7[0] - xc)
                    # J1 spot 6
                    stresses[j, 0, 5] = FzJ1 / A + MxJ1 / Ix * (p9[1] - yc) + MyJ1 / Iy * (p9[0] - xc)
                    # J1 spot 7
                    stresses[j, 0, 6] = FzJ1 / A + MxJ1 / Ix * (p10[1] - yc) + MyJ1 / Iy * (p10[0] - xc)
                    # J1 spot 8
                    stresses[j, 0, 7] = FzJ1 / A + MxJ1 / Ix * (p12[1] - yc) + MyJ1 / Iy * (p12[0] - xc)
                    # J2 spot 1
                    stresses[j, 1, 0] = FzJ2 / A + MxJ2 / Ix * (p1[1] - yc) + MyJ2 / Iy * (p1[0] - xc)
                    # J2 spot 2
                    stresses[j, 1, 1] = FzJ2 / A + MxJ2 / Ix * (p3[1] - yc) + MyJ2 / Iy * (p3[0] - xc)
                    # J2 spot 3
                    stresses[j, 1, 2] = FzJ2 / A + MxJ2 / Ix * (p4[1] - yc) + MyJ2 / Iy * (p4[0] - xc)
                    # J2 spot 4
                    stresses[j, 1, 3] = FzJ2 / A + MxJ2 / Ix * (p6[1] - yc) + MyJ2 / Iy * (p6[0] - xc)
                    # J2 spot 5
                    stresses[j, 1, 4] = FzJ2 / A + MxJ2 / Ix * (p7[1] - yc) + MyJ2 / Iy * (p7[0] - xc)
                    # J2 spot 6
                    stresses[j, 1, 5] = FzJ2 / A + MxJ2 / Ix * (p9[1] - yc) + MyJ2 / Iy * (p9[0] - xc)
                    # J2 spot 7
                    stresses[j, 1, 6] = FzJ2 / A + MxJ2 / Ix * (p10[1] - yc) + MyJ2 / Iy * (p10[0] - xc)
                    # J2 spot 8
                    stresses[j, 1, 7] = FzJ2 / A + MxJ2 / Ix * (p12[1] - yc) + MyJ2 / Iy * (p12[0] - xc)
    
                # apply the partial safety factor to the calculated stresses
                stresses[j, :, :] = stresses[j, :, :] * PSF
    
            # write stresses to new file
            stress_file = ['' for xx in range(4 + tstps)]
            stress_file[0] = 'This file contains self calculated stresses for certain members. It has been generated out of load time series file ' + All_files[i] + '\n'
            stress_file[1] = 'Abbreviations: sz = sigma_z, tV = transverse shear, tT = torsional shear, sM = sigma_mises, phi = main stress angle, Jxy = Joint x detail y\n'
            # generate channel names
            for cni in range(1 + 1):
                if cni == 0:
                    stress_file[2] = 'Time      \tWind1VelX \tWind1VelY \tWind1VelXY'
                    stress_file[3] = '(s)       \t(m/s)     \t(m/s)     \t(m/s)     '
                if cni > 0:
                    mn = '\tM' + str(Requested_Members[mi])
                    m = Requested_Members[mi]
                    stress_file[2] = stress_file[2] + mn + 'szJ11' + (3-mag(m))*' ' + mn + 'szJ12' + (3-mag(m))*' ' + mn + 'szJ13' + (3-mag(m))*' ' + mn + 'szJ14' + (3-mag(m))*' ' + mn + 'szJ15' + (3-mag(m))*' ' + mn + 'szJ16' + (3-mag(m))*' ' + mn + 'szJ17' + (3-mag(m))*' ' + mn + 'szJ18' + (3-mag(m))*' ' + mn + 'szJ21' + (3-mag(m))*' ' + mn + 'szJ22' + (3-mag(m))*' ' + mn + 'szJ23' + (3-mag(m))*' ' + mn + 'szJ24' + (3-mag(m))*' ' + mn + 'szJ25' + (3-mag(m))*' ' + mn + 'szJ26' + (3-mag(m))*' ' + mn + 'szJ27' + (3-mag(m))*' ' + mn + 'szJ28' + (3-mag(m))*' ' 
                    stress_file[3] = stress_file[3] + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' ' + '\t(N/m^2)' + 3*' '
            
            stress_file[2] = stress_file[2] + '\n'
            stress_file[3] = stress_file[3] + '\n'
            for ti in range(tstps):
                stress_file[4+ti] = stress_file[4+ti] + fs(new_data[ti, 0], 4) + (2)*' '
                stress_file[4+ti] = stress_file[4+ti] +  '\t' + fs(new_data[ti, 1], 4) + (2)*' '
                stress_file[4+ti] = stress_file[4+ti] +  '\t' + fs(new_data[ti, 2], 4) + (2)*' '
                stress_file[4+ti] = stress_file[4+ti] +  '\t' + fs(new_data[ti, 3], 4) + (2)*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' + fs(stresses[ti, 0, 0], 4) + (2-mag(stresses[ti, 0, 0]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 1], 4) + (2-mag(stresses[ti, 0, 1]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 2], 4) + (2-mag(stresses[ti, 0, 2]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 3], 4) + (2-mag(stresses[ti, 0, 3]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 4], 4) + (2-mag(stresses[ti, 0, 4]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 5], 4) + (2-mag(stresses[ti, 0, 5]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 6], 4) + (2-mag(stresses[ti, 0, 6]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 0, 7], 4) + (2-mag(stresses[ti, 0, 7]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 0], 4) + (2-mag(stresses[ti, 1, 0]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 1], 4) + (2-mag(stresses[ti, 1, 1]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 2], 4) + (2-mag(stresses[ti, 1, 2]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 3], 4) + (2-mag(stresses[ti, 1, 3]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 4], 4) + (2-mag(stresses[ti, 1, 4]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 5], 4) + (2-mag(stresses[ti, 1, 5]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 6], 4) + (2-mag(stresses[ti, 1, 6]))*' '
                stress_file[4+ti] = stress_file[4+ti] + '\t' +  fs(stresses[ti, 1, 7], 4) + (2-mag(stresses[ti, 1, 7]))*' '
                
                stress_file[4+ti] = stress_file[4+ti] + '\n'
            
            with open(ML_root + 'Load_Results/Member_Fatigue/Member_Stresses_V'+vers+'/SM'+str(mi+1)+All_files[i], 'w') as file:
                file.writelines(stress_file)
                file.close()
            
    
            stress_file = []
            stresses = []
            new_data = []
    
        data_i = []
    
    """
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
    mi = 0
    ax1.plot(new_data[:,0], stresses[:, mi * 5 * 2, 0, 0] / 1e6, color = 'r', label = 'normal stress')
    ax1.plot(new_data[:,0], stresses[:, mi * 5 * 2 + 1, 0, 0] / 1e6, color = 'b', label = 'transverse shear stress')
    ax1.plot(new_data[:,0], stresses[:, mi * 5 * 2 + 2, 0, 0] / 1e6, color = 'g', label = 'torsional shear stress')
    color_a = ['r','b','g','c']
    for i in range(4):
        ax2.plot(new_data[:,0], stresses[:, i * 5 * 2 + 3, 0, 0] / 1e6, color = color_a[i], label = 'mises stress - member '+str(i+1))
    ax3.plot(new_data[:,0], stresses[:, mi * 5 * 2 + 4, 0, 0] * 180 / np.pi, color = 'm', label = 'main stress angle')
    ax1.legend()
    ax1.grid('on')
    ax2.legend()
    ax2.grid('on')
    ax3.legend()
    ax3.grid('on')
    """