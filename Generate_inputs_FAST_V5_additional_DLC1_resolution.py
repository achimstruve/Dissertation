from __future__ import division
import numpy as np
import shutil
import os
from subprocess import *
from scipy.interpolate import *
import time
from random import random

# startet at 17:25 07.03.2018
# set root path
root_FAST="D:\Struve\Promotion\Calculation\FAST_8.16\FAST_Assessment\\"
root_TurbSim="D:\Struve\Promotion\Calculation\FAST_8.16\FAST_Assessment\InflowWind\Wind\\"
root_IECwind="D:\Struve\Promotion\Calculation\FAST_8.16\IEC_Wind\\"

# set available cores
cores=48
sleep_time=5.5
iteration=1 # do not touch this for normal studies

# Do MLife preparation?
Move_Out = 'False'
Calculate_bts = 'False'
Calculate_IEC = 'False'
Calculate_all_windfiles = 'False'
Calc_Rest_FAST = 'False'
Calc_Rest_TurbSim = 'False'


if Move_Out == 'False' and Calculate_bts == 'False' and Calculate_IEC == 'False' and Calculate_all_windfiles == 'False' and Calc_Rest_FAST == 'False' and Calc_Rest_TurbSim == 'False':
    ### set initial conditions ###
    OoPDefl=0                       # initial blade tip deflection of blade 1; (negativ = prebending direction) / m
    
    ### set output conditions ###
    SttsTime=10                     # time between screen status messages / s
    DT_out=0.05                  # time step for tabular output / s
    TStart=30                        # time to begin tabular output / s
    TStart_parked = 100             # time to begin tabular output for parked DLCs / s
    TStart_trans = 100             # time to begin tabular output for parked DLCs / s
    
    
    
    ### set wind turbine conditions ###
    turb_intens='B'
    wind_class=1
    pwr_law_exp = 0.2
    pwr_law_exp_idling = 0.11
    D_rot=126                       # rotor diameter / m
    towerHt=87.6                   # tower height / m
    hubHt=towerHt+2.4               # hub height / m
    overHang=-5.0191                # overhang is the distance between yaw axis and Hub C.M. negative for upwind turbine / m
    Vr=11.4                         # rated wind speed / m/s
    V_out=25                        # cut out wind speed / m/s
    inertia_nacelle_yaw = 2607890   # nacell yaw inertia / kg*m^2
    
    
    ### set FAST conditions ###
    DT_FAST=0.0125                           # time step for FAST simulation
    TMax_FAST=600 + TStart                   # duration for one FAST simulation
    TMax_FAST_trans=100  + TStart_trans            # maximum simulation time for transient DLCs in FAST
    TMax_idling = 600 + TStart_parked
    if DT_out<DT_FAST:
        DT_out=DT_FAST
    

    
    ### set DLC 3.1 and DLC 4.1 conditions ###
    SpdGenOn=100                    # generator speed where it goes on (interesting for start ups of DLC 3.1) from Definition of NREL 5MW wind turbine page 19 670 is given / rpm
    n_Vin31=1000
    n_Vr31=50
    n_Vout31=50
    n_Vin41=1000
    n_Vr41=50
    n_Vout41=50
    PitManRat=2       # PitManRat / deg/s
    
    ### set idling conditions
    GBoxEff = 98
    
    ### set TurbSim conditions ###
    DT_TurbSim=0.05                # time step for TurbSim velocity matricies
    AnalysisTime_TurbSim=600 + TStart        # duration for TurbSim timeseries
    AnalysisTime_TurbSim_idling= 600 + TStart_parked       # duration for TurbSim timeseries    
    
    GridHeight=np.ceil(2*(D_rot/2+np.abs(overHang))*1.15)  # took overHang instead of shaft length for the secure side
    GridWidth=np.ceil(2*(D_rot/2+np.abs(overHang))*1.15)   # took overHang instead of shaft length for the secure side
    
    # calculate "Turbulenzlaengenparameter" Lambda1 from IEC 61400-1:2011-08 equation 5
    if hubHt<60:
        Lambda1=hubHt*0.7
    if hubHt>=60:
        Lambda1=42
    
    # calculate maximum distance between TurbSim nodes
    ddn_TS_max1=D_rot*0.15
    ddn_TS_max2=Lambda1*0.25
    if ddn_TS_max1<=ddn_TS_max2:
        ddn_TS_max=ddn_TS_max1
    if ddn_TS_max2<ddn_TS_max1:
        ddn_TS_max=ddn_TS_max2
    # ddn_TS_max is the maximum diagonal distance between grid nodes of TurbSim from IEC 61400-1:2011-08 7.5 footnote 9)
    dn_TS_max=np.sin(np.pi/4)*ddn_TS_max
    
    NumGrid_Z=np.ceil(GridHeight/dn_TS_max)
    NumGrid_Y=np.ceil(GridWidth/dn_TS_max)
    
    
    
    ### set IECwind conditions
    DT_IECwind=0.0125
    IEC_trans_cond= 60
    
    ### generate files ###    
    # set standard filename
    filename11 = "Dlc_11"
    filename13 = "Dlc_13"
    filename14 = "Dlc_14"
    filename15 = "Dlc_15"
    filename21 = "Dlc_21"
    filename23 = "Dlc_23"
    filename31 = "Dlc_31"
    filename32 = "Dlc_32"
    filename33 = "Dlc_33"
    filename41 = "Dlc_41"
    filename42 = "Dlc_42"
    filename51 = "Dlc_51"
    filename61 = "Dlc_61"
    filename62 = "Dlc_62"
    filename63 = "Dlc_63"
    filename64 = "Dlc_64"
    filename71 = "Dlc_71"
    IEC_Input = "IEC"
    IEC_ending = ".ipt"
    TurbSim_Input = "TurbSim_Input"
    Input_ending = ".inp"
    
    # set filname of Batch files for Turb_Sim and FAST
    filename_batch_all_FAST="_Run_all_FAST"
    filename_batch_all_TS="_Run_all_TurbSim"
    filename_batch_all_IEC="_Run_all_IECwind"
    filename_batch_ending=".bat"
    
    # set filename of FAST files
    filename_FAST_end=".fst"
    
    # set filename of inflowWind files
    filename_inflow="NREL_5MW_InflowWind"
    filename_inflow_end=".dat"

    # set filename of ElastoDyn files
    filename_aero="NREL_5MW_AeroDyn15"
    filename_aero_end=".dat"

    # set filename of ElastoDyn files
    filename_elasto="NREL_5MW_ElastoDyn"
    filename_elasto_end=".dat"
    
    # set filename of ServoDyn files
    filename_servo="NREL_5MW_ServoDyn"
    filename_servo_end=".dat"

    # filesets for DLC 1.1
    seed_names_11=np.array(['S1'])
   
    # filesets for DLC 1.1
    wind_speeds_names_11=np.array(['V04','V06','V08','V10','V11.4','V12','V14','V16','V18','V20','V22','V24'])
    wind_speeds_11=np.array(['4.0','6.0','8.0','10.0','11.4','12.0','14.0','16.0','18.0','20.0','22.0','24.0'])
    wind_speeds_names_11=np.array(['V03','V05','V07','V09','V11','V13','V15','V17','V19','V21','V23','V25'])
    wind_speeds_11 = np.arange(3,26,2)

    # filesets for DLC 2.1
    seed_names_21=np.array(['S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12'])

    # filesets for DLC 2.1
    wind_speeds_names_21=np.array(['Vr','Vout'])
    wind_speeds_21=np.array([Vr,V_out])
    
    # filesets for EOG DLCs
    T_EOG = 10.5
    delta_EOG = np.array([0,T_EOG * 1/5,T_EOG * 2/5,T_EOG * 3/5,T_EOG * 4/5, T_EOG])
    delta_EOG_names = np.array([str(round(0,1)),str(round(T_EOG * 1/5,1)),str(round(T_EOG * 2/5,1)),str(round(T_EOG * 3/5,1)),str(round(T_EOG * 4/5,1)), str(round(T_EOG,1))])

    # filesets for EDC DLCs
    T_EDC = 6
    delta_EDC = np.array([0,T_EDC * 1/3,T_EDC * 2/3, T_EDC])
    delta_EDC_names = np.array([str(round(0,1)),str(round(T_EDC * 1/3,1)),str(round(T_EDC * 2/3,1)), str(round(T_EDC,1))])

    # filesets for DLC 1.3, DLC 5.1
    seed_names=np.array(['S1','S2','S3','S4','S5','S6'])
    nseeds_long = 6
    seed_namesl = []
    for sn in range(nseeds_long):
        seed_namesl.append('S'+str(sn+1))
    
    # filesets for DLC 1.3, DLC 1.5
    wind_speeds_names=np.array(['V04','V06','V08','V10','V11.4','V12','V14','V16','V18','V20','V22','V24'])
    wind_speeds=np.array(['4.0','6.0','8.0','10.0','11.4','12.0','14.0','16.0','18.0','20.0','22.0','24.0'])
    
    	# filesets for DLC 1.1, DLC 1.3, DLC 2.1, DLC 1.5, DLC 5.1
    directions_names=np.array(['+00deg'])
    directions=np.array(['0'])
    
    # filesets for DLC 1.4, DLC 1.5
    signs=np.array(['+','-'])
    wind_speeds_names_14=np.array(['Vr-2','Vr','Vr+2'])
    wind_speeds_14=np.array([Vr-2,Vr,Vr+2])
    DLC15_spec=np.array(['V','H'])

    # fileset for DLC 2.1, DLC 2.3, DLC 3.2, DLC 3.3, DLC 4.1, DLC 4.2
    pitch_angles = np.array([[0,11.4,12,13,14,15,16,17,18,19,20,21,22,23,24,25],[0,0,3.83,6.6,8.7,10.45,12.06,13.54,14.92,16.23,17.47,18.70,19.94,21.18,22.35,23.47]]) # initial pitch angles for wind speeds 11.4, 13, 15, 17, 19, 21, 23, 25
    pitch_angles_table=interp1d(pitch_angles[0],pitch_angles[1],kind='linear')

    # filesets for DLC 3.1, DLC 4.1
    wind_speeds_names_31_41=np.array(['V03','V11.4','V25'])
    wind_speeds_31_41=np.array([3,11.4,25])
    
    # filesets for DLC 2.3, DLC 3.2, DLC 3.3, DLC 4.2, DLC 5.1
    wind_speeds_names_2345=np.array(['Vr-2','Vr','Vr+2','Vout'])
    wind_speeds_2345=np.array([Vr-2,Vr,Vr+2,V_out])
    
    # filesets for DLC 6.1
    directions_names_61=np.array(['-08deg','+00deg','+08deg'])
    directions_61=np.array(['-8','0','8'])
    
    # filesets for DLC 6.2
    #directions_names_62=np.array(['-160deg','-140deg','-120deg','-100deg','-80deg','-60deg','-50deg','-40deg','-30deg','-20deg','-10deg','+00deg','+10deg','+20deg','+30deg','+40deg','+50deg','+60deg','+80deg','+100deg','+120deg','+140deg','+160deg','+180deg'])
    #directions_62=np.array(['-160','-140','-120','-100','-80','-60','-50','-40','-30','-20','-10','0','10','20','30','40','50','60','80','100','120','140','160','180'])
    directions_names_62=np.array(['-160deg','-140deg','-120deg','-100deg','-80deg','-60deg','+00deg','+60deg','+80deg','+100deg','+120deg','+140deg','+160deg','+180deg'])
    directions_62=np.array(['-160','-140','-120','-100','-80','-60','0','60','80','100','120','140','160','180'])
   

    # filesets for DLC 6.3
    directions_names_63=np.array(['+00deg'])
    directions_63=np.array(['0'])

    # filesets for DLC 6.4
    wind_speeds_names_64=np.array(['V02','V04','V06','V08','V10','V12','V14','V16','V18','V20','V22','V24','V26','V28','V30','V32','V34','V35'])
    wind_speeds_64=np.array(['2.0','4.0','6.0','8.0','10.0','12.0','14.0','16.0','18.0','20.0','22.0','24.0','26.0','28.0','30.0','32.0','34.0','35.0'])
    directions_names_64=np.array(['+00deg'])
    directions_64=np.array(['0'])

    # filesets for DLC 7.1
    directions_names_71=np.array(['-08deg','+00deg','+08deg'])
    directions_71=np.array(['-8','0','8'])

    # calculate number of overall files
    file_number_11 = len(seed_names_11)*len(wind_speeds_names_11)*len(directions)
    file_number_13 = 0
    file_number_14 = 0
    file_number_15 = 0
    file_number_21 = 0
    file_number_23 = 0
    file_number_31 = 0
    file_number_32 = 0
    file_number_33 = 0
    file_number_41 = 0
    file_number_42 = 0
    file_number_51 = 0
    file_number_61 = 0
    file_number_62 = 0
    file_number_63 = 0
    file_number_64 = 0
    file_number_71 = 0

    file_number_steady = file_number_31 + file_number_41
    file_number_IECwind = file_number_14 + file_number_15 + file_number_23 + file_number_32 + file_number_33 + file_number_42
    file_number_TurbSim = file_number_11 + file_number_13 + file_number_21 + file_number_51 + file_number_61 + file_number_62 + file_number_63 + file_number_64 + file_number_71
    
    file_number=file_number_TurbSim+file_number_IECwind+file_number_steady
    
    file_list=np.empty(file_number,dtype=list)
    file_list_raw=np.empty(file_number,dtype=list)
    
    
    # get V_ref and I_ref from IEC 61400-1
    if turb_intens == 'A':
        I_ref = 0.16
    if turb_intens == 'B':
        I_ref = 0.14
    if turb_intens == 'C':
        I_ref = 0.12
    if wind_class == 1:
        V_ref = 50
    if wind_class == 2:
        V_ref = 42.5
    if wind_class == 3:
        V_ref = 37.5
    
    Vm10_EWM50 = round(1.4 * V_ref * 1**pwr_law_exp_idling,3) # mean wind speed for 10 min time series of recurrence period 50 years
    Vm10_EWM1 = round(0.8 * Vm10_EWM50,3) # mean wind speed for 10 min time series of recurrence period 50 years
    
    # initialize file data
    commands_all_batch_FAST=["\n" for x in range(cores)]
    commands_all_batch_IEC=["\n" for x in range(cores)]
    commands_all_batch_TS=["\n" for x in range(cores)]
    
    ### generate batch files for execution of FAST and TurbSim ###
    for i in range(len(commands_all_batch_FAST)):
        string="start /B Call "+"_Run_FAST_"+str(i+1)+filename_batch_ending+"\n"
        commands_all_batch_FAST[i]=string
    for i in range(len(commands_all_batch_TS)):
        string="start /B Call "+"_Run_TurbSim_"+str(i+1)+filename_batch_ending+"\n"
        commands_all_batch_TS[i]=string
    for i in range(len(commands_all_batch_IEC)):
        string="start /B Call "+"_Run_Core"+str(i)+filename_batch_ending+"\n"
        commands_all_batch_IEC[i]=string
        
    filename_root_batch_FAST=root_FAST+filename_batch_all_FAST+filename_batch_ending
    with open(filename_root_batch_FAST, 'w') as file:
            file.writelines( commands_all_batch_FAST )
            file.close
    filename_root_batch_TS=root_TurbSim+filename_batch_all_TS+filename_batch_ending
    with open(filename_root_batch_TS, 'w') as file:
            file.writelines( commands_all_batch_TS )
            file.close
    filename_root_batch_IEC=root_IECwind+filename_batch_all_IEC+filename_batch_ending
    with open(filename_root_batch_IEC, 'w') as file:
            file.writelines( commands_all_batch_IEC )
            file.close

    # calculate distribution of jobs
    job_distr=int(np.ceil(file_number/cores))
    job_distr_rest = job_distr - file_number/cores
    job_distr_overflow = (file_number/cores)%1
    FAST_overflow_depot = 0
    job_distr_TurbSim=int(np.ceil(file_number_TurbSim/cores))
    job_distr_TurbSim_rest = job_distr_TurbSim - file_number_TurbSim/cores
    job_distr_TurbSim_overflow = (file_number_TurbSim/cores)%1
    TurbSim_overflow_depot = 0
    job_distr_IEC=int(np.ceil(file_number_IECwind/cores))
    
    # generate initialized batch files for distributed jobs
    commands_each_batch_FAST=["\n" for x in range(job_distr+2)]
    commands_each_batch_IEC=["\n" for x in range(2)]
    
    # create folders and batch files for each core to generate IEC wind files
    for i in range(cores):
        # generate folders
        dirname_each_batch_IEC = "Core" + str(i) + "\\"
        dirIEC=root_IECwind+dirname_each_batch_IEC
        if os.path.exists(dirIEC):
            shutil.rmtree(dirIEC)
        if not os.path.exists(dirIEC):
            os.makedirs(dirIEC)
        # copy initial IEC input file to each folder
        shutil.copyfile(root_IECwind+IEC_Input+IEC_ending, dirIEC+IEC_Input+IEC_ending)
        
        # generate batch file for each core
        filename_each_batch_IEC = "_Run_Core"
        filename_root_batch_IEC = root_IECwind + filename_each_batch_IEC + str(i) + filename_batch_ending
        commands_each_batch_IEC[0] = "D:\ncd "+root_IECwind+"Core" + str(i) + "\n"
        commands_each_batch_IEC[1] = "Call IECwind"
        with open(filename_root_batch_IEC, 'w') as file:
                file.writelines( commands_each_batch_IEC )
                file.close
    
    # read the default IEC input file
    filename_IEC=root_IECwind + IEC_Input + IEC_ending
    with open(filename_IEC, 'r') as file:
        data_IEC_raw = file.readlines()
    
    # initialize data matrix for IEC-windfiles
    data_IEC = np.empty([74,file_number_IECwind+cores], dtype = list)
    name_IEC = np.empty([file_number_IECwind], dtype = list)
    z_IEC = 0
    print "Generate "+str(file_number)+ " input files..."    
    z=0
    # DLC 1.1
    for j in range(len(wind_speeds_names_11)):
        seed_names_tmp = seed_names_11
        for i in range(len(seed_names_tmp)):
            for m in range(len(directions)):
                shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename11+"T_"+directions_names[m]+"_"+wind_speeds_names_11[j]+"_"+seed_names_tmp[i]+Input_ending)
                print filename11+"T_"+directions_names[m]+"_"+wind_speeds_names_11[j]+"_"+seed_names_tmp[i]+" was generated."
                filename11i=root_TurbSim+filename11+"T_"+directions_names[m]+"_"+wind_speeds_names_11[j]+"_"+seed_names_tmp[i]+Input_ending
                with open(filename11i, 'r') as file:
                    data = file.readlines()
                data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
                data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
                data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
                data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
                data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
                data[21] = '   '+str(AnalysisTime_TurbSim)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
                data[22] = '   '+str(AnalysisTime_TurbSim)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
                data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
                data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
                data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
                data[27] = '  '+str(directions[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
                data[33]= '"'+str(turb_intens)+'"            IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
                data[34]= '"NTM"          IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
                data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
                data[39] = str(wind_speeds_11[j])+'           URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
                data[41] = str(pwr_law_exp)+'          PLExp           - Power law exponent [-] (or "default")\n'
                with open(filename11i, 'w') as file:
                    file.writelines( data )
                    file.close
                file_list[z]=filename11+"T_"+directions_names[m]+"_"+wind_speeds_names_11[j]+"_"+seed_names_tmp[i]+Input_ending
                file_list_raw[z]=filename11+"T_"+directions_names[m]+"_"+wind_speeds_names_11[j]+"_"+seed_names_tmp[i]
                z=z+1
    """
    # DLC 1.3
    for j in range(len(wind_speeds_names)):
        seed_names_tmp = seed_names
        for i in range(len(seed_names_tmp)):    
            for m in range(len(directions)):
                shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename13+"T_"+directions_names[m]+"_"+wind_speeds_names[j]+"_"+seed_namesl[i]+Input_ending)
                print filename13+"T_"+directions_names[m]+"_"+wind_speeds_names[j]+"_"+seed_namesl[i]+" was generated."
                filename13i=root_TurbSim+filename13+"T_"+directions_names[m]+"_"+wind_speeds_names[j]+"_"+seed_namesl[i]+Input_ending
                with open(filename13i, 'r') as file:
                    data = file.readlines()
                data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
                data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
                data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
                data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
                data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
                data[21] = '   '+str(AnalysisTime_TurbSim)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
                data[22] = '   '+str(AnalysisTime_TurbSim)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
                data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
                data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
                data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
                data[27] = '  '+str(directions[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
                data[33]= '"'+str(turb_intens)+'"            IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
                data[34]= '"'+str(wind_class)+'ETM"         IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
                data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
                data[39] = wind_speeds[j]+'           URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
                data[41] = str(pwr_law_exp)+'          PLExp           - Power law exponent [-] (or "default")\n'
                with open(filename13i, 'w') as file:
                    file.writelines( data )
                    file.close
                file_list[z]=filename13+"T_"+directions_names[m]+"_"+wind_speeds_names[j]+"_"+seed_namesl[i]+Input_ending
                file_list_raw[z]=filename13+"T_"+directions_names[m]+"_"+wind_speeds_names[j]+"_"+seed_namesl[i]
                z=z+1
    
    # DLC 1.4
    data = data_IEC_raw
    for i in range(len(directions)):
        for j in range(len(wind_speeds_names_14)):
            for k in range(len(signs)):
                data = data_IEC_raw
                data[0]= 'ECD        CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
                data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
                data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
                data[4] = str(wind_speeds_14[j])+'       Vrated, m/s or ft/s\n'
                data[5] = '0          Vout, m/s or ft/s\n'
                data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
                data[7]= str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
                data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
                data[9]= str(signs[k])+'       SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
                data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
                data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
                data_IEC[0:, z_IEC] = data
                name_IEC[z_IEC] = filename14+"_"+directions_names[i]+"_"+signs[k]+"_"+wind_speeds_names_14[j]
                file_list[z]=filename14+"_"+directions_names[i]+"_"+signs[k]+"_"+wind_speeds_names_14[j]+IEC_ending
                file_list_raw[z]=filename14+"_"+directions_names[i]+"_"+signs[k]+"_"+wind_speeds_names_14[j]
                print filename14+"_"+directions_names[i]+"_"+signs[k]+"_"+wind_speeds_names_14[j]+" was generated."
                z=z+1
                z_IEC = z_IEC + 1
    
    
    # DLC 1.5
    for i in range(len(wind_speeds)):
        for j in range(len(DLC15_spec)):
            for k in range(len(signs)):
                for m in range(len(directions)):
                    data = data_IEC_raw
                    data[0]= 'EWS'+str(DLC15_spec[j])+'       CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
                    data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
                    data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
                    data[4] = str(wind_speeds[i])+'       Vrated, m/s or ft/s\n'
                    data[5] = '0          Vout, m/s or ft/s\n'
                    data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
                    data[7]=  str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
                    data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
                    data[9]= str(signs[k])+'          SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
                    data[10]= str(hubHt)+'        Wind turbine hub-height, m or ft\n'
                    data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
                    data_IEC[0:, z_IEC] = data
                    name_IEC[z_IEC] = filename15+"_"+directions_names[m]+"_"+DLC15_spec[j]+"_"+signs[k]+"_"+wind_speeds_names[i]
                    file_list[z]=filename15+"_"+directions_names[m]+"_"+DLC15_spec[j]+"_"+signs[k]+"_"+wind_speeds_names[i]+IEC_ending
                    file_list_raw[z]=filename15+"_"+directions_names[m]+"_"+DLC15_spec[j]+"_"+signs[k]+"_"+wind_speeds_names[i]
                    print filename15+"_"+directions_names[m]+"_"+DLC15_spec[j]+"_"+signs[k]+"_"+wind_speeds_names[i]+" was generated."
                    z=z+1
                    z_IEC = z_IEC + 1
    
    # DLC 2.1
    for i in range(len(seed_names_21)):
        for j in range(len(wind_speeds_names_21)):
            for m in range(len(directions)):
                shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename21+"T_"+directions_names[m]+"_"+wind_speeds_names_21[j]+"_"+seed_names_21[i]+Input_ending)
                print filename21+"T_"+directions_names[m]+"_"+wind_speeds_names_21[j]+"_"+seed_names_21[i]+" was generated."
                filename21i=root_TurbSim+filename21+"T_"+directions_names[m]+"_"+wind_speeds_names_21[j]+"_"+seed_names_21[i]+Input_ending
                with open(filename21i, 'r') as file:
                    data = file.readlines()
                data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
                data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
                data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
                data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
                data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
                data[21] = '   '+str(AnalysisTime_TurbSim)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
                data[22] = '   '+str(AnalysisTime_TurbSim)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
                data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
                data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
                data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
                data[27] = '  '+str(directions[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
                data[33]= '"'+str(turb_intens)+'"            IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
                data[34]= '"NTM"          IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
                data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
                data[39] = str(wind_speeds_21[j])+'           URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
                data[41] = str(pwr_law_exp)+'          PLExp           - Power law exponent [-] (or "default")\n'
                with open(filename21i, 'w') as file:
                    file.writelines( data )
                    file.close
                file_list[z]=filename21+"T_"+directions_names[m]+"_"+wind_speeds_names_21[j]+"_"+seed_names_21[i]+Input_ending
                file_list_raw[z]=filename21+"T_"+directions_names[m]+"_"+wind_speeds_names_21[j]+"_"+seed_names_21[i]
                z=z+1

    # DLC 2.3
    data = data_IEC_raw
    for i in range(len(delta_EOG)):
        for j in range(len(wind_speeds_names_2345)):
            data = data_IEC_raw
            data[0]= 'EOG50      CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
            data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
            data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
            data[4] = str(wind_speeds_2345[j])+'       Vrated, m/s or ft/s\n'
            data[5] = '0          Vout, m/s or ft/s\n'
            data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
            data[7]= str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
            data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
            data[9]= '+       SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
            data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
            data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
            data_IEC[0:, z_IEC] = data
            name_IEC[z_IEC] = filename23+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            file_list[z]=filename23+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+IEC_ending
            file_list_raw[z]=filename23+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            print filename23+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+" was generated."
            z=z+1
            z_IEC = z_IEC + 1

    # DLC 3.1
    for i in range(len(wind_speeds_names_31_41)):
        file_list[z]=filename31+"_"+wind_speeds_names_31_41[i]
        file_list_raw[z]=filename31+"_"+wind_speeds_names_31_41[i]
        z=z+1

    # DLC 3.2
    data = data_IEC_raw
    for i in range(len(delta_EOG)):
        for j in range(len(wind_speeds_names_2345)):
            data = data_IEC_raw
            data[0]= 'EOG50      CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
            data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
            data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
            data[4] = str(wind_speeds_2345[j])+'       Vrated, m/s or ft/s\n'
            data[5] = '0          Vout, m/s or ft/s\n'
            data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
            data[7]= str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
            data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
            data[9]= '+       SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
            data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
            data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
            data_IEC[0:, z_IEC] = data
            name_IEC[z_IEC] = filename32+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            file_list[z]=filename32+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+IEC_ending
            file_list_raw[z]=filename32+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            print filename32+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+" was generated."
            z=z+1
            z_IEC = z_IEC + 1
    
    # DLC 3.3
    data = data_IEC_raw
    for i in range(len(delta_EDC)):
        for j in range(len(wind_speeds_names_2345)):
            for k in range(len(signs)):
                data = data_IEC_raw
                data[0]= 'EDC50      CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
                data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
                data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
                data[4] = str(wind_speeds_2345[j])+'       Vrated, m/s or ft/s\n'
                data[5] = '0          Vout, m/s or ft/s\n'
                data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
                data[7]= str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
                data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
                data[9]= str(signs[k])+'       SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
                data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
                data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
                data_IEC[0:, z_IEC] = data
                name_IEC[z_IEC] = filename33+"_"+delta_EDC_names[i]+"_"+signs[k]+"_"+wind_speeds_names_2345[j]
                file_list[z]=filename33+"_"+delta_EDC_names[i]+"_"+signs[k]+"_"+wind_speeds_names_2345[j]+IEC_ending
                file_list_raw[z]=filename33+"_"+delta_EDC_names[i]+"_"+signs[k]+"_"+wind_speeds_names_2345[j]
                print filename33+"_"+delta_EDC_names[i]+"_"+signs[k]+"_"+wind_speeds_names_2345[j]+" was generated."
                z=z+1
                z_IEC = z_IEC + 1
    
    # DLC 4.1
    for i in range(len(wind_speeds_names_31_41)):
        file_list[z]=filename41+"_"+wind_speeds_names_31_41[i]
        file_list_raw[z]=filename41+"_"+wind_speeds_names_31_41[i]
        z=z+1

    # DLC 4.2
    for i in range(len(delta_EOG)):
        for j in range(len(wind_speeds_names_2345)):
            data = data_IEC_raw
            data[0]= 'EOG50      CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
            data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
            data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
            data[4] = str(wind_speeds_2345[j])+'       Vrated, m/s or ft/s\n'
            data[5] = '0          Vout, m/s or ft/s\n'
            data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
            data[7]= str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
            data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
            data[9]= '+       SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
            data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
            data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
            data_IEC[0:, z_IEC] = data
            name_IEC[z_IEC] = filename42+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            file_list[z]=filename42+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+IEC_ending
            file_list_raw[z]=filename42+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]
            print filename42+"_"+delta_EOG_names[i]+"_"+wind_speeds_names_2345[j]+" was generated."
            z=z+1
            z_IEC = z_IEC + 1
    
    # DLC 5.1
    for i in range(len(wind_speeds_names_2345)):
        for j in range(len(seed_names)):
            for k in range(len(directions_names)):
                shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename51+"T_"+directions_names[k]+"_"+wind_speeds_names_2345[i]+"_"+seed_names[j]+Input_ending)
                print filename51+"T_"+directions_names[k]+"_"+wind_speeds_names_2345[i]+"_"+seed_names[j]+" was generated."
                filename51i=root_TurbSim+filename51+"T_"+directions_names[k]+"_"+wind_speeds_names_2345[i]+"_"+seed_names[j]+Input_ending
                with open(filename51i, 'r') as file:
                    data = file.readlines()
                data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
                data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
                data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
                data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
                data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
                data[21] = '   '+str(TMax_FAST_trans)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
                data[22] = '   '+str(TMax_FAST_trans)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
                data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
                data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
                data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
                data[27] = '  '+str(directions[k])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
                data[33]= '"'+str(turb_intens)+'"            IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
                data[34]= '"NTM"          IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
                data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
                data[39] = str(wind_speeds_2345[i])+'           URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
                data[41] = str(pwr_law_exp)+'          PLExp           - Power law exponent [-] (or "default")\n'                
                with open(filename51i, 'w') as file:
                    file.writelines( data )
                    file.close
                file_list[z]=filename51+"T_"+directions_names[k]+"_"+wind_speeds_names_2345[i]+"_"+seed_names[j]+Input_ending
                file_list_raw[z]=filename51+"T_"+directions_names[k]+"_"+wind_speeds_names_2345[i]+"_"+seed_names[j]
                z=z+1
    """
    """
    # DLC 6.1
    data = data_IEC_raw
    for i in range(len(directions_61)): 
        data = data_IEC_raw
        data[0]= 'EWM50     CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
        data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
        data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
        data[4] = str(Vr)+'       Vrated, m/s or ft/s\n'
        data[5] = str(V_out)+'          Vout, m/s or ft/s\n'
        data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
        data[7]=  str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
        data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
        data[9]= '+          SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
        data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
        data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
        data_IEC[0:, z_IEC] = data
        name_IEC[z_IEC] = filename61+"_"+directions_names_61[i]
        file_list[z]=filename61+"_"+directions_names_61[i]+IEC_ending
        file_list_raw[z]=filename61+"_"+directions_names_61[i]
        print filename61+"_"+directions_names_61[i]+" was generated."
        z=z+1
        z_IEC = z_IEC + 1
    
    # DLC 6.2
    data = data_IEC_raw
    for i in range(len(directions_62)): 
        data = data_IEC_raw
        data[0]= 'EWM50     CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
        data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
        data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
        data[4] = str(Vr)+'       Vrated, m/s or ft/s\n'
        data[5] = str(V_out)+'          Vout, m/s or ft/s\n'
        data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
        data[7]=  str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
        data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
        data[9]= '+          SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
        data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
        data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
        data_IEC[0:, z_IEC] = data
        name_IEC[z_IEC] = filename62+"_"+directions_names_62[i]
        file_list[z]=filename62+"_"+directions_names_62[i]+IEC_ending
        file_list_raw[z]=filename62+"_"+directions_names_62[i]
        print filename62+"_"+directions_names_62[i]+" was generated."
        z=z+1
        z_IEC = z_IEC + 1
    
    # DLC 6.3
    data = data_IEC_raw
    for i in range(len(directions_63)): 
        data = data_IEC_raw
        data[0]= 'EWM01     CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
        data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
        data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
        data[4] = str(Vr)+'       Vrated, m/s or ft/s\n'
        data[5] = str(V_out)+'          Vout, m/s or ft/s\n'
        data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
        data[7]=  str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
        data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
        data[9]= '+          SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
        data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
        data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
        data_IEC[0:, z_IEC] = data
        name_IEC[z_IEC] = filename63+"_"+directions_names_63[i]
        file_list[z]=filename63+"_"+directions_names_63[i]+IEC_ending
        file_list_raw[z]=filename63+"_"+directions_names_63[i]
        print filename63+"_"+directions_names_63[i]+" was generated."
        z=z+1
        z_IEC = z_IEC + 1    

    # DLC 7.1
    data = data_IEC_raw
    for i in range(len(directions_71)): 
        data = data_IEC_raw
        data[0]= 'EWM01     CONDITION FOR WHICH WIND WILL BE GENERATED,(see Table 1 below)\n'
        data[2] = str(wind_class)+'          IEC WIND TURBINE CLASS\n'
        data[3] = str(turb_intens)+'          WIND TURBULENCE CATEGORY (A or B)\n'
        data[4] = str(Vr)+'       Vrated, m/s or ft/s\n'
        data[5] = str(V_out)+'          Vout, m/s or ft/s\n'
        data[6] = '0         Slope of the wind inflow (IEC specifies between -8 and +8), deg\n'
        data[7]=  str(IEC_trans_cond)+'        START OF IEC TRANSIENT CONDITION, sec\n'
        data[8]= str(DT_IECwind)+'       WIND FILE TIME STEP FOR TRANSIENT PORTION, sec\n'
        data[9]= '+          SIGN FOR WIND DIRECTION AND HORIZONTAL SHEAR (+ OR - OR "BOTH")\n'
        data[10]= str(hubHt)+'         Wind turbine hub-height, m or ft\n'
        data[11]= str(D_rot)+'        Wind turbine rotor diameter, m or ft\n'
        data_IEC[0:, z_IEC] = data
        name_IEC[z_IEC] = filename71+"_"+directions_names_71[i]
        file_list[z]=filename71+"_"+directions_names_71[i]+IEC_ending
        file_list_raw[z]=filename71+"_"+directions_names_71[i]
        print filename71+"_"+directions_names_71[i]+" was generated."
        z=z+1
        z_IEC = z_IEC + 1
    """
    """
    # DLC 6.1
    for i in range(len(seed_names)):
        for m in range(len(directions_61)):
            shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename61+"T_"+directions_names_61[m]+"_"+seed_names[i]+Input_ending)
            print filename61+"T_"+directions_names_61[m]+"_"+seed_names[i]+" was generated."
            filename61i=root_TurbSim+filename61+"T_"+directions_names_61[m]+"_"+seed_names[i]+Input_ending
            with open(filename61i, 'r') as file:
                data = file.readlines()
            data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
            data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
            data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
            data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
            data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
            data[21] = '   '+str(AnalysisTime_TurbSim_idling)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
            data[22] = '   '+str(AnalysisTime_TurbSim_idling)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
            data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
            data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
            data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
            data[27] = '  '+str(directions_61[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
            data[33]= '"'+turb_intens+'"             IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
            data[34]= '"'+str(wind_class)+'EWM50"         IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
            data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
            data[39] = str(Vm10_EWM50) + '          URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
            data[41] = str(pwr_law_exp_idling)+'          PLExp           - Power law exponent [-] (or "default")\n'
            with open(filename61i, 'w') as file:
                file.writelines( data )
                file.close
            file_list[z]=filename61+"T_"+directions_names_61[m]+"_"+seed_names[i]+Input_ending
            file_list_raw[z]=filename61+"T_"+directions_names_61[m]+"_"+seed_names[i]
            z=z+1

    # DLC 6.2
    for i in range(len(seed_names)):
        for m in range(len(directions_62)):
            shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename62+"T_"+directions_names_62[m]+"_"+seed_names[i]+Input_ending)
            print filename62+"T_"+directions_names_62[m]+"_"+seed_names[i]+" was generated."
            filename62i=root_TurbSim+filename62+"T_"+directions_names_62[m]+"_"+seed_names[i]+Input_ending
            with open(filename62i, 'r') as file:
                data = file.readlines()
            data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
            data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
            data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
            data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
            data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
            data[21] = '   '+str(AnalysisTime_TurbSim_idling)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
            data[22] = '   '+str(AnalysisTime_TurbSim_idling)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
            data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
            data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
            data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
            data[27] = '  '+str(directions_62[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
            data[33]= '"'+turb_intens+'"             IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
            data[34]= '"'+str(wind_class)+'EWM50"         IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
            data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
            data[39] = str(Vm10_EWM50) + '          URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
            data[41] = str(pwr_law_exp_idling)+'          PLExp           - Power law exponent [-] (or "default")\n'
            with open(filename62i, 'w') as file:
                file.writelines( data )
                file.close
            file_list[z]=filename62+"T_"+directions_names_62[m]+"_"+seed_names[i]+Input_ending
            file_list_raw[z]=filename62+"T_"+directions_names_62[m]+"_"+seed_names[i]
            z=z+1
    
    # DLC 6.3
    for i in range(len(seed_names)):
        for m in range(len(directions_63)):
            shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename63+"T_"+directions_names_63[m]+"_"+seed_names[i]+Input_ending)
            print filename63+"T_"+directions_names_63[m]+"_"+seed_names[i]+" was generated."
            filename63i=root_TurbSim+filename63+"T_"+directions_names_63[m]+"_"+seed_names[i]+Input_ending
            with open(filename63i, 'r') as file:
                data = file.readlines()
            data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
            data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
            data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
            data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
            data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
            data[21] = '   '+str(AnalysisTime_TurbSim_idling)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
            data[22] = '   '+str(AnalysisTime_TurbSim_idling)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
            data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
            data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
            data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
            data[27] = '  '+str(directions_63[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
            data[33]= '"'+turb_intens+'"             IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
            data[34]= '"'+str(wind_class)+'EWM1"         IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
            data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
            data[39] = str(Vm10_EWM1) + '          URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
            data[41] = str(pwr_law_exp_idling)+'          PLExp           - Power law exponent [-] (or "default")\n'
            with open(filename63i, 'w') as file:
                file.writelines( data )
                file.close
            file_list[z]=filename63+"T_"+directions_names_63[m]+"_"+seed_names[i]+Input_ending
            file_list_raw[z]=filename63+"T_"+directions_names_63[m]+"_"+seed_names[i]
            z=z+1

    # DLC 6.4
    for i in range(len(wind_speeds_names_64)):
        for j in range(len(seed_names)):
            for k in range(len(directions_names_64)):
                shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename64+"T_"+directions_names_64[k]+"_"+wind_speeds_names_64[i]+"_"+seed_names[j]+Input_ending)
                print filename64+"T_"+directions_names_64[k]+"_"+wind_speeds_names_64[i]+"_"+seed_names[j]+" was generated."
                filename64i=root_TurbSim+filename64+"T_"+directions_names_64[k]+"_"+wind_speeds_names_64[i]+"_"+seed_names[j]+Input_ending
                with open(filename64i, 'r') as file:
                    data = file.readlines()
                data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
                data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
                data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
                data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
                data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
                data[21] = '   '+str(AnalysisTime_TurbSim_idling)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
                data[22] = '   '+str(AnalysisTime_TurbSim_idling)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
                data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
                data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
                data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
                data[27] = '  0            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
                data[33]= '"'+str(turb_intens)+'"            IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
                data[34]= '"NTM"          IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
                data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
                data[39] = str(wind_speeds_64[i])+'           URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
                data[41] = str(pwr_law_exp)+'          PLExp           - Power law exponent [-] (or "default")\n'  
                with open(filename64i, 'w') as file:
                    file.writelines( data )
                    file.close
                file_list[z]=filename64+"T_"+directions_names_64[k]+"_"+wind_speeds_names_64[i]+"_"+seed_names[j]+Input_ending
                file_list_raw[z]=filename64+"T_"+directions_names_64[k]+"_"+wind_speeds_names_64[i]+"_"+seed_names[j]
                z=z+1
    
    # DLC 7.1
    for i in range(len(seed_names)):
        for m in range(len(directions_71)):
            shutil.copyfile(root_TurbSim+TurbSim_Input+Input_ending, root_TurbSim+filename71+"T_"+directions_names_71[m]+"_"+seed_names[i]+Input_ending)
            print filename71+"T_"+directions_names_71[m]+"_"+seed_names[i]+" was generated."
            filename71i=root_TurbSim+filename71+"T_"+directions_names_71[m]+"_"+seed_names[i]+Input_ending
            with open(filename71i, 'r') as file:
                data = file.readlines()
            data[4]= str(int(random()*2147483647))+'          RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'
            data[5]= str(int(random()*2147483647))+'          RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'
            data[18] = ' '+str(NumGrid_Z)+'          NumGrid_Z       - Vertical grid-point matrix dimension\n'
            data[19] = ' '+str(NumGrid_Y)+'          NumGrid_Y       - Horizontal grid-point matrix dimension\n'
            data[20] = ' '+str(DT_TurbSim)+'         TimeStep        - Time step [seconds]\n'
            data[21] = '   '+str(AnalysisTime_TurbSim_idling)+'         AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'
            data[22] = '   '+str(AnalysisTime_TurbSim_idling)+'         UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'
            data[23]= ' '+str(hubHt)+'           HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'
            data[24] = str(GridHeight)+'          GridHeight      - Grid height [m]\n'
            data[25] = str(GridWidth)+'          GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'
            data[27] = '  '+str(directions_71[m])+'            HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'
            data[33]= '"'+turb_intens+'"             IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)\n'
            data[34]= '"'+str(wind_class)+'EWM1"         IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'
            data[38] = ' '+str(hubHt)+'           RefHt           - Height of the reference wind speed [m]\n'
            data[39] = str(Vm10_EWM1) + '          URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'
            data[41] = str(pwr_law_exp_idling)+'          PLExp           - Power law exponent [-] (or "default")\n'
            with open(filename71i, 'w') as file:
                file.writelines( data )
                file.close
            file_list[z]=filename71+"T_"+directions_names_71[m]+"_"+seed_names[i]+Input_ending
            file_list_raw[z]=filename71+"T_"+directions_names_71[m]+"_"+seed_names[i]
            z=z+1
    """
    # generate batch files for each file
    z1=0
    z2=0
    for i in range(cores):
        # FAST
        filename_each_batch_FAST=root_FAST+"_Run_FAST_"+str(i+1)+filename_batch_ending
        commands_each_batch_FAST=["\n" for x in range(job_distr+2)]
        if FAST_overflow_depot <= 0:
            job_distr_current = job_distr
            FAST_overflow_depot = FAST_overflow_depot + job_distr_rest
        else:
            job_distr_current = job_distr - 1
            FAST_overflow_depot = FAST_overflow_depot - job_distr_overflow
        for j in range(job_distr_current):
            if z1<file_number:
                string="FAST_x64 "+str(file_list_raw[z1])+filename_FAST_end+"\n"
                commands_each_batch_FAST[j]=string
            if z1>=file_number:
                string="\n"
                commands_each_batch_FAST[j]=string
            z1=z1+1
        with open(filename_each_batch_FAST, 'w') as file:
                file.writelines( commands_each_batch_FAST )
                file.close
        # TurbSim
        filename_each_batch_TS=root_TurbSim+"_Run_TurbSim_"+str(i+1)+filename_batch_ending
        commands_each_batch_TS=["\n" for x in range(job_distr_TurbSim+2)]
        if TurbSim_overflow_depot <= 0:
            job_distr_TurbSim_current = job_distr_TurbSim
            TurbSim_overflow_depot = TurbSim_overflow_depot + job_distr_TurbSim_rest
        else:
            job_distr_TurbSim_current = job_distr_TurbSim - 1
            TurbSim_overflow_depot = TurbSim_overflow_depot - job_distr_TurbSim_overflow
        for k in range(z2,len(file_list_raw)):        
            actual_DLC_tmp = file_list_raw[k]
            actual_DLC_type = actual_DLC_tmp[6]
            print job_distr_TurbSim_current
            if actual_DLC_type != 'T':
                z2 = k + 1
                print actual_DLC_type, z2, k
            print "test"
            if job_distr_TurbSim_current > 0:
                if actual_DLC_type == 'T':
                    string="TurbSim_x64 "+file_list[k]+"\n"
                    print actual_DLC_type, string, k, z2, len(commands_each_batch_TS), k - z2, commands_each_batch_TS[k - z2]
                    commands_each_batch_TS[k - z2]=string                    
                    job_distr_TurbSim_current = job_distr_TurbSim_current - 1
            if job_distr_TurbSim_current <= 0:
                z2 = k + 1
                break
            #if k >= file_number:
                #string="\n"
                #commands_each_batch_TS[j]=string      
        with open(filename_each_batch_TS, 'w') as file:
                file.writelines( commands_each_batch_TS )
                file.close
        
    # generate inflow files for each job
    print "Generate "+str(file_number)+" inflow Wind files..."
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        filename_each_inflow=filename_inflow+"_"+file_list_raw[i]+filename_inflow_end
        shutil.copyfile(root_FAST+filename_inflow+filename_inflow_end, root_FAST+"\InflowWind\\"+filename_inflow+"_"+file_list_raw[i]+filename_inflow_end)
        print filename_each_inflow+" was generated."
        filename=root_FAST+"InflowWind\\"+filename_each_inflow
        with open(filename, 'r') as file:
            data = file.readlines()
        data[4]='          3   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n'
        data[9]='         '+str(hubHt)+'   WindVziList    - List of coordinates in the inertial Z direction (m)\n'
        data[12]='         '+str(hubHt)+'   RefHt          - Reference height for horizontal wind speed\n'
        data[19]='"Wind/'+file_list_raw[i]+'.bts"    Filename       - Name of the Full field wind file to use (.bts)\n'
        if actual_DLC == '14' or actual_DLC == '15' or actual_DLC == '23' or actual_DLC == '32' or actual_DLC == '33' or actual_DLC == '42':
            data[4]='          2   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n'
            data[15]='"Wind/'+file_list_raw[i]+'.WND"    Filename       - Filename of time series data for uniform wind field.\n'
            data[16]='    '+str(hubHt)+'   RefHt          - Reference height for horizontal wind speed\n'
            data[17]='     '+str(D_rot)+'   RefLength      - Reference length for linear horizontal and vertical sheer\n'
        if actual_DLC =='31':
            V_wind=float(actual_DLC_tmp[8:])
            data[4]='          1   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n'
            data[11]='        '+str(V_wind)+'   HWindSpeed     - Horizontal windspeed\n'
            data[12]='        '+str(hubHt)+'   RefHt          - Reference height for horizontal wind speed\n'
        if actual_DLC =='41':
            V_wind=float(actual_DLC_tmp[8:])
            data[4]='          1   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n'
            data[11]='        '+str(V_wind)+'   HWindSpeed     - Horizontal windspeed\n'
            data[12]='        '+str(hubHt)+'   RefHt          - Reference height for horizontal wind speed\n'
        with open(filename, 'w') as file:
            file.writelines( data )
            file.close
    print "All inflow Wind files were generated."

    # generate AeroDyn files for each job
    """
    print "Generate "+str(file_number)+" AeroDyn files..."
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        filename_each_aero=filename_aero+"_"+file_list_raw[i]+filename_aero_end
        shutil.copyfile(root_FAST+filename_aero+filename_aero_end, root_FAST+"\AeroDyn\\"+filename_aero+"_"+file_list_raw[i]+filename_aero_end)
        print filename_each_aero+" was generated."
        filename=root_FAST+"AeroDyn\\"+filename_each_aero
        with open(filename, 'r') as file:
            data = file.readlines()
        #if iteration>1:
            #data[12]='"AeroStructure/Ellipse/NREL_5MW_AeroDyn_Tower.dat"    TwrFile        - Tower drag file name (quoted string)\n'
        #data[12]='"AeroStructure/Conventional/NREL_5MW_AeroDyn_Tower.dat"    TwrFile        - Tower drag file name (quoted string)\n'
        data[4]='0.02479       DTAero             - Time interval for aerodynamic calculations {or "default"} (s)\n'
        data[7]='"          1   TwrPotent          - Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}\n'
        data[16]='          1   SkewMod            - Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]\n'
        data[17]='True          TipLoss            - Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]\n'
        data[18]='True          HubLoss            - Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]\n'
        data[19]='False         TanInd             - Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]\n'
        data[20]='True          AIDrag             - Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]\n'
        data[21]='False         TIDrag             - Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]\n'
        data[22]='"Default"     IndToler           - Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [used only when WakeMod=1]\n'
        data[23]='        100   MaxIter            - Maximum number of iteration steps (-) [used only when WakeMod=1]\n'
        data[25]='          1   UAMod              - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]\n'
        data[34]='"Airfoils/Cylinder1.dat"    FoilNm      - Names of the airfoil files [NumFoil lines] (quoted strings)\n'
        data[34]='"Airfoils/Cylinder1.dat"    FoilNm      - Names of the airfoil files [NumFoil lines] (quoted strings)\n'
        data[35]='"Airfoils/Cylinder2.dat"\n'
        data[36]='"Airfoils/DU40_A17.dat"\n'
        data[37]='"Airfoils/DU35_A17.dat"\n'
        data[38]='"Airfoils/DU30_A17.dat"\n'
        data[39]='"Airfoils/DU25_A17.dat"\n'
        data[40]='"Airfoils/DU21_A17.dat"\n'
        data[41]='"Airfoils/NACA64_A17.dat"\n'
        if actual_DLC[0] == "6" or actual_DLC[0] == "7":
            data[5]='          0   WakeMod            - Type of wake/induction model (switch) {0=none, 1=BEMT}\n'
            data[6]='          1   AFAeroMod          - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model}\n'
        # in the following line is a fault concerning the right position of the files of the file_list
        with open(filename, 'w') as file:
            file.writelines( data )
            file.close
    print "All AeroDyn files were generated."
   """
   
    # generate ElastoDyn files for each job
    print "Generate "+str(file_number)+" ElastoDyn files..."
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        filename_each_elasto=filename_elasto+"_"+file_list_raw[i]+filename_elasto_end
        shutil.copyfile(root_FAST+filename_elasto+filename_elasto_end, root_FAST+"\ElastoDyn\\"+filename_elasto+"_"+file_list_raw[i]+filename_elasto_end)
        print filename_each_elasto+" was generated."
        filename=root_FAST+"ElastoDyn\\"+filename_each_elasto
        with open(filename, 'r') as file:
            data = file.readlines()
        data[27]='          '+str(OoPDefl)+'   OoPDefl     - Initial out-of-plane blade-tip displacement (meters)\n'
        data[65]='       '+str(towerHt)+'   TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)\n'
        
        data[82] = str(inertia_nacelle_yaw / 2.0) +"   PtfmRIner   - Platform inertia for roll tilt rotation about the platform CM (kg m^2)\n"
        data[83] = str(inertia_nacelle_yaw) + "   PtfmPIner   - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)\n"
        data[87]='"Blade/NREL_5MW_Blade.dat"    BldFile(1)  - Name of file containing properties for blade 1 (quoted string)\n'
        data[88]='"Blade/NREL_5MW_Blade_light.dat"    BldFile(2)  - Name of file containing properties for blade 2 (quoted string)\n'
        data[89]='"Blade/NREL_5MW_Blade_heavy.dat"    BldFile(3)  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]\n'
        data[109]='"Tower/NREL_5MW_ElastoDyn_Tower.dat"    TwrFile     - Name of file containing tower properties (quoted string)\n'
        # for DLC 11 and 13
        if ((actual_DLC == '11') or (actual_DLC == '13')) and ((file_list_raw[i][16] == '1') or (file_list_raw[i][16] == '2')):
            tmp_str_vw = float(file_list_raw[i][16:18])
            if file_list_raw[i][17] == '1':
                tmp_str_vw = 11.4
            data[29]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
            data[30]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
            data[31]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for DLC 21
        if actual_DLC == '21':
            if actual_DLC_tmp[19] == '_' :
                tmp_str_vw = float(V_out)
            if actual_DLC_tmp[17] == '_' :
                tmp_str_vw = float(Vr)
            data[29]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
            data[30]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
            data[31]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for DLC 23
        if actual_DLC == '23':
            if actual_DLC_tmp[12] == 'r':
                tmp_str_vw = float(Vr)
            if actual_DLC_tmp[12] == 'o':
                tmp_str_vw = float(V_out)
            if len(actual_DLC_tmp) > 13:
                if actual_DLC_tmp[13] == 'r' :
                    tmp_str_vw = float(Vr)
                if actual_DLC_tmp[13] == 'o' :
                    tmp_str_vw = float(V_out)
            data[29]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
            data[30]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
            data[31]='     '+str(pitch_angles_table(tmp_str_vw))+'   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for DLC 31
        if actual_DLC == '31' or actual_DLC == '32' or actual_DLC == '33':  
            data[29]='         90   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
            data[30]='         90   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
            data[31]='         90   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            data[34]='          0   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for DLC 41 and DLC 42
        if actual_DLC == '41' or actual_DLC == '42':
            if (actual_DLC_tmp =='Dlc_41_V03') or (actual_DLC_tmp =='Dlc_41_V11.4') or ((actual_DLC == '42') and not (actual_DLC_tmp[-2:] == 'ut')) or ((actual_DLC == '42') and not (actual_DLC_tmp[-2:] == '+2')):
                data[29]='          0   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                data[30]='          0   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                data[31]='          0   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
                data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
            if (actual_DLC == '42') and (actual_DLC_tmp[-2:] == '+2'):
                data[29]='      '+str(pitch_angles_table(Vr+2))+'   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                data[30]='      '+str(pitch_angles_table(Vr+2))+'   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                data[31]='      '+str(pitch_angles_table(Vr+2))+'   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
                data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
            if (actual_DLC_tmp =='Dlc_41_V25') or ((actual_DLC == '42') and (actual_DLC_tmp[-2:] == 'ut')):
                data[29]='      23.47   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                data[30]='      23.47   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                data[31]='      23.47   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
                data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for DLC 51
        if actual_DLC == '51':
            if (actual_DLC_tmp =='Dlc_51_V03') or (actual_DLC_tmp =='Dlc_51_V11.4'):
                data[29]='          0   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                data[30]='          0   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                data[31]='          0   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
                data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
            if actual_DLC_tmp =='Dlc_51_V25':
                data[29]='      23.47   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                data[30]='      23.47   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                data[31]='      23.47   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
                data[34]='       12.1   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
        # for parking/idling DLCs 61, 62, 63, 64, and 71
        if actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '64' or actual_DLC == '71':
            data[14]='False         GenDOF      - Generator DOF (flag)\n' # for parking
            #data[14]='True          GenDOF      - Generator DOF (flag)\n' # for idling
            data[29]='         90   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
            data[30]='         90   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
            data[31]='         90   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            data[34]='          0   RotSpeed    - Initial or fixed rotor speed (rpm)\n'
            
            if actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '71':
                if actual_DLC_tmp[12] == 'd':
                    WindDir = float(actual_DLC_tmp[8:12])
                if actual_DLC_tmp[12] != 'd':
                    WindDir = float(actual_DLC_tmp[8:11])

                # calculate the blades pitch angle such, that the trailing edge of the top blade points downwind.
                BlPitch = 90 - WindDir
                if WindDir < -90:
                    BlPitch = -270 - WindDir
                data[29]='         '+str(BlPitch)+'   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'
                #data[30]='         '+str(BlPitch)+'   BlPitch(2)  - Blade 2 initial pitch (degrees)\n'
                #data[31]='         '+str(BlPitch)+'   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n'
            
            if actual_DLC == '71':
                data[29]='          0   BlPitch(1)  - Blade 1 initial pitch (degrees)\n'

            data[100]='         '+str(GBoxEff)+'   GBoxEff     - Gearbox efficiency (%)\n' # for idling
        data[111]='True          SumPrint    - Print summary data to "<RootName>.sum" (flag)\n'
        with open(filename, 'w') as file:
            file.writelines( data )
            file.close
    print "All ElastoDyn files were generated."
    
    # generate ServoDyn files for each job
    print "Generate "+str(file_number)+" ServoDyn files..."
    z23 = 0
    z23_ = 0
    z32 = 0
    z32_ = 0
    z33 = 0
    z33_ = 0    
    z42 = 0
    z42_ = 0
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        filename_each_servo=filename_servo+"_"+file_list_raw[i]+filename_servo_end
        shutil.copyfile(root_FAST+filename_servo+filename_servo_end, root_FAST+"\ServoDyn\\"+filename_servo+"_"+file_list_raw[i]+filename_servo_end)
        print filename_each_servo+" was generated."
        filename=root_FAST+"ServoDyn\\"+filename_each_servo
        with open(filename, 'r') as file:
            data = file.readlines()         
        # for DLC 21
        if actual_DLC == '21':
            data[6]='          5   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[22]='False         GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'
            data[24]='          0   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[25]='     9999.9   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]\n'
            data[8]='         '+str(IEC_trans_cond)+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond+0.2)+' TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond+0.2)+' TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='          0   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
        # for DLC 23
        if actual_DLC == '23':
            data[6]='          5   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[22]='True          GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'
            data[24]='          0   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[25]='        '+str(IEC_trans_cond + delta_EOG[z23])+'   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]\n'
            data[8]='         '+str(IEC_trans_cond + delta_EOG[z23] + 0.2)+' TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond + delta_EOG[z23] + 0.2)+' TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond + delta_EOG[z23] + 0.2)+' TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='         90   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            z23_ += 1
            if z23_ == len(wind_speeds_names_2345):
                z23_ = 0
                z23 += 1        
        # for DLC 31
        if actual_DLC == '31':
            data[6]='          0   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[21]='True         GenTiStr     - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n'
            data[23]='        '+str(SpdGenOn)+'   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n'
            data[24]='          '+str(0)+'   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[8]='         '+str(IEC_trans_cond)+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond)+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond)+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='          0   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='          0   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='          0   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            if actual_DLC_tmp == 'Dlc_31_V25':
                data[14]='      23.47   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
                data[15]='      23.47   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
                data[16]='      23.47   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
        # for DLC 32
        if (actual_DLC == '32'):
            data[6]='          0   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[21]='True         GenTiStr     - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n'
            data[23]='        '+str(SpdGenOn)+'   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n'
            data[24]='          '+str(0)+'   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[8]='         '+str(IEC_trans_cond + delta_EOG[z32])+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond + delta_EOG[z32])+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond + delta_EOG[z32])+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='          0   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='          0   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='          0   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            if actual_DLC_tmp[-2:] == 'ut':
                data[14]='      23.47   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
                data[15]='      23.47   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
                data[16]='      23.47   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            if actual_DLC_tmp[-2:] == '+2':
                data[14]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
                data[15]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
                data[16]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            z32_ += 1
            if z32_ == len(wind_speeds_names_2345):
                z32_ = 0
                z32 += 1
        # for DLC 33        
        if (actual_DLC == '33'):
            data[6]='          0   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[21]='True         GenTiStr     - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n'
            data[23]='        '+str(SpdGenOn)+'   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n'
            data[24]='          '+str(0)+'   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[8]='         '+str(IEC_trans_cond + delta_EDC[z33])+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond + delta_EDC[z33])+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond + delta_EDC[z33])+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='          0   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='          0   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='          0   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            if actual_DLC_tmp[-2:] == 'ut':
                data[14]='      '+str(20.99)+'   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
                data[15]='      '+str(20.99)+'   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
                data[16]='      '+str(20.99)+'   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            if actual_DLC_tmp[-2:] == '+2':
                data[14]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
                data[15]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
                data[16]='      '+str(pitch_angles_table(Vr+2))+'   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            z33_ += 1
            if z33_ == len(wind_speeds_names_2345)*2:
                z33_ = 0
                z33 += 1 
        # for DLC 41
        if actual_DLC == '41':
            data[6]='          5   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[22]='False         GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'
            data[24]='          0   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[8]='         '+str(IEC_trans_cond)+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond)+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond)+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='          '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='          '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='          '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='         90   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
        # for DLC 42
        if actual_DLC == '42':
            data[6]='          5   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[22]='False         GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'
            data[24]='          0   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[8]='         '+str(IEC_trans_cond + delta_EOG[z42])+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond + delta_EOG[z42])+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond + delta_EOG[z42])+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='          '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='          '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='          '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='         90   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            z42_ += 1
            if z42_ == len(wind_speeds_names_2345):
                z42_ = 0
                z42 += 1
        # for DLC 51
        if actual_DLC == '51':
            data[6]='          5   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[22]='False         GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'
            data[24]='          0   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'
            data[25]='     9999.9   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]\n'
            data[8]='         '+str(IEC_trans_cond)+'   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='         '+str(IEC_trans_cond)+'   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='         '+str(IEC_trans_cond)+'   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='         90   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            data[46]='          1   HSSBrMode    - HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[47]='         '+str(IEC_trans_cond)+'   THSSBrDp     - Time to initiate deployment of the HSS brake (s)\n'
        # for parking/idling DLCs
        if actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '64' or actual_DLC == '71':
            data[6]='          0   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'            
            data[8]='     9999.9   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'
            data[9]='     9999.9   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'
            data[10]='     9999.9   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'
            data[11]='        '+str(PitManRat)+'   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'
            data[12]='        '+str(PitManRat)+'   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'
            data[13]='        '+str(PitManRat)+'   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'
            data[14]='         90   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n'
            data[15]='         90   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n'
            data[16]='         90   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'
            data[24]='          '+str(TMax_idling*2)+'   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n' # for idling
            data[46]='          1   HSSBrMode    - HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'
            data[47]='          0   THSSBrDp     - Time to initiate deployment of the HSS brake (s)\n'
            data[49]='    '+str(28116.2)+'   HSSBrTqF     - Fully deployed HSS-brake torque (N-m)\n' # 10 percent brake torque for idling

        with open(filename, 'w') as file:
            file.writelines( data )
            file.close
    print "All ServoDyn files were generated."
    
    # generate FAST .fst files for each job
    print "Generate "+str(file_number)+" FAST inputs..."
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        filename_each_FAST=file_list_raw[i]+filename_FAST_end
        shutil.copyfile(root_FAST+'DLC_standard'+filename_FAST_end, root_FAST+filename_each_FAST)
        print filename_each_FAST+" was generated."
        filename=root_FAST+filename_each_FAST
        with open(filename, 'r') as file:
            data = file.readlines()
        data[5]='        '+str(TMax_FAST)+'   TMax            - Total run time (s)\n'
        # for DLC 14, DLC 15, DLC 23, DLC 61, DLC 63
        if actual_DLC == '14' or actual_DLC == '15' or actual_DLC == '23' or actual_DLC == '31' or actual_DLC == '32' or actual_DLC == '33' or actual_DLC == '41' or actual_DLC == '42' or actual_DLC == '51':
            data[5]='        '+str(TMax_FAST_trans)+'   TMax            - Total run time (s)\n'
        if actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '64' or actual_DLC == '71':
            data[5]='        '+str(TMax_idling)+'   TMax            - Total run time (s)\n'
        data[6]='       '+str(DT_FAST)+'   DT              - Recommended module time step (s)\n'
        data[14]='          1   CompAero        - Compute aerodynamic loads (switch) {0=None; 1=AeroDyn v14; 2=AeroDyn v15}\n'
        data[21]='"ElastoDyn/NREL_5MW_ElastoDyn_'+file_list_raw[i]+'.dat"    EDFile          - Name of file containing ElastoDyn input parameters (quoted string)\n'
        data[25]='"InflowWind/NREL_5MW_InflowWind_'+file_list_raw[i]+'.dat"    InflowFile      - Name of file containing inflow wind input parameters (quoted string)\n'
        data[26]='"NREL_5MW_AeroDyn.dat"    AeroFile        - Name of file containing aerodynamic input parameters (quoted string)\n'
        if actual_DLC[0] == "6" or actual_DLC[0] == "7":
            data[26]='"NREL_5MW_AeroDyn_idle.dat"    AeroFile        - Name of file containing aerodynamic input parameters (quoted string)\n'
        #data[26]='"AeroDyn/NREL_5MW_AeroDyn15_'+file_list_raw[i]+'.dat"    AeroFile        - Name of file containing aerodynamic input parameters (quoted string)\n'
        data[27]='"ServoDyn/NREL_5MW_ServoDyn_'+file_list_raw[i]+'.dat"    ServoFile       - Name of file containing control and electrical-drive input parameters (quoted string)\n'
        data[34]='          '+str(SttsTime)+'   SttsTime        - Amount of time between screen status messages (s)\n'
        data[36]='       '+str(DT_out)+'   DT_Out          - Time step for tabular output (s) (or "default")\n'
        data[37]='	  '+str(TStart)+'   TStart          - Time to begin tabular output (s)\n'
        if iteration>1:
            data[26]='"AeroDyn/NREL_5MW_AeroDyn_'+file_list_raw[i]+'.dat"    AeroFile        - Name of file containing aerodynamic input parameters (quoted string)\n'
        # for parking/idling DLCs
        if actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '64' or actual_DLC == '71':
            data[37]='	  '+str(TStart_parked)+'   TStart          - Time to begin tabular output (s)\n'
            
        with open(filename, 'w') as file:
            file.writelines( data )
            file.close
    print "All FAST .fst files were generated."



    # set filename of Crunch-input-file
    crunch_filename="Input_Crunch.cru"
        
    print "Manipulate Crunch input file..."
    # manipulate Crunch input-file by considering all copied Dlc_11 files        
    with open(crunch_filename, 'r') as file:
            data = file.readlines()
    
    # initialize new file data
    new_data=[""]*(104+file_number_11)
    new_data[0:99]=data[0:99    ]
    
    # constant entries
    new_data[99]=str(file_number_11)+'              NumFiles:           The number of input files to read.\n'
    
    z=100
    n_discevent=0
    for i in range(int(file_number)):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        if actual_DLC == '11':
            z=z+1
            new_data[z]='"'+"Out_Files/"+actual_DLC_tmp+".out"+'"\n'
    with open(crunch_filename, 'w') as file:
        file.writelines( new_data )
        file.close
    print "Crunch file was prepared."
    
    
    
    # prepare Mlife files
    print "Prepare files for MLife fatigue analysis."
    # set filename of Mlife-input-file
    mlife_filename="Input_MLife.mlif"
        
    print "Manipulate Mlife input file..."
    # manipulate Mlife input-file by considering all copied Dlc_11 files        
    with open(mlife_filename, 'r') as file:
            data = file.readlines()
    
    # initialize new file data
    new_data=[""]*(134+file_number_11 + file_number_64 + file_number_31 + file_number_41)
    new_data[0:130]=data[0:130]
    
    # constant entries
    new_data[128]=str(file_number_11)+'  1.1   1.3   1.5   1.7    (Weibull-Weighted Normal Operation: NumNormFiles, PSF1, PSF2, PSF3, PSF4)\n'
    new_data[131+file_number_11]=str(file_number_64)+'  1.1   1.3   1.5   1.7    (Weibull-Weighted Idling: NumIdleFiles, PSF1, PSF2, PSF3, PSF4)\n'
    new_data[132+file_number_11 + file_number_64]=str(int(file_number_31 + file_number_41))+"  1.2   1.3   1.4   1.6    (Discrete Events: NumDiscFiles, PSF1, PSF2, PSF3, PSF4)\n"
    new_data[133+file_number_11+file_number_64+file_number_31 + file_number_41]="==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.\n"
    
    z=128
    n_discevent=0
    for i in range(int(file_number)):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        if actual_DLC == '11':
            z=z+1
            new_data[z]='"'+"Out_Files/"+actual_DLC_tmp+".out"+'"\n'
        if actual_DLC == '64':
            z=z+1
            new_data[z]='"'+"Out_Files/"+actual_DLC_tmp+".out"+'"\n'
        if actual_DLC == '31' or actual_DLC == '41':
            z=z+1
            if actual_DLC_tmp == 'Dlc_31_V03':
                n_discevent=n_Vin31
            if actual_DLC_tmp == 'Dlc_31_V11.4':
                n_discevent=n_Vr31
            if actual_DLC_tmp == 'Dlc_31_V25':
                n_discevent=n_Vout31
            if actual_DLC_tmp == 'Dlc_41_V03':
                n_discevent=n_Vin41
            if actual_DLC_tmp == 'Dlc_41_V11.4':
                n_discevent=n_Vr41
            if actual_DLC_tmp == 'Dlc_41_V25':
                n_discevent=n_Vout41
            new_data[z+4]=str(n_discevent)+' "'+"Out_Files/"+actual_DLC_tmp+".out"+'"\n'
    with open(mlife_filename, 'w') as file:
        file.writelines( new_data )
        file.close
    print "MLife files were prepared."


    # set filename of Mlife_Roses-input-file
    mlife_roses_filename="Input_MLife_Roses.mlif"
        
    print "Manipulate Mlife_Roses input file..."
    # manipulate Mlife_Roses input-file by considering all copied Dlc_11 files        
    with open(mlife_roses_filename, 'r') as file:
            data = file.readlines()
    
    # initialize new file data
    new_data=[""]*(80+file_number_11)
    new_data[0:76]=data[0:76]
    
    # constant entries
    new_data[75]=str(file_number_11)+'  1.1   1.3   1.5   1.7    (Weibull-Weighted Normal Operation: NumNormFiles, PSF1, PSF2, PSF3, PSF4)\n'
    new_data[77+file_number_11]="0  1.1   1.3   1.5   1.7    (Weibull-Weighted Idling: NumIdleFiles, PSF1, PSF2, PSF3, PSF4)\n"
    new_data[78+file_number_11]="0  1.2   1.3   1.4   1.6    (Discrete Events: NumDiscFiles, PSF1, PSF2, PSF3, PSF4)\n"
    new_data[79+file_number_11]="==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.\n"
    
    z=76
    n_discevent=0
    for i in range(int(file_number)):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        if actual_DLC == '11':
            z=z+1
            new_data[z]='"'+"Out_Files/"+actual_DLC_tmp+".mod"+'"\n'
    with open(mlife_roses_filename, 'w') as file:
        file.writelines( new_data )
        file.close
    print "MLife_Roses files were prepared."
    
    
    # prepare Mextremes files
    print "Prepare files for Mextreme analysis."
    # set filename of Mextremes-input-file
    mextr_filename="Input_MExtremes.mext"
        
    print "Manipulate Mextremes input file..."
    # manipulate Mextremes input-file by considering all copied Dlc_13 files        
    with open(mextr_filename, 'r') as file:
            data = file.readlines()
    
    # initialize new file data
    Extreme_DLCs = 14
    n_MExtremes_files = file_number_11 + file_number_13 + file_number_14 + file_number_15 + file_number_21 + file_number_23 + file_number_32 + file_number_33 + file_number_42 + file_number_51 + file_number_61 + file_number_62 + file_number_63 + file_number_71
    new_data=[""]*(175 + n_MExtremes_files + Extreme_DLCs)
    new_data[0:174]=data[0:174]
    
    # constant entries
    new_data[175 + n_MExtremes_files + Extreme_DLCs -1]="==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.\n"
    
    z=174
    DLC_11 = 'False'
    DLC_13 = 'False'
    DLC_14 = 'False'
    DLC_15 = 'False'
    DLC_21 = 'False'
    DLC_23 = 'False'
    DLC_32 = 'False'
    DLC_33 = 'False'
    DLC_42 = 'False'
    DLC_51 = 'False'
    DLC_61 = 'False'
    DLC_62 = 'False'
    DLC_63 = 'False'
    DLC_71 = 'False'
    new_data[z-1]=str(Extreme_DLCs)+'                 NumDLCs           The number of Design Load Cases\n'
    for i in range(file_number):
        actual_DLC_tmp = file_list_raw[i]
        actual_DLC = actual_DLC_tmp[4:6]
        if actual_DLC == '11' or actual_DLC == '13' or actual_DLC == '14' or actual_DLC == '15' or actual_DLC == '21' or actual_DLC == '23' or actual_DLC == '32' or actual_DLC == '33' or actual_DLC == '42' or actual_DLC == '51' or actual_DLC == '61' or actual_DLC == '62' or actual_DLC == '63' or actual_DLC == '71':
            if actual_DLC == '11' and DLC_11 == 'False':
                new_data[z]=str(file_number_11)+'  '+str(1.25 * 1.2)+'   1.35   1.5   1.7    "DLC 1.1"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_11 = 'True'            
            if actual_DLC == '13' and DLC_13 == 'False':
                z = z + 1
                new_data[z]=str(file_number_13)+'  1.35   1.3   1.5   1.7    "DLC 1.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_13 = 'True'
            if actual_DLC == '14' and DLC_14 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_14)+'  1.35   1.3   1.5   1.7    "DLC 1.4"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_14 = 'True'
            if actual_DLC == '15' and DLC_15 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_15)+'  1.35   1.3   1.5   1.7    "DLC 1.5"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_15 = 'True'
            if actual_DLC == '21' and DLC_21 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_21)+'  1.35   1.3   1.5   1.7    "DLC 2.1"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_21 = 'True'
            if actual_DLC == '23' and DLC_23 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_23)+'  1.10   1.3   1.5   1.7    "DLC 2.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_23 = 'True'
            if actual_DLC == '32' and DLC_32 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_32)+'  1.35   1.3   1.5   1.7    "DLC 3.2"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_32 = 'True'
            if actual_DLC == '33' and DLC_33 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_33)+'  1.35   1.3   1.5   1.7    "DLC 3.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_33 = 'True'
            if actual_DLC == '42' and DLC_42 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_42)+'  1.35   1.3   1.5   1.7    "DLC 3.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_42 = 'True'
            if actual_DLC == '51' and DLC_51 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_51)+'  1.35   1.3   1.5   1.7    "DLC 3.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_51 = 'True'
            if actual_DLC == '61' and DLC_61 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_61)+'  1.35   1.3   1.5   1.7    "DLC 6.1"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_61 = 'True'
            if actual_DLC == '62' and DLC_62 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_62)+'  1.10   1.3   1.5   1.7    "DLC 6.2"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_62 = 'True'
            if actual_DLC == '63' and DLC_63 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_63)+'  1.35   1.3   1.5   1.7    "DLC 6.3"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_63 = 'True'
            if actual_DLC == '71' and DLC_71 == 'False':
                z = z + 1                
                new_data[z]=str(file_number_71)+'  1.10   1.3   1.5   1.7    "DLC 7.1"         (Normal Operation: NumDLCfiles(1), PSF1, PSF2, PSF3, PSF4)\n'
                DLC_71 = 'True'
            z = z + 1
            new_data[z]='"'+"Out_Files/"+actual_DLC_tmp+".out"+'"\n'
    with open(mextr_filename, 'w') as file:
        file.writelines( new_data )
        file.close
    print "Mextreme files were prepared."
    print "finished!"
    
# move output files to Out_Files folder
if Move_Out=='True':
    print "Move "+str(file_number)+ " files to Out_Files directory..."
    for i in range(len(file_list_raw)):
        if os.path.exists(root_FAST+file_list_raw[i]+".out"):
            shutil.move(root_FAST+file_list_raw[i]+".out", root_FAST+"/Out_Files//"+file_list_raw[i]+".out")
            print file_list_raw[i]+".out was moved...\n"
    print "All .out files were moved to /Out_Files directory."

# generate TurbSim wind files
if Calculate_bts == 'True' or Calculate_all_windfiles == 'True':
    print "Calculate all "+str(file_number_TurbSim)+ " TurbSim wind files..."
    p = Popen([os.path.join(root_TurbSim, filename_batch_all_TS + filename_batch_ending)], cwd=root_TurbSim)
    print "All "+str(file_number_TurbSim)+ " TurbSim wind files were calculated."

# generate IEC wind files
z61=0
z62=0
z63=0
z71=0
if Calculate_IEC == 'True' or Calculate_all_windfiles == 'True':
    print "Calculate all "+str(file_number_IECwind)+ " IEC wind files..."
    if file_number_IECwind % cores > 0:
        IEC_calculations = file_number_IECwind/cores + 1
    if file_number_IECwind % cores == 0:
        IEC_calculations = file_number_IECwind/cores
    
    z_IEC2 = 0
    for i in range(int(IEC_calculations)):
        for j in range(cores):
            print data_IEC[0,z_IEC2]
            print z_IEC2
            if data_IEC[0,z_IEC2] is not None and z_IEC2 < file_number_IECwind:
                with open(root_IECwind + "\Core" + str(j) + "\IEC.ipt", 'w') as file:
                    file.writelines( data_IEC[0:,z_IEC2] )
                    file.close
            z_IEC2 = z_IEC2 + 1
        p = Popen([os.path.join(root_IECwind, "_Run_all_IECwind.bat")], cwd = root_IECwind)
        time.sleep(sleep_time)
        z_IEC2 = z_IEC2 - cores
        for j in range(cores):
            if data_IEC[0,z_IEC2] is not None and z_IEC2 < file_number_IECwind:
                print data_IEC[0,z_IEC2]
                IEC_files = os.listdir(root_IECwind + "\Core" + str(j))
                for k in range(len(IEC_files)):
                    IEC_file_tmp = IEC_files[k]
                    if IEC_file_tmp[len(IEC_file_tmp)-4:len(IEC_file_tmp)] == ".WND":
                        while not os.path.isfile(root_IECwind + "\Core" + str(j) + "\\" + IEC_files[k]):
                           time.sleep(sleep_time) 
                        os.rename(root_IECwind + "\Core" + str(j) + "\\" + IEC_files[k], root_IECwind + "\Core" + str(j) + "\\" + name_IEC[z_IEC2] + '.WND')
                        shutil.move(root_IECwind + "\Core" + str(j) + "\\" + name_IEC[z_IEC2] + '.WND', root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND')
                        name_IEC_tmp = name_IEC[z_IEC2]
                        print name_IEC_tmp
                        
                        # use this in the case of deterministic wind for idling DLCs
                        # manipulate DLC 61 and 63 wind files afterwards
                        if name_IEC_tmp[4:6] == '61' or name_IEC_tmp[4:6] == '62' or name_IEC_tmp[4:6] == '63' or name_IEC_tmp[4:6] == '71':
                            with open(root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND', 'r') as file:
                                data_manipulate = file.readlines()
                            data_man_temp=["\n" for x in range(2)]
                            if name_IEC_tmp[4:6] == '61':
                                direc = directions_61[z61]
                                z61 += 1
                            if name_IEC_tmp[4:6] == '62':
                                direc = directions_62[z62]
                                z62 += 1
                            if name_IEC_tmp[4:6] == '63':
                                direc = directions_63[z63]
                                z63 += 1
                            if name_IEC_tmp[4:6] == '71':
                                direc = directions_71[z71]
                                z71 += 1
                            for m in range(2):
                                data_man_temp[m]=str(data_manipulate[10+m])
                                string=data_man_temp[m]
                                string=string[0:23] + str(direc) + '.000\t' + string[30:]
                                if m == 1:
                                    string = '  '+str(TMax_idling)+'\t' + string[10:]
                                data_man_temp[m]=string
                                data_manipulate[10+m]=data_man_temp[m]
                            with open(root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND', 'w') as file:
                                file.writelines( data_manipulate )
                                file.close
                        
                        
                        # manipulate power law exponent in each IEC file
                        if name_IEC_tmp[4:6] == '14' or name_IEC_tmp[4:6] == '15' or name_IEC_tmp[4:6] == '23' or name_IEC_tmp[4:6] == '32' or name_IEC_tmp[4:6] == '33' or name_IEC_tmp[4:6] == '42':
                            with open(root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND', 'r') as file:
                                    data_manipulate = file.readlines()
                                    file.close()
                            if name_IEC_tmp[4:6] == '14' or name_IEC_tmp[4:6] == '15':
                                data_man_temp=["\n" for x in range(len(data_manipulate)-15)]
                                for m in range(len(data_man_temp)):
                                    data_man_temp[m]=str(data_manipulate[15+m])
                                    string=data_man_temp[m]
                                    string=string.replace("0.200",str(pwr_law_exp))
                                    string = string[0:20] + str(round(float(string[20:30]) + float(name_IEC_tmp[7:10]),3)) + '\t' + string[30:]
                                    data_man_temp[m]=string
                                    data_manipulate[15+m]=data_man_temp[m]
                            if name_IEC_tmp[4:6] == '23' or name_IEC_tmp[4:6] == '32' or name_IEC_tmp[4:6] == '33' or name_IEC_tmp[4:6] == '42':
                                data_man_temp=["\n" for x in range(len(data_manipulate)-16)]
                                for m in range(len(data_man_temp)):
                                    data_man_temp[m]=str(data_manipulate[16+m])
                                    string=data_man_temp[m]
                                    string=string.replace("0.200",str(pwr_law_exp))
                                    data_man_temp[m]=string
                                    data_manipulate[16+m]=data_man_temp[m]
                            with open(root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND', 'w') as file:
                                file.writelines( data_manipulate )
                                file.close
                        
                        # copy IEC wind files to destination folder in the FAST directory
                        shutil.copyfile(root_IECwind + 'All_Windfiles\\' + name_IEC[z_IEC2] + '.WND', root_FAST+'InflowWind\Wind\\' + name_IEC[z_IEC2] + '.WND')
            z_IEC2 = z_IEC2 + 1
    print "All "+str(file_number_IECwind)+ " IEC wind files were calculated."

# Get uncalculated .fst files and try to calculate them again
z_exist = 0
remaining_files = []
if Calc_Rest_FAST == 'True':
    for i in range(len(file_list_raw)):
        if not os.path.exists(root_FAST+"/Out_Files//"+file_list_raw[i]+".out"):
            z_exist = z_exist + 1
            remaining_files.append(file_list_raw[i])
            print "File ",file_list_raw[i]+".out", " does not exist."
    
    print "Try to calculate these ",str(z_exist)," files again? [y/n]"
    re_calc_answ = raw_input()
    if re_calc_answ == 'y':
        file_number_rest = len(remaining_files)
        # calculate distribution of jobs
        job_distr=int(np.ceil(file_number_rest/cores))
        job_distr_rest = job_distr - file_number_rest/cores
        job_distr_overflow = (file_number_rest/cores)%1
        FAST_overflow_depot = 0
        
        # generate initialized batch files for distributed jobs
        commands_each_batch_FAST=["\n" for x in range(job_distr+2)]
        z1=0
        for i in range(cores):
            # FAST
            filename_each_batch_FAST=root_FAST+"_Run_FAST_"+str(i+1)+filename_batch_ending
            commands_each_batch_FAST=["\n" for x in range(job_distr+2)]
            if FAST_overflow_depot <= 0:
                job_distr_current = job_distr
                FAST_overflow_depot = FAST_overflow_depot + job_distr_rest
            else:
                job_distr_current = job_distr - 1
                FAST_overflow_depot = FAST_overflow_depot - job_distr_overflow
            for j in range(job_distr_current):
                if z1<file_number_rest:
                    string="FAST_x64 "+str(remaining_files[z1])+filename_FAST_end+"\n"
                    commands_each_batch_FAST[j]=string
                if z1>=file_number_rest:
                    string="\n"
                    commands_each_batch_FAST[j]=string
                z1=z1+1
            with open(filename_each_batch_FAST, 'w') as file:
                    file.writelines( commands_each_batch_FAST )
                    file.close
        
        print "New batch files for FAST execution created."

# Get uncalculated TurbSim .inp files and try to calculate them again
z_exist = 0
remaining_files = []
if Calc_Rest_TurbSim == 'True':
    for i in range(len(file_list)):
        if (not os.path.exists(root_FAST+"/InflowWind/Wind//"+file_list_raw[i]+".bts")) and (file_list[i][-3:] == 'inp'):
            z_exist = z_exist + 1
            remaining_files.append(file_list_raw[i])
            print "File ",file_list_raw[i]+".bts", " does not exist."
    
    print "Try to calculate these ",str(z_exist)," files again? [y/n]"
    re_calc_answ = raw_input()
    if re_calc_answ == 'y':
        file_number_rest = len(remaining_files)
        # calculate distribution of jobs
        job_distr=int(np.ceil(file_number_rest/cores))
        job_distr_rest = job_distr - file_number_rest/cores
        job_distr_overflow = (file_number_rest/cores)%1
        TurbSim_overflow_depot = 0
        
        # generate initialized batch files for distributed jobs
        commands_each_batch_TurbSim=["\n" for x in range(job_distr+2)]
        z1=0
        for i in range(cores):
            # TurbSim
            filename_each_batch_TurbSim=root_TurbSim+"_Run_TurbSim_"+str(i+1)+filename_batch_ending
            commands_each_batch_TurbSim=["\n" for x in range(job_distr+2)]
            if TurbSim_overflow_depot <= 0:
                job_distr_current = job_distr
                TurbSim_overflow_depot = TurbSim_overflow_depot + job_distr_rest
            else:
                job_distr_current = job_distr - 1
                TurbSim_overflow_depot = TurbSim_overflow_depot - job_distr_overflow
            for j in range(job_distr_current):
                if z1<file_number_rest:
                    string="TurbSim_x64 "+str(remaining_files[z1])+Input_ending+"\n"
                    commands_each_batch_TurbSim[j]=string
                if z1>=file_number_rest:
                    string="\n"
                    commands_each_batch_TurbSim[j]=string
                z1=z1+1
            with open(filename_each_batch_TurbSim, 'w') as file:
                    file.writelines( commands_each_batch_TurbSim )
                    file.close
        
        print "New batch files for TurbSim execution created."
    
        
    