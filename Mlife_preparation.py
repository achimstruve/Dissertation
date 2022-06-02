import numpy as np
import shutil

# set root path
root="C:\Users\onrok\Desktop\Achim\Promotion\Calculation\Tower_Studies\Fatigue_Calculation\\"

# set filename of Mlife-input-file
mlife_filename="Dlc12.mlif"


# define all files to be considered by Mlife #
# DLC 1.1
directories=np.array(['S1','S2','S3','S4','S5','S6'])
wind_speeds=np.array(['03','05','07','09','11','13','15','17','19','21','23','25'])

print "copy "+str(len(directories)*len(wind_speeds))+ " files to Outs directory..."
for i in directories:
    for j in wind_speeds:
        shutil.copyfile(root+i+"\\"+i+"_"+j+".out", root+"Outs\\"+i+"_"+j+".out")
        print i+"_"+j+" was copied...\n"

print "Manipulate Mlife-file..."
# manipulate Mlife input-file by considering all copied files        
with open(mlife_filename, 'r') as file:
        data = file.readlines()

# initialize new file data
new_data=[""]*(134+len(directories)*len(wind_speeds))
new_data[0:130]=data[0:130]

# constant entries
new_data[131+len(directories)*len(wind_speeds)]="0  1.1   1.3   1.5   1.7    (Weibull-Weighted Idling: NumIdleFiles, PSF1, PSF2, PSF3, PSF4)\n"
new_data[132+len(directories)*len(wind_speeds)]="0  1.2   1.3   1.4   1.6    (Discrete Events: NumDiscFiles, PSF1, PSF2, PSF3, PSF4)\n"
new_data[133+len(directories)*len(wind_speeds)]="==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.\n"

z=130
for i in directories:
    for j in wind_speeds:
        z=z+1
        new_data[z]='"'+"Outs/"+i+"_"+j+".out"+'"\n'
        with open(mlife_filename, 'w') as file:
            file.writelines( new_data )
            file.close
print "finished!"