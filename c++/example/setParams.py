import json
import numpy as np
import pickle
import sys
''' This file sets the colossus parameters and writes the param files for pr-calc. It requires a global parameter ini file and is run as:

python setParams.py globalParams.ini '''


inifile=open(sys.argv[1], 'r')

# Reading the inifile
lines = inifile.readlines()

params =[]
names = []
# Pulling lines from the line and saving the names and parameter values
for line in lines:
    if '#' in line:
        line=''
    else:
        if '=' in line:
            line=line.strip('\n')
            name = line.split('=')[0]
            param=line.split('=')[1]
            params.append(param)
            names.append(name)


print('The names of all parameters in the global file are:')
print(names, params)


#Set cosmology
cosmoparams = {names[i]:params[i] for i in range(7)}
with open('parameterfiles/param.col', 'wb') as handle:
    pickle.dump(cosmoparams, handle)
    
# Set ics
ics = open('parameterfiles/tempParam.ics', 'r')
icslines = ics.readlines()
ics_out = open('parameterfiles/param.ics', 'w')
paramnums_ics=[2]

# Write out the parameter file for ics
for line in icslines:
    if '#' in line:
        ics_out.write(line)
    else:
        splitline=line.split()
        if(len(splitline)==2):
            if(splitline[0]=='h'):
                line=' '.join([splitline[0], '  ', params[2]])
            if(splitline[0]=='BoxSize'):
                line=' '.join([splitline[0], '  ', params[8]])
            if(splitline[0]=='N'):
                line=' '.join([splitline[0], '  ', params[9]])
        else:
            line=line   
        ics_out.write(line)
        ics_out.write('\n')
        
ics_out.close()

# Set ics
d2z = open('parameterfiles/tempParam.d2z', 'r')
d2zlines = d2z.readlines()
d2z_out = open('parameterfiles/param.d2z', 'w')

# Write out the parameter file for d2z
for line in d2zlines:
    if '#' in line:
        d2z_out.write(line)
    else:
        splitline=line.split()
        if(len(splitline)==2):
            if(splitline[0]=='Omegab'):
                line=' '.join([splitline[0], '  ', params[4]])
            if(splitline[0]=='Omegal'):
                line=' '.join([splitline[0], '  ', str(1-float(params[3]))])
            if(splitline[0]=='Omegam'):
                line=' '.join([splitline[0], '  ', params[3]])
            if(splitline[0]=='h'):
                line=' '.join([splitline[0], '  ', params[2]])
            if(splitline[0]=='sigma8'):
                line=' '.join([splitline[0], '  ', params[5]])
            if(splitline[0]=='ns'):
                line=' '.join([splitline[0], '  ', params[6]])
            if(splitline[0]=='w'):
                line=' '.join([splitline[0], '  ', params[7]])
            if(splitline[0]=='BoxSize'):
                line=' '.join([splitline[0], '  ', params[8]])
            if(splitline[0]=='N'):
                line=' '.join([splitline[0], '  ', params[9]])
            if(splitline[0]=='zInit'):
                line=' '.join([splitline[0], '  ', params[10]])

            if(splitline[0]=='Rmin'):
                line=' '.join([splitline[0], '  ', params[19]])
            if(splitline[0]=='Nscales'):
                line=' '.join([splitline[0], '  ', params[20]])

        else:
            line=line   
        d2z_out.write(line)
        d2z_out.write('\n')
        
d2z_out.close()

# Set ics
fsm = open('parameterfiles/tempParam.fsm', 'r')
fsmlines = fsm.readlines()
fsm_out = open('parameterfiles/param.fsm', 'w')

# Write out the parameter file for d2z
for line in fsmlines:
    if '#' in line:
        fsm_out.write(line)
    else:
        splitline=line.split()
        if(len(splitline)==2):
            if(splitline[0]=='Omegab'):
                line=' '.join([splitline[0], '  ', params[4]])
            if(splitline[0]=='Omegal'):
                line=' '.join([splitline[0], '  ', str(1-float(params[3]))])
            if(splitline[0]=='Omegam'):
                line=' '.join([splitline[0], '  ', params[3]])
            if(splitline[0]=='h'):
                line=' '.join([splitline[0], '  ', params[2]])
            if(splitline[0]=='Sigma8'):
                line=' '.join([splitline[0], '  ', params[28]])
            if(splitline[0]=='ns'):
                line=' '.join([splitline[0], '  ', params[6]])
            if(splitline[0]=='w'):
                line=' '.join([splitline[0], '  ', params[7]])
            if(splitline[0]=='BoxSize'):
                line=' '.join([splitline[0], '  ', params[8]])
            if(splitline[0]=='NPixels'):
                line=' '.join([splitline[0], '  ', params[17]])
            if(splitline[0]=='zInit'):
                line=' '.join([splitline[0], '  ', params[10]])
            if(splitline[0]=='N'):
                line=' '.join([splitline[0], '  ', params[9]])
            if(splitline[0]=='NRedshifts'):
                line=' '.join([splitline[0], '  ', params[12]])
            if(splitline[0]=='fov'):
                line=' '.join([splitline[0], '  ', params[14]])
            if(splitline[0]=='theta'):
                line=' '.join([splitline[0], '  ', params[15]])
            if(splitline[0]=='phi'):
                line=' '.join([splitline[0], '  ', params[16]])
            if(splitline[0]=='InitialRedshift'):
                line=' '.join([splitline[0], '  ', params[29]])
            if(splitline[0]=='FinalRedshift'):
                line=' '.join([splitline[0], '  ', params[30]])

        else:
            line=line
        fsm_out.write(line)
        fsm_out.write('\n')

fsm_out.close()


# Set ics
lmb = open('parameterfiles/tempParam.lmb', 'r')
lmblines = lmb.readlines()
lmb_out = open('parameterfiles/param.lmb', 'w')

# Write out the parameter file for d2z
for line in lmblines:
    if '#' in line:
        lmb_out.write(line)
    else:
        splitline=line.split()
        if(len(splitline)==2):
            if(splitline[0]=='kmin'):
                line=' '.join([splitline[0], '  ', params[21]])
            if(splitline[0]=='kmax'):
                line=' '.join([splitline[0], '  ', params[22]])
            if(splitline[0]=='nk'):
                line=' '.join([splitline[0], '  ', params[23]])
            if(splitline[0]=='zmin'):
                line=' '.join([splitline[0], '  ', params[24]])
            if(splitline[0]=='zmax'):
                line=' '.join([splitline[0], '  ', params[25]])
            if(splitline[0]=='nz'):
                line=' '.join([splitline[0], '  ', params[26]])
            if(splitline[0]=='BoxSize'):
                line=' '.join([splitline[0], '  ', params[8]])
            if(splitline[0]=='N'):
                line=' '.join([splitline[0], '  ', params[9]])
            if(splitline[0]=='zInit'):
                line=' '.join([splitline[0], '  ', params[10]])
            if(splitline[0]=='nl'):
                line=' '.join([splitline[0], '  ', params[27]])

        else:
            line=line
        lmb_out.write(line)
        lmb_out.write('\n')

lmb_out.close()
