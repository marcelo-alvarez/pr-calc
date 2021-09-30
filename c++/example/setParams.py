import json
import numpy as np
import pickle
import sys


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


print(names)


#Set cosmology
cosmoparams = {names[i]:params[i] for i in range(6)}
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
        for paramnum in paramnums_ics:
            line.replace('XX', params[paramnum])
            print(line)
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
        print(line)
        d2z_out.write(line)
        d2z_out.write('\n')
        
d2z_out.close()

    
