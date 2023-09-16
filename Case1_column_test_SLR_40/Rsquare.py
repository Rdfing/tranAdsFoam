# import numpy library as np
import numpy as np

# some parameter
C_initial = 0.0005
U_super = 0.000666666667 # m/s
bed_length = 4.33 /100 # m
unit_time = bed_length/U_super # time for single bed length/volume

# experiment data file
filename="SLR_40.csv"
# use skiprows to skip rows
expData = np.loadtxt(filename, delimiter=",", skiprows=1)
expData[:,0] = expData[:,0]*unit_time # conver to the write time unit

# CFD data file
filename="postProcessing/surfaceMonitor_outlet_T/0/surfaceFieldValue.dat"
# use skiprows to skip rows
CFDData = np.loadtxt(filename, skiprows=6)
CFDData[:,1] = CFDData[:,1]/C_initial # normalized the concentration

# interpolate the data to the uniform grid
xref = np.linspace(min(expData[:,0]),max(expData[:,0]),100)
expData_yref = np.interp(xref, expData[:,0], expData[:,1])
CFDData_yref = np.interp(xref, CFDData[:,0], CFDData[:,1])

# compute the Error and Rsqure
Error = np.nanmean(np.power((expData_yref-CFDData_yref),2))
R_square = 1-np.sum(np.power((expData_yref-CFDData_yref),2))/np.sum(np.power((expData_yref-np.nanmean(expData_yref)),2))

# print Error
# print R_square

#
file1 = open('Rsquare.dat', 'w') 
file1.writelines("Error R_square\n") 
file1.write('%0.5f\n' % Error)    
file1.write('%0.5f\n' % R_square)
file1.close() 


import matplotlib.pyplot as plt
plt.plot(xref, expData_yref, 'o',label='exp')
plt.plot(CFDData[:,0], CFDData[:,1], linestyle='-', marker='^',label='CFD (Error=%0.3f, Rsqure=%0.3f)' % (Error,R_square))
# plt.plot(xref, CFDData_yref, linestyle='-', marker='^',label='CFD (Error=%0.3f, Rsqure=%0.3f)' % (Error,R_square))
# plt.show()
plt.legend(loc='upper left', frameon=False)
plt.xlabel('Time (s)')
plt.ylabel('C/C0')
plt.savefig('comparsion.png')

