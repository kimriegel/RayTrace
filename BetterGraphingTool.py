#Cleaned Up Graphing tool using functions instead of endless loops and copy paste

import numpy as np
import matplotlib.pyplot as plt

c = 'OrigFor_radiosity_50_006.dat' #Putting Fortran file into a variable for function

a = 'ISO_Comp_Build_PD0.50.06.txt'  #putting Python file into a variable for function

#Python Reader Function 
def PythonReader(a):

    time_arrayP = np.array([])
    pressure1_arrayP = np.array([])
    pressure2_arrayP = np.array([])
    pressure3_arrayP = np.array([])
    pressure4_arrayP = np.array([])
    pressure5_arrayP = np.array([])

    with open(a, 'r', encoding='utf8') as f:
        for line in f:
            if 'ZONE T=" Single Point "' in line:
                time = f.readline()
                time = time.split(' ')
                time = float(time[2])
                time_arrayP = np.append(time_arrayP, time)
    f.close()

    with open(a, 'r', encoding='utf8') as f:
        for line in f:
            if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
                # Receiver1
                pressure = f.readline()
                pressure1 = pressure.split('\t')
                pressure1 = float(pressure1[4])
                pressure1_arrayP = np.append(pressure1_arrayP, pressure1)
                # Receiver2
                pressure = f.readline()
                pressure2 = pressure.split('\t')
                pressure2 = float(pressure2[4])
                pressure2_arrayP = np.append(pressure2_arrayP, pressure2)
                # Receiver3
                pressure = f.readline()
                pressure3 = pressure.split('\t')
                pressure3 = float(pressure3[4])
                pressure3_arrayP = np.append(pressure3_arrayP, pressure3)
                # Receiver4
                pressure = f.readline()
                pressure4 = pressure.split('\t')
                pressure4 = float(pressure4[4])
                pressure4_arrayP = np.append(pressure4_arrayP, pressure4)
                # Receiver5
                pressure = f.readline()
                pressure5 = pressure.split('\t')
                pressure5 = float(pressure5[4])
                pressure5_arrayP = np.append(pressure5_arrayP, pressure5)
    f.close()
    return time_arrayP, pressure1_arrayP, pressure2_arrayP, pressure3_arrayP, pressure4_arrayP, pressure5_arrayP

#Fortran Reader Function

def FortranReader(c):
    # Fortran Reader Code

    # Creates the neccesary blank arrays for the program.

    time_arrayF = np.array([])
    pressure1_arrayF = np.array([])
    pressure2_arrayF = np.array([])
    pressure3_arrayF = np.array([])
    pressure4_arrayF = np.array([])
    pressure5_arrayF = np.array([])

    # Mirrors the code that appends values to the time array, but is tailored for the different spacing and readout of the
    # Fortran text file as opposed to the Python text file.

    with open(c, 'r', encoding='utf8') as f:
        for line in f:
            if 'ZONE T="Single Point        "' in line:
                time = f.readline()
                time = time.split('  ')
                time = float(time[1])
                time_arrayF = np.append(time_arrayF, time)
    f.close()

    # Same thing as above, but for the pressure receivers.

    with open(c, 'r', encoding='utf8') as f:
        for line in f:
            if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
                # Receiver 1
                pressure = f.readline()
                pressure1 = pressure.split('     ')
                pressure1 = float(pressure1[3])
                pressure1_arrayF = np.append(pressure1_arrayF, pressure1)
                # Receiver2
                pressure = f.readline()
                pressure2 = pressure.split('     ')
                pressure2 = float(pressure2[3])
                pressure2_arrayF = np.append(pressure2_arrayF, pressure2)
                # Receiver3
                pressure = f.readline()
                pressure3 = pressure.split('     ')
                pressure3 = float(pressure3[3])
                pressure3_arrayF = np.append(pressure3_arrayF, pressure3)
                # Receiver4
                pressure = f.readline()
                pressure4 = pressure.split('     ')
                pressure4 = float(pressure4[3])
                pressure4_arrayF = np.append(pressure4_arrayF, pressure4)
                # Receiver5
                pressure = f.readline()
                pressure5 = pressure.split('     ')
                pressure5 = float(pressure5[3])
                pressure5_arrayF = np.append(pressure5_arrayF, pressure5)
    f.close()
    return time_arrayF, pressure1_arrayF, pressure2_arrayF, pressure3_arrayF, pressure4_arrayF, pressure5_arrayF

#FortranRead

############### Fortran Plotter ####################
timef, fr1, fr2, fr3, fr4, fr5 = FortranReader(c)

ReceiverListF = fr1, fr2, fr3, fr4, fr5 

i = 1 

for x in ReceiverListF:
    #Plotting Receivers Fortran
    plt.figure(layout='constrained')
    plt.plot(timef, x, color='blue', alpha=0.5)
    plt.title('Fortran Pressure vs Time Receiver ' + str(i))
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.xlim(0.0, 0.8)
    plt.grid()
    plt.show()
    i = i + 1 

############ End Fortran Plotter ###############

############ Python Plotter ################

timep, pr1, pr2, pr3, pr4, pr5 = PythonReader(a)

ReceiverListP = pr1, pr2, pr3, pr4, pr5

i = 1

for x in ReceiverListP:
    #Plotting Receivers Python
    plt.figure(layout='constrained')
    plt.plot(timep, x, color='red', alpha=0.5)
    plt.title('Python Pressure vs Time Receiver ' + str(i))
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.xlim(0.0, 0.8)
    plt.grid()
    plt.show()
    i = i + 1 

############### End Python Plotter ##############


################ Fortran vs Python Plotter ################

i = 1

for x, z in zip(ReceiverListP, ReceiverListF):
    plt.figure(layout = 'constrained')
    plt.plot(timep, x, color = 'red', alpha = 0.5)
    plt.plot(timef, z, color = 'blue', alpha = 0.5)
    plt.title('Python vs Fortran, Pressure vs Time Receiver ' + str(i))
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.xlim(0.0, 0.8)
    plt.grid()
    plt.show()
    i = i + 1 

############################################################



