import numpy as np
import matplotlib.pyplot as plt

from Parameterfile import Fs


with open("OutputTest.txt") as ipFile:
    pressure = np.loadtxt(ipFile)       # y axis

K = len(pressure)
time = np.arange(K) /Fs                 # x axis

plt.plot(time,pressure,'r--')
    # Labeling axes
plt.xlabel('Time')
plt.ylabel('Pressure')
plt.title('Pressure vs Time of Receiver 2')
    # Setting boundaries
#plt.axis([0,(6 / (1000)**3),0,30])
plt.savefig('TestGraph.png')
plt.show()

print('Done with this file')