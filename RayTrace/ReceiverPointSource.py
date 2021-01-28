# Receiver Point 
# Python 3.7.0 

import numpy as np


class Receiver:

    """
    Receiver class with attributes of position (x,y,z)
    Functionality for pressure not yet added.
    """
    planenum = 1
    planename1 = 'Single Point'
    arraysize = 1     # Number of receivers. Supposed to be 5 for this test
    # I use it only for printing receiver label on graphs
    sizex = 2
    sizey = 2
    sizez = 1
    initial_frequency = None    # Gives this value to all receivers
    rList = []  # See append_list

    # I personally prefer writing Receiver.Array to Receiver.receiverarray
    # Array = np.array([None])
    # This is now completely overlapped by rlist

    def __init__(self, position):
        """
        Create and defines position of receiver
        
        Works automatically when class is called
        """
        self.position = np.array(position)
        self.recNumber = int(Receiver.arraysize)         # planned for debugging but we don't seem to use it
        Receiver.rList.append(self)  # See append_list
        # Initial values
        self.pressure = 0
        self.magnitude = 0
        self.direction = 0

        Receiver.arraysize += 1

    def on_hit(self, amplitude, phase):
        """ 
        My version of old receiver hit function. 
        Modifies direction and magnitude of rays with respect to each receiver
        """
        xj = complex(0, 1)
        # print('initiating hit function')

        temp1 = abs(self.magnitude) * np.exp(xj*self.direction)
        temp2 = abs(amplitude[:]) * np.exp(xj*phase[:])
        # print(temp2.shape)
        # print(list(temp2[-20:]))
        # print(list(temp2[:]))
        # print(list(phase))
        temp3 = temp1 + temp2 

        self.magnitude = abs(temp3)
        self.direction = np.arctan2(np.imag(temp3), np.real(temp3))
        # See bug log 3/13 for what happened with positions checks

    def sphere_check(self, sr_2, f, veci):

        # This function performs a check whether a ray hits a sphere.  If
        # it does hit then the function returns the distance to the sphere

        # Sc is receiver position
        # sr2 is radius2
        huge = 1000000.0

        s_c = self.position
        oc = s_c - veci
        # L2OC = np.sum( (oc*oc),axis=1)
        l2_oc = np.dot(oc, oc)    # Equivalent?
        tca = np.dot(oc, f)
        t2hc = sr_2 - l2_oc + (tca**2)
        if l2_oc == sr_2:
            dx = huge
        elif tca < 0.0:
            dx = huge
        elif t2hc < 0.0:
            dx = huge
        else:
            dx = tca - (t2hc**(1/2))
        # Hey future me, remember to delete one of these later

        self.dx = dx
        return dx

    def time_reconstruct(self, sizefft):

        # This Function computes the timesignal from a given fft.  It writes the
        # time signal to an array
        # Timearray is now defined in here. Do not call it.
        xj = complex(0, 1)

        # Create the complex array to feed into the inverse fft function
        # Create complex array and compute inverse fft first attempt Python
        if self.magnitude[0] == 0.0:            # If first magnitude is zero then all timesignal is zero
            self.signal = np.zeros(sizefft) 

        else:                                   # If not then calculate the timesignal 
            tempfft = abs(self.magnitude[:]) * np.exp(xj*self.direction)
            tempfft = np.append(0, tempfft)
            # print(list(np.real(tempfft)))
            # Compute ifft using numpy
            self.signal = np.fft.irfft(tempfft, int(sizefft)) * sizefft
            # print('time signal @receiver 2: \n',list(self.signal[:100]))

    # def timeheader(cls,f,time,sizex,sizey,sizez,planename):
    @classmethod
    def time_header(cls, f, time, w):

        # This function prints the header between each time set
        # Time is the real time that the event happens
        # omega (w) is the signal in that receiver at the specified time

        # time = pass
        # This function prints the header between each time set
        f.write('ZONE T=" %s "\n' % (planename, ))  # this worked in the command line
        f.write('STRANDID=1, SOLUTIONTIME= %f \n' % time)                     # timearray tiempo/PF.Fs
        f.write('I= %d\t J= %d\t K=%d\t ZONETYPE=Ordered\n' % (sizex, sizey, sizez))
        f.write('DATAPACKING=POINT\n')
        f.write('DT=(SINGLE SINGLE SINGLE SINGLE )\n')

        for R in cls.rList:
            print('\t%f\t%f\t%f\t%f' % (R.position[0], R.position[1], R.position[2], R.signal[w]), file=f)
            # time signal
        return 

    @classmethod
    def initialize(cls, ipfile):

        # Reads in receiver points from the txt file and translates them to be used in our receiver method
        # Receivers are automatically initialized from given inputfile

        with open(ipfile) as vertex:        # Read in from file
            rho = np.genfromtxt(vertex)
        cls.Array = rho
        point_no = rho.shape[0]   # Create an itterater, Number of Receiver Points
        for r in range(point_no):  # Translate to usable format
            Receiver(rho[r, :])             # Come back later to change the needed inputs for method
        print('initialized receivers')


# Initializing receivers
def initialize_receivers():     

    # Create Individual receivers
    # For debugging, do not use

    r_1 = Receiver((93.4428213, 28.8397178, 0.151))
    r_2 = Receiver((64.5832, 28.5998, 0.151))
    r_3 = Receiver((64.5832, 28.5998, 7.9423))
    r_4 = Receiver((-2.40793, 31.5003401, 0.151))
    r_5 = Receiver((75.11005, 28.4945787, 0.151))

    return 


# For backwards compatibility
if __name__ != "__main__":      # only runs if opened outside this file
    planenum = 1
    planename1 = 'Single Point'
    planename = 'Single Point'
    arraysize1 = Receiver.arraysize
    sizex = 2
    sizey = 2
    sizez = 1
    sizex1 = sizex
    sizey1 = sizey
    sizez1 = sizez

    pass
