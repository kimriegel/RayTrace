# Receiver Point 
# Python 3.7.0 

import numpy as np

class Receiver:
    """
    Receiver class with attributes of position (x,y,z)
    Functionality for pressure not yet added.
    """
    planenum=1
    planename1='Single Point'
    arraysize=0     #Number of receivers. Supposed to be 5 for this test    # I use it differently now, and only in one function
    sizex=2
    sizey=2
    sizez=1
    initial_frequency = None    # Gives this value to all receivers

    rList = [] #See append_list

    # I personally prefer writing Receiver.Array to Receiver.receiverarray
    #Array = np.array([None])    
    # This is now completely overlapped by rlist

    def __init__(self,position):
        """
        Create and defines position of receiver
        
        Works automatically when class is called
        """
        self.position=np.array(position)
        self.recNumber = Receiver.arraysize         #planned for debugging but we don't seem to use it
        Receiver.rList.append(self) #See append_list
            # Initial values 
        self.pressure = 0
        self.magnitude = 0
        self.direction = 0

        Receiver.arraysize += 1

    def on_Hit(self,amplitude,phase):
        """ 
        My version of old receiver hit function. 
        Modifies direction and magnitude of rays with respect to each receiver
        """

        temp1 = abs(self.magnitude) * np.exp(XJ*self.direction)
        temp2 = abs(amplitude[:])   * np.exp(XJ*phase[:])
        temp3 = temp1 + temp2 

        self.magnitude =  abs(temp3)                                 
        self.direction =  np.arctan2(np.imag(temp3) , np.real(temp3))
        # See bug log 3/13 for what happened with positions checks

    def SphereCheck(self,Sr2,F,veci):
        """
        This function performs a check whether a ray hits a sphere.  If
        it does hit then the function returns the distance to the sphere
        """
        # Sc is receiver position
        # sr2 is radius2
        HUGE=1000000.0

        Sc = self.position
        OC = Sc - veci 
        L2OC = np.dot(OC,OC)    #Equivalent?
        tca = np.dot(OC,F)
        t2hc = Sr2 - L2OC + (tca**2)
        if L2OC == Sr2:    
            dx = HUGE
        elif tca < 0.0:
            dx = HUGE
        elif t2hc < 0.0:
            dx = HUGE
        else:
            dx = tca - (t2hc**(1/2))
        #Hey future me, remember to delete one of these later

        self.dx = dx
        self.dxreceiver = 0
        return  dx 

    def timeReconstruct(self,sizefft):
        """
        This Function computes the timesignal from a given fft.  It writes the
        time signal to an array
        Timearray is now defined in here. Do not call it.
        """
        # Create the complex array to feed into the inverse fft function
        if self.magnitude[0] == 0.0:            # If first magnitude is zero then all timesignal is zero
            self.signal = np.zeros(sizefft) 

        else:                                   # If not then calculate the timesignal 
            tempfft = abs(self.magnitude[:]) * np.exp(XJ*self.direction)
            tempfft = np.append(0,tempfft)
            self.signal=np.fft.irfft(tempfft,int(sizefft)) * sizefft

    @classmethod
    def timeHeader(cls,f,time,w):
        """
        This function prints the header between each time set
        Time is the real time that the event happens
        omega (w) is the signal in that receiver at the specified time
        """
        #This function prints the header between each time set
        f.write('ZONE T=" %s "\n' %(planename, ) )  #this worked in the command line
        f.write('STRANDID=1, SOLUTIONTIME= %d \n'%(time,) )                                     #timearray tiempo/PF.Fs 
        f.write('I= %d\t J= %d\t K=%d\t ZONETYPE=Ordered\n' %(sizex, sizey, sizez))
        f.write('DATAPACKING=POINT\n')
        f.write('DT=(SINGLE SINGLE SINGLE SINGLE )\n')

        for R in cls.rList:
          print('\t%f\t%f\t%f\t%f' %(R.position[0],R.position[1],R.position[2],R.signal[w]),file=f)    #time signal

        return 

    @classmethod
    def initialize(cls,ipfile):
        """
        Reads in receiver points from the txt file and translates them to be used in our receiver method
        Receivers are automatically initialized from given inputfile
        """
        with open(ipfile) as vertex:        #Read in from file
            rho = np.genfromtxt(vertex)
        cls.Array = rho
        pointNo = rho.shape[0]   #Create an itterater, Number of Receiver Points
        for r in range(pointNo): #Translate to usable format
            Receiver(rho[r,:])             # Come back later to change the needed inputs for method
        print('initialized receivers')


if __name__ == "__main__" :
    #Only exists for checking bugs with specific parameters now
    import RayTrace
    #Receiver.initialize("PointReceivers.txt")
    #print(Receiver.rList)
    #print(Receiver.arraysize)
    pass

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
    XJ=complex(0,1)

