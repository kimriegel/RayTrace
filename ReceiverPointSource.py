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
        XJ = complex(0,1)
        #print('initiating hit function')

        temp1 = abs(self.magnitude) * np.exp(XJ*self.direction)
        temp2 = abs(amplitude[:])   * np.exp(XJ*phase[:])
        #print(temp2.shape)
        #print(list(temp2[-20:]))
        #print(list(temp2[:]))
        #print(list(phase))
        temp3 = temp1 + temp2 

        self.magnitude =  abs(temp3)                                 
        self.direction =  np.arctan2(np.imag(temp3) , np.real(temp3))
        # See bug log 3/13 for what happened with positions checks

    def SphereCheck(self,Sr2,F,veci):
        '''
        This function performs a check whether a ray hits a sphere.  If
        it does hit then the function returns the distance to the sphere
        '''
        # Sc is receiver position
        # sr2 is radius2
        HUGE=1000000.0

        Sc = self.position
        OC = Sc - veci 
        #L2OC = np.sum( (OC*OC),axis=1)    
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
        '''
        This Function computes the timesignal from a given fft.  It writes the
        time signal to an array
        Timearray is now defined in here. Do not call it.
        '''

        XJ=complex(0,1)
        # Create the complex array to feed into the inverse fft function
        # Create complex array and compute inverse fft first attempt Python
        if self.magnitude[0] == 0.0:            # If first magnitude is zero then all timesignal is zero
            self.signal = np.zeros(sizefft) 

        else:                                   # If not then calculate the timesignal 
            tempfft = abs(self.magnitude[:]) * np.exp(XJ*self.direction)
            tempfft = np.append(0,tempfft)
            #print(list(np.real(tempfft)))
            # Compute ifft using numpy
            self.signal=np.fft.irfft(tempfft,int(sizefft)) * sizefft
            #print('time signal @receiver 2: \n',list(self.signal[:100]))

    #def timeheader(cls,f,time,sizex,sizey,sizez,planename):
    @classmethod
    def timeHeader(cls,f,time,w):
        """
        This function prints the header between each time set
        Time is the real time that the event happens
        omega (w) is the signal in that receiver at the specified time

        """
        #time = pass
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

# Initializing receivers
def initialize_receivers():     
    """
    Create Individual receivers
    For debugging, do not use
    """
    R1 = Receiver((93.4428213,28.8397178,0.151))
    R2 = Receiver((64.5832,28.5998,0.151))
    R3 = Receiver((64.5832,28.5998,7.9423))
    R4 = Receiver((-2.40793,31.5003401,0.151))
    R5 = Receiver((75.11005,28.4945787,0.151))

    return 

if __name__ == "__main__" :
    #Only exists for checking bugs with specific parameters now
    #initialize_receivers()
    #ears = Receiver.rList
    #for R in ears:
    #    print(R)
    import RayTrace
    #Receiver.initialize("PointReceivers.txt")
    #print(Receiver.rList)
    #print(Receiver.arraysize)
    pass

#def from input():
#
#ip = 'PointReceivers.txt'
#with open(ip) as vertex:        #Read in from file
#    #print(points)
#    rho = np.genfromtxt(vertex)
#    #rho=np.loadtxt(points)
#Receiver.arraysize = rho.shape[0]   #Create an itterater
##print (Receiver.arraysize)
##print(rho)
#for r in range(Receiver.arraysize): #Translate to usable format
##    print(rho[r,:])
#    #Receiver(rho[r,:])             # Come back later to change the needed inputs for method
#    Receiver(rho[r,0],rho[r,0],rho[r,0])             # Placeholder
#    
#
#K=len(inputsignal)
#print(rho.shape)
#points = (rho.shape)
#print (rho.shape[0])
#for r in points:
#    print(vertex[r,:])


#R1 = Receiver((93.4428213,28.8397178,0.151))
#R2 = Receiver((64.5832,28.5998,0.151))
#R3 = Receiver((64.5832,28.5998,7.9423))
#R4 = Receiver((-2.40793,31.5003401,0.151))
#R5 = Receiver((75.11005,28.4945787,0.151))
#R1.SphereCheck(0.0225,    [-0.94808515,-0.29638514,-0.11528397],    [138.9117595,50.99787506,8.93999318])

#print(R1.recNumber)
#print(R4.recNumber)

#test = Receiver(7,94,3)
#Receiver.frequency = 83
#test.frequency += 7
#
#print(test.frequency)
#print(Receiver.rList[0].frequency)

#print(Receiver.rList)
#print(Receiver.Array)

#Receiver.rList[1].foo = 3
#print(Receiver.rList[1].foo)
#print(Receiver.rList)
#print(Receiver.Array)

## Test List
#Receiver.append_list(R1)
#Receiver.append_list(R2)
#Receiver.append_list(R3)
#Receiver.append_list(R4)
#Receiver.append_list(R5)

#print(R1)
#print(Receiver.rList)              # A list of objects. Probably not helpful in hindsight 
#print(Receiver.rList[0])
#print(Receiver.rList[0].__dict__)   #Show attributes for receiver

#Receiver.rList[0].pressure = 100    #Add pressure directly from list

#print(R1.pressure)                  #Access pressure from specific receiver

#print(Receiver.rList[0].__dict__)   #Show attributes for receiver
# The receiver can now be accessed from two separate places. This may eliminate the need for Receiver.Array completely

## Adding receivers to Array
#Receiver.create_receiverarray()
#Receiver.from_receiver(R1)
#Receiver.from_receiver(R2)
#Receiver.from_receiver(R3)
#Receiver.from_receiver(R4)
#Receiver.from_receiver(R5)

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

#print(R1.position)
#print(R2.position)
#print(Receiver.Array)
#print(Receiver.arraysize)
