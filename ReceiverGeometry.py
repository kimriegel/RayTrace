# Receiver Point 
# Python 3.7.0 

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This is a module that converts all points of a given geometry file to 
receivers while maintaining their functionality as a building. For 
using singular points please refer to Receiver Point Source.

Will function as plane receiver and be a subclass of Geometry, but 
as that has not been integrated yet, this is here
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


import numpy as np
import pywavefront as pwf

class Receiver:

    planename1='Single Point'
    arraysize=0     #Number of receivers. Supposed to be 5 for this test    # I use it differently now, and only in one function
    sizex=2
    sizey=2
    sizez=1
    initial_frequency = None    # Gives this value to all receivers

    rList = [] #See append_list

    def __init__(self,position):
        """
        Create and defines position of receiver
        
        Works automatically when class is called
        """
        pos = (position[0],position[2],position[1])
        self.position=np.array(pos)
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
            self.timesignal = np.zeros(sizefft) 

        else:                                   # If not then calculate the timesignal 
            tempfft = abs(self.magnitude[:]) * np.exp(XJ*self.direction)
            tempfft = np.append(0,tempfft)
            # Compute ifft using numpy
            self.timesignal=np.fft.ifft(tempfft,sizefft)

    @classmethod
    def initialize(cls,ipname):
        """
        Reads in receiver points from the txt file and translates them to be used in our receiver method
        Receivers are automatically initialized from given inputfile
        """
        ipfile = pwf.Wavefront(ipname)    # Read in geometry file 
        env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
                create_materials=True, collect_faces=True, parse=True, cache=False)
        vertices = env.wavefront.vertices   # Sample the vertices
        for v in vertices:
            Receiver(v)        
            print('initialized receivers')

if __name__ == "__main__" :
    ears = Receiver.rList
    Receiver.initialize('SingleBuilding.obj')
    for r in ears:
        print(r.position)
    # Turning a geometry file into a receiver cluster   -Works 5/16/19

    # For backwards compatibility
if __name__ != "__main__":      # only runs if opened outside this file
    planenum = 1
    planename1 = 'Single Point'
    arraysize1 = Receiver.arraysize
    sizex = 2
    sizey = 2
    sizez = 1
    sizex1 = sizex
    sizey1 = sizey
    sizez1 = sizez
