#Table of Contents:
#1.    ABSORPTION
#2.    TIMERECONSTRUCT
#3.    RECEIVERHITFUNC
#4.    Header
#5.    TimeHeader
#6.    Grid
#7.    InitialGrid
#8.    SPHERECHECK
#9.    CROSS
#10.   POLYGON (Messy array definitions, I will return later)
#11.   INSIDECHECK (I just realized I skipped this one 
#12.   PLANE (Finished, but should double and triple check for typos)
#13.   BOX

#Notes:
#Eventually get to seeing how much memory it eats
import numpy as np
import math as m
import ReceiverPointSource as RPS 

import time
#Using to check how long functions take
HUGE=1000000.0
XJ=complex(0,1)

def ABSORPTION(ps,freq,hr,Temp):
    ''' This function computes the air absorption for a given frequency, 
    ambient pressure, relative humidity and temperature.
    '''
    #t=time.time()   #Start time counter
# Define all variables and reference values
    ps0=1.0
    hr=20.0
    T0=293.15
    T01=273.16
    F=freq/ps

# Compute all relevant parameters
    psat=ps0*10**(-6.8346*(T01/Temp)**1.261+4.6151)
    h=ps0*(hr/ps)*(psat/ps0)
    FrN=1/ps0*(T0/Temp)**(1/2)*(9+280*h*m.exp(-4.17*((T0/Temp)**(1/3)-1)))
    FrO=1/ps0*(24+4.04*10**4*h*((.02+h)/(.391+h)))
    term1=0.01275*(m.exp((-2239.1/Temp))/(FrO+F**2/FrO))
    term2=0.1068*(m.exp(-3352/Temp)/(FrN+F**2/FrN))
    ABSORPTION=ps0*F**2*((1.84*10**(-11.0)*(Temp/T0)**(0.5)*ps0)+(Temp/T0)**(-5.0/2.0)*(term1+term2))
    return ABSORPTION

# This doesn't get used. But it keeps getting changed when I forget >:{
#def TIMERECONSTRUCT(sizefft,timearray,arraysize,temparray):
#    '''
#    This Function computes the timesignal from a given fft.  It writes the
#    time signal to an array
#    '''
#    print('timeconstruct has been called')
#    #temparray and timetemparray are three dimensional arrays
#    #will create 3d arrays using numpy zeros function
#    timetemparray= np.zeros((arraysize,sizefft,5))
#    # defining tempfft as a 1 dimensional array of size sizefft/2+1
#    tempfft = np.zeros((sizefft//2+1))
#    for D in range(0,arraysize):
#        for W in range(0,sizefft):
#            timetemparray[D,W,0]=temparray[D,0,0]
#            timetemparray[D,W,1]=temparray[D,0,1]
#            timetemparray[D,W,2]=temparray[D,0,2]
#            timetemparray[D,W,3]=timearray[W]
#            timetemparray[D,W,4]=0.0
#    print('timetemparray has been initialized')
#    # Create the complex array to feed into the inverse fft function
#    for D in range(0,arraysize) :
#        if temparray[D,0,4] == 0.0:
#            for W in range(0,sizefft):
#                timetemparray[D,W,4]= 0.0
#        else:
#            for W in range(int(sizefft/2)+1) :
#                if W == 1:
#                    tempfft[W]=complex([0])
#                else:
#                    tempfft[W]=abs(temparray[D,W-1,4])*m.exp(XJ*temparray[D,W-1,5])
#    # use nummpy to compute inverse fft
#    # Original fftw function
#    #             call dfftw_plan_dft_c2r_1d(invplan,sizefft,tempfft,
#    #     *        timesignal, FFTW_ESTIMATE)
#        timesignal=np.fft.irfft(tempfft,sizefft)
#        for W in range(0,sizefft) :
#            timetemparray[D,W,4]=timesignal[W]
#    return timetemparray



def receiverHITFUNC(sizefft,outputarray,arraysize,temparray):
    '''
    This Function adds the pressures from a ray when it hits the receiver.
    '''
    # Define arrays with numpy zeros function
    XJ=complex(0,1)
    #print('everything seems to initiate')
    
    # Add new pressures to existing pressures in temparray 
    # First Look for the correct location.
    for D in range(0,arraysize):
        #print('output: ', outputarray[0,1:3],'\ntemp: ',temparray[D,0,0:2]) #bugfixes
        #if (outputarray[0,1] == temparray[D,0,0] and 
        #outputarray[0,2] == temparray[D,0,1] and 
        #outputarray[0,3] == temparray[D,0,2]):
            #print('first If statement passed')
#        else:
#            print('statement did not pass') #more bug fixes
    # If the location is the same, loop through the frequency and add current values with new values.
        for W in range(0,(sizefft//2)):
            #print(temparray[D,W,5])
            #temp1=complex(np.multiply(abs(temparray[D,W,4]),m.e**(XJ*temparray[D,W,5])))
            temp1=abs(temparray[D,W,4])*m.e**(XJ*temparray[D,W,5])
            #if D==0 and W ==0:
                #print(temp1) 
                #print('abs( ',temparray[D,W,4],')*m.e**( ',XJ,'*',temparray[D,W,5],')')
            #print(temp1)
#            if (W == 0):
# #               print('temp1 fine')
           # temp2=complex(np.multiply(abs(outputarray[W,4]),m.e**(XJ*outputarray[W,5])))
            temp2=abs(outputarray[W,4])*m.e**(XJ*outputarray[W,5])
#            if (W == 0):
  #              print('temp2 fine')
            temp3=temp1+temp2
   #             print('temp3 fine')
            temparray[D,W,4]=abs(temp3) #magnitude
    #            print('temparray 5 fine')
            temparray[D,W,5]=np.arctan2(np.imag(temp3),np.real(temp3))  #direction
    # imagpart and realpart in original code fortran functions
    # using numpy .real and .imag function to acheive same result (Probably)
    #        if (W == 1):
    #            print('temparray 6 fine')
    #            print(temparray[1,W,5])
            #print(temparray[D,W,4:5])
    #print(temp2,temp1)
    #print('final temparray: ',temparray[1,0,4], ' and ', temparray[1,0,5])
    #print('Got Through the end')
    #print(temparray[1,1,4])
    return temparray

def Header(outputfile):
    '''
    This function prints the header for the tecplot data

    '''
    #this function prints the header for the tecplot data
    f=open(outputfile,"w")

    f.write('TITLE = "Pressure at earlevel"\n' )
    f.write('VARIABLES = "X[m]" "Y[m]" "Z[m]" "P[Pa]"\n' )
    f.write('TEXT\n' )
    f.write('CS=FRAME\n' )
    f.write('X=71.9660948264,Y=82.9866270431\n' )
    f.write('C=BLACK\n' )
    f.write('S=LOCAL\n' )
    f.write('HU=POINT\n' )
    f.write('LS=1 AN=MIDCENTER\n' )
    f.write('BX=Filled BXM=60 LT=0.1 BXO=BLACK BXF=WHITE\n' )
    f.write('F=HELV\n' )
    f.write('H=20 A=0\n' )
    f.write('MFC=""\n' )
    f.write('CLIPPING=CLIPTOVIEWPORT\n' )
    f.write('T="Time = &(SOLUTIONTIME%4f)" \n'  )
    Header=0
    #f.close()
    return Header


def TimeHeader(f,time,sizex,sizey,sizez,planename):
    """
    This function prints the header between each time set
    """
    #This function prints the header between each time set
    f.write('ZONE T=" %s "\n' %(planename, ) )  #this worked in the command line
    f.write('STRANDID=1, SOLUTIONTIME= %d \n'%(time,) )
    f.write('I= %d\t J= %d\t K=%d\t ZONETYPE=Ordered\n' %(sizex, sizey, sizez))
    f.write('DATAPACKING=POINT\n')
    f.write('DT=(SINGLE SINGLE SINGLE SINGLE )\n')
    header=0
    return TimeHeader


#*****Grid is Unused*****************    
#************************************
#def Grid(radius,A,B,C,D,xmin,ymin,zmin,xmax,ymax,zmax,receiverarray, arraysize,step):
#    # this function creates an equally spaced grid of size step apart for a receiver plane array
#    # following line defines the array, unsure what ,s outside of parentheses does
#    #real receiverarray(arraysize,3),s
#    
#    receiverarray=np.zeros(arraysize,3)
#    s=step/radius
#    if xmin == xmax:
#        count = 0
#        while i < int((zmax-zmin)/step):
#            while j < int((ymax-ymin)/(step)):
#                receiverarray[count,0]=(D-B*ymin+(s*j+(1-s))*radius)-C*(zmin+(s*i+(1-s))*radius)/A
#                receiverarray[count,1]=ymin+(s*j+(1-s))*radius
#                receiverarray[count,2]=zmin+(s*i+(1-s))*radius
#                count=count+1
#                j += 1
#            i += 1
#        sizex=int((ymax-ymin)/(step))
#        sizey=int((zmax-zmin)/step)
#        sizez=1
#    if ymin == ymax:
#        count = 0
#        while i < int((xmax-xmin)/(step)):
#            while j < int((zmax-zmin)/(step)):
#                receiverarray[count,0]=xmin+(s*i+(1-s))*radius
#                receiverarray[count,1]=(D-A*(xmin+(s*i+(1-s))*radius)-C*(zmin+(s*j+(1-s))*radius))/B
#                receiverarray[count,2]=zmin+(s*j+(1-s))*radius
#                count=count+1
#                j += 1
#            i += 1
#        sizex=int((zmax-zmin)/step)
#        sizey=int((xmax-xmin)/(step))
#        sizez=1
#    if zmin == zmax:
#        count = 0
#        while i < int((xmax-xmin)/(step)):
#            while j < int((ymax-ymin)/(step)): 
##        for i in range(int((xmax-xmin)/(step))):
##            for j in range(int((ymax-ymin)/(step))):
#                receiverarray[count,1]=xmin+(s*i+(1-s))*radius
#                receiverarray[count,2]=ymin + (s*j+(1-s))*radius-B*(ymin+(s*j+(1-s))*radius)/C
#                count=count+1
#                j += 1
#            i += 1
#        sizex=int((ymax-ymin)/(step))
#        sizey=int((xmax-xmin)/step)
#        sizez=1
#    return receiverarray,sizex,sizey,sizez
#********************************************
    
def InitialGrid(radius,A,B,C,D,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,arraysize):
    """This function creates an equally spaced grid of size "step" apart"""
    receiverarray=np.zeros((arraysize,3))
    yspace=radius*abs(m.cos(phi))
    zspace=radius*abs(m.sin(theta))
    count = 0
    i=0
    j=0
    #this works for bug fixes. But it's not that great.
    if xmin == xmax:
        for i in range(0,int((zmax-zmin)//zspace)):
            for j in range(0,int((ymax-ymin)//yspace)):
                receiverarray[count,0]=(D-B*(ymin+(j+1)*yspace)-C*(zmin+(i+1)*zspace))/A
                receiverarray[count,1]=ymin+(j+1)*yspace
                receiverarray[count,2]=zmin+(i+1)*zspace
                count=count+1
        # // by itself does not convert to an int
        sizex=int((ymax-ymin)//(yspace))
        sizey=int((zmax-zmin)//zspace)
        sizez=1
    if ymin == ymax:
        for i in range(0,int(xmax-xmin)//(xspace)):
            for j in range(0,int((zmax-zmin)//(zspace))):
                receiverarray[count,0]=xmin+(i+1)*xspace
                receiverarray[count,1]=(D-A*(xmin+(i+1)*xspace)-C*(zmin+(j+1)*zspace))/B
                receiverarray[count,2]=zmin+j*zspace
                count=count+1
        sizex=int((zmax-zmin)/zspace)
        sizey=int((xmax-xmin)/(xspace))
        sizez=1
    if zmin == zmax:
        for i in range(int((xmax-xmin)/(xspace))):
            for j in range(int((ymax-ymin)/(yspace))):
                receiverarray[count,0]=xmin+(i+1)*xspace  
                receiverarray[count,1]=ymin+(j+1)*yspace
                receiverarray[count,2]=(D-A*(xmin+(i+1)*xspace)-B*(ymin+(j+1)*yspace))/C
                count=count+1
        sizex=int((xmax-xmin)/(xspace))
        sizey=int((ymax-ymin)/yspace)
        sizez=1
    return receiverarray, sizex, sizey, sizez

def SPHERECHECK(Sc,Sr2,F,veci):
    '''
    This function performs a check whether a ray hits a sphere.  If
    it does hit the function returns the distance to the sphere
    '''
    HUGE=1000000.0
    OC=np.zeros(3)      #put a pin in this

    OC[0]=Sc[0]-veci[0]
    OC[1]=Sc[1]-veci[1]
    OC[2]=Sc[2]-veci[2]
    L2OC=np.dot(OC,OC)
    tca=np.dot(OC,F)
    #Takes dot product of OC and OC (Dot square?)
    t2hc=Sr2-L2OC+tca**2
    if L2OC == Sr2:
        dx = HUGE
    elif tca < 0.0:
        dx = HUGE
    elif t2hc < 0.0:
        dx = HUGE
    else:
        dx = tca - (t2hc**(1/2))
    return dx

def CROSS(A, B):
    '''
    This function calculates a cross product of A and B and returns normal
    '''
    normal=np.zeros(3)
    normal[0]=A[1]*B[2]-A[2]*B[1]
    normal[1]=A[2]*B[0]-A[0]*B[2]
    normal[2]=A[0]*B[1]-A[1]*B[0]
    length=( normal[0]**2.0 + normal[1]**2 +normal[2]**2)**(1/2)
    if (length != 0.0):
        normal=normal/length
    #print('normal is ',normal)
    return normal

#def tri(veci,F,Q,Number,PointNumbers,PolyArray,BuildingPoints,normal,FaceNormalNo,vn,dxbuilding,behind):
#def tri(veci,F,Q,Number,PointNumbers,PolyArray,vertices,normal,FaceNormalNo,vn,dxbuilding,behind):
def tri(veci,F,Q,Number,PointNumbers,PolyArray,v,normal,FaceNormalNo,vn,dxbuilding,behind):
    """
    Lab notebook 5/16
    This is an attempt to merge box and polygon into one function since
    we are working entirely in triangular meshes now

    ********************************Untested***********************************
    A 1:1 translation was made from Fortran. This is the closest match to
    what we are trying to do with triangle geometry. However there is no 
    readily available geometry file to test this. 
    [No Description given in Fortran]
    """
    size = 3
    G = np.zeros((size,2))  # (3,2)
    # inits
    NC=0
    behind=0
    normal = vn[PolyArray[Q,1],:]
    #normal[0]=FaceNormals[int(PolyArray[Q,1]),0]
    #normal[1]=FaceNormals[int(PolyArray[Q,1]),1]
    #normal[2]=FaceNormals[int(PolyArray[Q,1]),2]

    d=-np.dot(normal,ValueError(PolyArray[Q,2]))
    Vd=np.dot(normal,F)

    if Vd >= 0.0:
        dxbuilding = HUGE
    V0= -(np.dot(normal,veci)+d)
    t=V0/Vd
    if(t < 0.0):
        dxbuilding=HUGE
        behind = 1
        #Stage 1
    intersection = veci + F*t
    maximum = max(abs(normal))
        # G: What if two normal values are the same? Anyway:
    if(maximum == abs(normal[0])):
        for P in range(size):
            G[P,:] = (intersection[1]-v[int(PolyArray[Q,1+P]),1]
                      ,intersection[2]-v[int(PolyArray[Q,1+P]),2])
    elif (maximum == normal[1]):
        for P in range(size):
            G[P,:] = (intersection[0]-v[int(PolyArray[Q,1+P]),0]
                      ,intersection[2]-v[int(PolyArray[Q,1+P]),2])
    elif (maximum == normal[2]):
        for P in range(size):
            G[P,:] = (intersection[0]-v[int(PolyArray[Q,1+P]),0]
                      ,intersection[1]-v[int(PolyArray[Q,1+P]),1])
    #Stage 2
    for P in range(size):
        if P == size:
            if G[P,1] < 0.0:
                SH = -1
            else:
                SH = 1
            if G[0,1] < 0.0:
                NSH = -1
            else:
                NSH = 1 
        else:
            if G[P,1] < 0.0:
                SH = -1
            else:
                SH = 1
            if G[P+1,2] < 0.0:
                NSH = -1
            else:
                NSH = 1
        if SH != NSH:
            if (P == size):
                if (G[P,0] > 0.0) and (G[0,0]>0.0):
                    NC += 1
                elif (G[P,0]> 0.0) or (G[0,0] > 0.0):
                    if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
                        NC += 1
            else:
                if (G[P,0] > 0.0) and (G[P+1,0] > 0.0):
                    NC += 1
                elif (G[P,0] > 0.0) or (G[P+1,1] > 0.0):
                    if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
                        NC += 1
        odd = NC % 2    #get remainder to find if odd or not
        # This was this way in original fortran
        if odd:
            dxbuilding = t
        else:
            dxbuilding = HUGE

        return dxbuilding,behind

def POLYGON(Vecip1,F,Q,size,Number,PointNumbers,PolyArray,BuildingPoints,normal,FaceNormalNo,FaceNormals,dxbuilding):
    '''
    ********************************Untested***********************************
    A 1:1 translation was made from Fortran. This is the closest match to
    what we are trying to do with triangle geometry. However there is no 
    readily available geometry file to test this.
    
    [No Description given in Fortran]
    '''

    G = np.zeros((size,2))
    HUGE=1000000.0
    NC=0
    behind=0
    normal[1]=FaceNormals[int(PolyArray[Q,1]),0]
    normal[2]=FaceNormals[int(PolyArray[Q,1]),1]
    normal[3]=FaceNormals[int(PolyArray[Q,1]),2]
    #An array defined as an array from a function of an array and a point.
    #I will recheck syntax, just getting through everything now

    # This is what Kory originally had written, returning syntax error-Will
    # d=-np.dot(normal,BuildingPoints(int(PolyArray[Q,2]),1:3))
    #^^^unsure how to translate this part

    # Will's attempt (To avoid error):
    d=-np.dot(normal,BuildingPoints(PolyArray[Q,2]))
    # it ran, doesn't mean it's right. Will check back.
    Vd=np.dot(normal,F)

    if Vd >= 0.0:
        dxbuilding = HUGE
    V0= -(np.dot(normal,Vecip1)+d)
    t=V0/Vd
    if(t < 0.0):
        dxbuilding=HUGE
        behind = 1
    
    intersection = Vecip1 + F*t
    #intersection[1]=Vecip1[1]+F[1]*t
    #intersection[2]=Vecip1[2]+F[2]*t
    #intersection[3]=Vecip1[3]+F[3]*t
    maximum = max(abs(normal))
    if(maximum == abs(normal[0])):
        for P in range(size):
            G[P,:2] = (intersection[1]-BuildingPoints[int(PolyArray[Q,1+P]),1]
                      ,intersection[2]-BuildingPoints[int(PolyArray[Q,1+P]),2])
    elif (maximum == normal[1]):
        for P in range(size):
            G[P,:2] = (intersection[0]-BuildingPoints[int(PolyArray[Q,1+P]),0]
                      ,intersection[2]-BuildingPoints[int(PolyArray[Q,1+P]),2])
    elif (maximum == normal[2]):
        for P in range(size):
            G[P,:2] = (intersection[0]-BuildingPoints[int(PolyArray[Q,1+P]),0]
                      ,intersection[1]-BuildingPoints[int(PolyArray[Q,1+P]),1])
    for P in range(size):
        if P == size:
            if G[P,1] < 0.0:
                SH = -1
            else:
                SH = 1
            if G[0,1] < 0.0:
                NSH = -1
            else:
                NSH = 1 
        else:
            if G[P,1] < 0.0:
                SH = -1
            else:
                SH = 1
            if G[P+1,2] < 0.0:
                NSH = -1
            else:
                NSH = 1
        if SH != NSH:
            if (P == size):
                if (G[P,0] > 0.0) and (G[0,0]>0.0):
                    NC += 1
                elif (G[P,0]> 0.0) or (G[0,0] > 0.0):
                    if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
                        NC += 1
            else:
                if (G[P,0] > 0.0) and (G[P+1,0] > 0.0):
                    NC += 1
                elif (G[P,0] > 0.0) or (G[P+1,1] > 0.0):
                    if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
                        NC += 1
        odd = NC % 2    #get remainder to find if odd or not
        # This was this way in original fortran
        #if odd== 0:     
        #    dxbuilding = HUGE
        #else:
        #    dxbuilding = t
        # Setting it this way myself:
        if odd:
            dxbuilding = t
        else:
            dxbuilding = HUGE
        return dxbuilding,behind

def PLANE(Vecip1, B1, B2, planehit):
    '''
    This function calculates the normal at the hitpoint of a box.
    '''
#George:It sure would be a mess if there was a typo anywhere in here
    #This function calculates the normal at the hitpoint of a box.
    #import numpy as np
    #print('planehit is ', planehit)
    #global nbox
    nbox=[0,0,0]
    if planehit == 1:
        #print('vecip1',Vecip1)
        #print('B1 is ',B1)
        if Vecip1[0] == B1[0]:
            Point2=[B1[0],B1[1],B2[2]]  
            Point3=[B1[0],B2[1],B1[2]] 
            nbox=CROSS(np.subtract(Point2,B1),np.subtract(Point3,B1))

        elif (Vecip1[0] == B2[0]) :
            Point1=(B2[0],B1[1],B1[2])
            Point2=(B2[0],B1[1],B2[2])
            Point3=(B2[0],B2[1],B1[2])
            #print('points: ',Point1,Point2,Point3)
            #nbox=CROSS((Point3-Point1),(Point2-Point1))
            nbox=CROSS(np.subtract(Point3,Point1),np.subtract(Point2,Point1))
            #print('nbox works. It is', nbox)
    if planehit == 2:
        #print('this happens 2')

        #print(Vecip1[1],B1[1])
 #****************************************************************
                    #This is not a good solution
        if Vecip1[1] == B1[1]:
            Point2=(B2[0], B1[1], B1[2])  
            Point3=(B1[0], B1[1], B2[2]) 
            nbox=CROSS(np.subtract(Point2,B1),np.subtract(Point3,B1))
            #print('test one')
        elif Vecip1[1] == B2[1]: 
            #is this really correct???
            #print('test two')
            Point1=(B1[0],B2[1],B1[0])
            Point2=(B1[0],B2[1],B2[2])
            Point3=(B2[0],B2[1],B1[2])
            nbox=CROSS(np.subtract(Point2,Point1),np.subtract(Point3,Point1))
    if planehit == 3:
         if Vecip1[2] == B1[2]:
            Point2=(B2[0], B1[1], B1[2])  
            Point3=(B1[0], B2[1], B1[2]) 
            nbox=CROSS(np.subtract(Point3,B1),np.subtract(Point2,B1))
         elif Vecip1[2] == B2[2]:
            Point2=(B1[0],B2[1],B2[2])
            Point3=(B2[0],B1[1],B2[2])
            nbox=CROSS(np.subtract(Point2,B2),np.subtract(Point3,B2))
    #print('Here is that nbox it keeps saying you\'re missing', nbox )
    return nbox
    #return

# Call: for Q in range(0,BG.Boxnumber):
#           dxnear, dxfar, hit, planehit=fun.BOX(BG.Boxarraynear[Q], BG.Boxarrayfar[Q],veci,F)

def BOX(B1,B2,Vecip1,F):
    '''
    This function checks to see if the ray hits a box.  It determines which
    plane the ray hits
        T1x is the distance to the close side
        T2x is th distance to the far side
    '''
    hit=5
    HUGE=1000000.0
    dxnear=-HUGE        
    dxfar=HUGE
    tempF=F
    planehit=0
    if ((F[0] == 0.0) or (F[1] == 0.0) or (F[2] == 0.0)):
        if (F[0] == 0.0):
            if((vecip1[0] < B1[0]) or (vecip1[0] > B2[0])):
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
        if (F[1] == 0.0):
            if((vecip1[1] < B1[1]) or (vecip1[1] > B2[1])):
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
        if (F[2] == 0.0):
            if((vecip1[2] < B1[2]) or (vecip1[2] > B2[2])):
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
#    print('hit',hit)
    if hit != 0.0 :

        if F[0] == 0.0:
            tempF[0]=1.0
        if F[1] == 0.0:
            tempF[1]=1.0
        if F[2] == 0.0:
            tempF[2]=1.0
#        print(B1[1],B2[1], Vecip1[1], F[1])
        if (F[0] != 0.0):
            T1X=(B1[0]-Vecip1[0])/tempF[0]
            T2X=(B2[0]-Vecip1[0])/tempF[0]
#            print('does this ever even happen')     #yes
#        print('T1X,T2x',T1X,T2X)
#        print(T1X, T2X)
            if T1X > T2X :
                tmp =T1X
                T1X =T2X
                T2X =tmp
#                print('how about now')      #nope
            if T1X > dxnear:
#                print('now? ********')      #yup
                dxnear = T1X
            if T2X < dxfar:
                dxfar = T2X
            if dxnear > dxfar :
                hit = 0
                dxnear = HUGE
                return dxnear, dxfar,hit, planehit
            elif dxfar < 0.0:
                hit =0
                dxnear= HUGE
                return dxnear, dxfar,hit, planehit
        if F[1] != 0.0 :
            T1Y=(B1[1]-Vecip1[1])/tempF[1]
            T2Y=(B2[1]-Vecip1[1])/tempF[1]
#            print('T1Z, T2Z',T1Y,T2Y)
            if T1Y > T2Y:
               tmp = T1Y
               T1Y = T2Y
               T2Y = tmp
            if T1Y > dxnear : 
                dxnear=T1Y
            if T2Y < dxfar : 
                dxfar=T2Y
            if dxnear > dxfar :
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
            elif (dxfar < 0.0):
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit

#"break outside loop"
#Look into it later -G

                #goto 100
        #print(B1(3), Vecip1(3), tempF(3))
        #print(B2(3), Vecip1(3), tempF(3))
        if F[2] != 0.0:
            T1Z=(B1[2]-Vecip1[2])/tempF[2]
            T2Z=(B2[2]-Vecip1[2])/tempF[2]
#            print('T1Z, T2Z',T1Z,T2Z)
            if T1Z > T2Z:
                tmp=T1Z
                T1Z=T2Z
                T2Z=tmp
            if T1Z>dxnear:
                dxnear=T1Z
            if T2Z < dxfar:
                dxfar=T2Z
            if dxnear > dxfar:
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
            elif dxfar < 0:
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
            elif dxnear < 0:
                hit=0
                dxnear=HUGE
                return dxnear, dxfar,hit, planehit
#                break
#            break
#        print('dxnear ',dxnear,'dxfar ', dxfar, ' and hit is ', hit)
    if hit != 0:
        if dxnear < dxfar:
            hit =1
#            print(T1X, T1Y,T1Z, ' is this right?')
            if dxnear == T1X:
                planehit = 1
            if dxnear == T1Y:
                planehit = 2
            if dxnear == T1Z:
                planehit = 3
                #planehit is 0,1,2,3 We don't need to change the indices        10/2/18
    return dxnear, dxfar,hit, planehit

# William Costa: Function 'BOX 2' in Functions.f appears to be another version of 'BOX'
# 'BOX2' does not appear to be called in any instance in RayTrace.f
# Will pass on porting for now, will port if needed

def ROTATION(axis, angle, rotationmatrix):
    '''
    [No Description given in Fortran]
    '''
    #axis = np.zeros(3)
    #rotationmatrix= np.zeros(3,3)
    rotationmatrix[1,1]= axis[1]**2+(1-axis[1]**2)*m.cos(angle)
    rotationmatrix[1,2]=axis[1]*axis[2]*(1-cos(angle))+axis[3]*m.sin(angle)
    rotationmatrix[1,3]=axis[1]*axis[3]*(1-cos(angle))-axis[2]*m.sin(angle)                
    rotationmatrix[2,1]=axis[1]*axis[2]*(1-cos(angle))-axis[3]*m.sin(angle)                 
    rotationmatrix[2,2]=axis[2]**2+(1-axis[2]**2)*m.cos(angle)              
    rotationmatrix[2,3]=axis[2]*axis[3]*(1-cos(angle))+axis[1]*m.sin(angle)
    rotationmatrix[3,1]=axis[1]*axis[3]*(1-cos(angle))+axis[2]*m.sin(angle)
    rotationmatrix[3,2]=axis[2]*axis[3]*(1-cos(angle))-axis[1]*m.sin(angle)
    rotationmatrix[3,3]=axis[3]**2+(1-axis[3]**2)*m.cos(angle)
    return rotationmatrix

#print(__name__ ,"__main__")
#We finished!?
#ARE YOU NOT ENTERTAINED?! -G