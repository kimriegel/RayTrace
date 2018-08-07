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


def ABSORPTION(ps,freq,hr,Temp):
# This function computes the air absorption for a given frequency, 
# ambient pressure, relative humidity and temperature.
    import math as m
# Define all variables and reference values
#      real ps0, ps, freq, hr, Temp, T0, T01, F, FrN, FrO
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


#******************Unfinished*************

def TIMERECONSTRUCT(sizefft,timearray,arraysize,temparray,timetemparray):

    import numpy as np
    XJ=(0,1)

    print('timeconstruct has been called')

    #temparray and timetemparray are three dimensional arrays
    #will create 3d arrays using numpy zeros function

    temparray = np.zeros((arraysize,sizefft/2,6))
    timetemparray= np.zeros((arraysize,sizefft,5))

    # Author: Will
    # For loop iterating through three dimensional arrays
    for D in range(1,arraysize):
        for W in range(1,sizefft):
            timetemparray[D,W,1]=temparray[D,1,1]
            timetemparray[D,W,2]=temparray[D,1,2]
            timetemparray[D,W,3]=temparray[D,1,3]
            timetemparray[D,W,4]=timearray(W)
            timetemparray[D,W,5]=0.0
    print('timetemparray has been initialized')

    # Create the complex array to feed into the inverse fft function

    #original fortran code

        #C Create the complex array to feed into the inverse fft function
        #      DO 37 D=1, arraysize
        #         if (temparray(D,1,5).eq.0.0)then
        #            DO 40 W=1,sizefft
        #               timetemparray(D,W,5)=0.0
        # 40         CONTINUE
        #         else
        #         DO 38 W=1,sizefft/2+1
        #            if (W.eq.1) then
        #               tempfft(W)=cmplx(0.0)
        #            else
        #            tempfft(W)=cmplx(abs(temparray(D,W-1,5))*exp(XJ*
        #     *           temparray(D,W-1,6)))
        #            endif
        # 38      CONTINUE
        #         print*, 'created temparray'

    #Up to this point converted to Python
    # Find equivalent fft functions in Numpy

        #C use fftw to compute the inverse fft.
        #         call dfftw_plan_dft_c2r_1d(invplan,sizefft,tempfft,
        #     *        timesignal, FFTW_ESTIMATE)
        #         call dfftw_execute(invplan, tempfft, timesignal)
        #         call dfftw_destroy_plan(invplan)
        #         print*, 'created time signature'
        #         DO 39 W=1,sizefft         
        #            timetemparray(D,W,5)=timesignal(W)
        # 39      CONTINUE
        #         endif
        # 37   CONTINUE
        #      return
        #      end

    # Author: Will, Create complex array and compue inverse fft first attempt Python

    for D in range(1,arraysize):
        if temparray[D,1,5]== 0.0:
            for W in range(1,sizefft):
                timetemparray[D,W,5]= 0.0
        else:
            for W in range(1,(sizefft/2)+1):
                if W == 1:
                    tempfft(W)=complex(0.0)
                else:
                    tempfft(W)=complex(abs(temparray[D,W-1,5])*m.exp(XJ*temparray[D,W-1,6]))
        print('Created temparray')

    # use nummpy to compute inverse fft
    # use ifft numpy function with tempfft and sizefft as input
    # use timesignal as output

    # Original fftw function
    #             call dfftw_plan_dft_c2r_1d(invplan,sizefft,tempfft,
    #     *        timesignal, FFTW_ESTIMATE)

        timesignal=np.ifft(tempfft,sizefft)
        print('Created time signature')       
        for W in range(1,sizefft):
            timetemparray[D,W,5]=timesignal(W)
    return TIMERECONSTRUCT




def RECEIVERHITFUNC(sizefft,outputarray,arraysize,temparray):
    ##All print commands commented out were from original code but stayed for consistency
    import math as 
    
    import numpy as np
    ##This Function adds the pressures from a ray when it hits the receiver.

    ## Define all variables
   
    # original fortran variables

    #      INTEGER sizefft,D,arraysize,W
    #      real outputarray(sizefft/2,7)
    #      real temparray(arraysize,sizefft/2,6)
    #      double complex temp1, temp2,temp3
    #      COMPLEX XJ

    # Define arrays with numpy zeros function

    outputarray=np.zeros(sizefft/2,7)
    temparray=np.zeros(arraysize,sizefft/2,6)
    XJ=(0,1)
    print('everything seems to initiate')
    
    # Add new pressures to existing pressures in temparray 
    # First Look for the correct location.

    for D in range(1,arraysize):
        print('outputarray2', outputarray[1,2])
        print('outputarray3', outputarray[1,3])
        print('outputarray4', outputarray[1,4])
        print('temparray1',temparray[D,1,1])
        print('temparray2',temparray[D,1,2])
        print('temparray3',temparray[D,1,3])
        if (outputarray[1,2] == temparray[D,1,1] and outputarray[1,3] == temparray[D,1,2] and outputarray[1,4] == temparray[D,1,3]):
            print('first If statement passed')
    # If the location is the same, loop through the frequency and add current values with new values.

        for W in range(1,sizefft/2):
            temp1=complex(abs(temparray[D,W,5])*m.exp(XJ*temparray[D,W,6]))
            if (W == 1):
                print('temp1 fine')
                temp2=cmplx(abs(outputarray[W,5])*m.exp(XJ*outputarray[W,6]))
            if (W == 1):
                print('temp2 fine')
                temp3=temp1+temp2
            if (W == 1):
                print('temp3 fine')
                temparray[D,W,5]=abs(temp3)
            if (W == 1):
                print('temparray 5 fine')
                temparray[D,W,6]=ATAN2(np.imag(temp3),np.real(temp3))
    # imagpart and realpart in original code fortran functions
    # using numpy .real and .imag function to acheive same result (Probably)

            if (W == 1):
                print('temparray 6 fine')
                print(temparray[1,W,5])
    # combine if statements?
    print('Got Through the end')
    return RECEIVERHITFUNC

def Header(fileid):
    #this function prints the header for the tecplot data
    f=open(fileid,"w")

    f.write('TITLE = "Pressure at earlevel"/n' )
    f.write('VARIABLES = "X[m]" "Y[m]" "Z[m]" "P[Pa]"/n' )
    f.write('TEXT/n' )
    f.write('CS=FRAME/n' )
    f.write('X=71.9660948264,Y=82.9866270431/n' )
    f.write('C=BLACK/n' )
    f.write('S=LOCAL/n' )
    f.write('HU=POINT/n' )
    f.write('LS=1 AN=MIDCENTER/n' )
    f.write('BX=Filled BXM=60 LT=0.1 BXO=BLACK BXF=WHITE/n' )
    f.write('F=HELV/n' )
    f.write('H=20 A=0/n' )
    f.write('MFC=""/n' )
    f.write('CLIPPING=CLIPTOVIEWPORT/n' )
    f.write('T="Time = &(SOLUTIONTIME%4f)" /n'  )
    f.close(fileid)
    Header=0
    return Header


def TimeHeader(fileid,time,sizex,sizey,sizez,planename):
    #This function prints the header between each time set
    f=open(fileid,"w")
    f.write('ZONE',' T="',planename,'"/n')
    f.write('STRANDID=1, SOLUTIONTIME=','time, /n' )
    f.write('I=',sizex,'J=',sizey,'K=',sizez,'ZONETYPE=Ordered/n')
    f.write('DATAPACKING=POINT/n')
    f.write('DT=(SINGLE SINGLE SINGLE SINGLE )/n')
    f.close(fileid)
    header=0
    return TimeHeader

def Grid(radius,A,B,C,D,xmin,ymin,zmin,xmax,ymax,zmax,receiverarray, arraysize,sizex,sizey,sizez,step):
    # this function creates an equally spaced grid of size step apart for a receiver plane array

    # following line defines the array, unsure what ,s outside of parentheses does
    #real receiverarray(arraysize,3),s

    np.zeros(arraysize,3)
    s=step/radius
    if xmin == xmax:
        count = 1
        for i in range(1,int((zmax-zmin)/step)):
            for j in range(1,innt((ymax-ymin)/(step))
                receiverarray(count,1)=(D-B*ymin+(s*j+(1-s))*radius)-C*(zmin+(s*i+(1-s))*radius))/A
                receiverarray(count,2)=ymin+(s*j+(1-s))*radius
                receiverarray(count,3)=zmin+(s*i+(1-s))*radius
                count=count+1
        sizex=int((ymax-ymin)/(step))
        sizey=int((zmax-zmin)/step)
        sizez=1
    if ymin == ymax:
        count = 1
        for i in range(1,int((xmax-xmin)/(step))):
            for j in range(1,int((zmax-zmin)/(step))):
                receiverarray(count,1)=xmin+(s*i+(1-s))*radius
                receiverarray(count,2)=(D-A*(xmin+(s*i+(1-s))*radius)-C*(zmin+(s*j+(1-s))*radius))/B
                receiverarray*count,3)=zmin+(s*j+(1-s))*radius
                count=count+1
        sizex=int((zmax-zmin)/step)
        sizey=int((xmax-xmin)/(step))
        sizez=1
    if zmin == zmax:
        count = 1
        for i in range(int((xmax-xmin)/(step))):
            for j in range(int((ymax-ymin)/(step)):
                receiverarray(count,1)=xmin+(s*i+(1-s))*radius
                receiverarray*count,2)=ymin + (s*j+(1-s))*radius)-B*(ymin+(s*j+(1-s))*radius))/C
                count=count+1
        sizex=int((ymax-ymin)/(step))
        sizey=int((xmax-xmin)/step)
        sizez=1
    return Grid
    
def InitialGrid(radius,A,B,C,D,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,receiverarray, arraysize,sizex,sizey,sizez):
    #This function creates an equally spaced grid of size step apart
    yspace=radius*abs(cos(phi))
    zspace=radius*abs(sin(theta))
    if xmin == xmax:
        count=1      
#        DO 1 i=1,int((zmax-zmin)/zspace),1
#            Do 2 j=1,int((ymax-ymin)/(yspace)),1
#does the final 1 mean that it counts in increments of one?
        for i in range(1,int((zmax-zmin)/zpace)):
            for j in range(1,int((ymax-ymin/yspace))):
                receiverarray[count,1]=(D-B*(ymin+j*yspace)-C*(zmin+i*zspace))/A
                receiverarray[count,2]=ymin+j*yspace
                receiverarray[count,3]=zmin+i*zspace
                count=count+1
        sizex=int((ymax-ymin)/(yspace))
        sizey=int((zmax-zmin)/zspace)
        sizez=1
    if ymin == ymax:
        count=1
        for i in range(1,int(xmax-xmin)/(xspace)):
            for j in range(1,int(zmax-zmin)/(zspace)):
                receiverarray[count,1]=xmin+i*xspace
                receiverarray[count,2]=(D-A*(xmin+i*xspace)-C*(zmin+j*zspace))/B
                receiverarray[count,3]=zmin+j*zspace
                count=count+1
        sizex=int((zmax-zmin)/zspace)
        sizey=int((xmax-xmin)/(xspace))
        sizez=1
    if zmin == zmax:
        count = 1
        for i in range(1,int(xmax-xmin)/(xspace)):
            for j in range(1,int(ymax-ymin)/(yspace)):
                receiverarray[count,1]=xmin+i*xspace  
                receiverarray[count,2]=ymin+j*yspace
                receiverarray[count,3]=(D-A*(xmin+i*xspace)-B*(ymin+j*yspace))/C
                count=count+1
        sizex=int((xmax-xmin)/(xspace))
        sizey=int((ymax-ymin)/yspace)
        sizez=1
    return InitialGrid



def SPHERECHECK(Sc,Sr2,F,veci1):
    import numpy as np
    #This function performs a check whether a ray hits a sphere.  If
    #it does hit the function returns the distance to the sphere
    HUGE=1000000.0
    OC[1]=Sc[1]-veci[1]
    OC[2]=Sc[2]-veci[2]
    OC[3]=Sc[3]-veci[3]
    L2OC=np.dot(OC,OC)
    tca=np.dot(OC,F)
    #Takes dot product of OC and OC (Dot square?)
    #F may be de[f]ined elsewhere, but [f]inding it is going to be a mess.
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
    #This function calculates a cross product of A and B and returns normal

    normal[1]=A[2]*B[3]-A[3]*B[2]
    normal[2]=A[3]*B[1]-A[1]*B[3]
    normal[3]=A[1]*B[2]-A[2]*B[1]
    length=(normal[1]**2.0+normal[2]**2+normal[3]**2)**(1/2)
    if (length != 0.0):
        normal=normal/length
    return normal

def POLYGON(Vecip1,F,Q,size,Number,PointNumbers,PolyArray,BuildingPoints,normal,FaceNormalNo,FaceNormals,dxbuilding,behind):
    #********************************Unfinished***********************************

    import numpy as np
    HUGE=1000000.0
    NC=0
    behind=0
    normal[1]=FaceNormals[int(PolyArray[Q,1]),1]
    normal[2]=FaceNormals[int(PolyArray[Q,1]),2]
    normal[3]=FaceNormals[int(PolyArray[Q,1]),3]
    #An array defined as an array from a function of an array and a point.
    #I will recheck syntax, just getting through everything now

    d=-np.dot(normal,BuildingPoints(int(PolyArray(Q,2)),1:3))
    #^^^unsure how to translate this part
    Vd=np.dot(normal,F)

    if Vd >= 0.0:
        dxbuilding = HUGE
    V0=-[np.dot(normal,Vecip1)+d]
    t=V0/Vd
    if(t < 0.0):
        dxbuilding=HUGE
        behind=1

    intersection[1]=Vecip1[1]+F[1]*t
    intersection[2]=Vecip1[2]+F[2]*t
    intersection[3]=Vecip1[3]+F[3]*t
    maximum=max(abs(normal[1]),abs(normal[2]),abs(normal[3]))
    if(maximum == abs(normal[1])):
        for P in range(1,size):
            G[P,1:2]= (intersection[2]-BuildingPoints[int(PolyArray[Q,1+P]),2],intersection[3]-BuildingPoints[int(PolyArray(Q,1+P)),3])
            #syntax check
    #*****************************************************************************



def PLANE(Vecip1, B1, B2, planehit,nbox):
#George:It sure would be a mess if there was a typo anywhere in here
    #This function calculates the normal at the hitpoint of a box.
    if planehit == 1:
        if Vecip1[1] == B1[1]:
            Point2=(B1[1],B1[2],B2[3])  
            Point3=(B1[1],B2[2],B1[3]) 
            CROSS((Point2-B1),(Point3-B1),nbox)
        elif Vecip1[1] == B2(1) :
            Point1=(B2[1],B1[2],B1[3])
            Point2=(B2[1],B1[2],B2[3])
            Point3=(B2[1],B2[2],B1[3])
            CROSS((Point3-Point1),(Point2-Point1),nbox)
    if planehit == 2:
        if Vecip1[2] == B1[2]:
            Point2=(B2[1], B1[2], B1[3])  
            Point3=(B1[1], B1[2], B2[3]) 
            CROSS((Point2-B1),(Point3-B1),nbox)
        elif Vecip1[2] == B2[2]: 
            Point1=(B1[1],B2[2],B1[1])
            Point2=(B1[1],B2[2],B2[3])
            Point3=(B2[1],B2[2],B1[3])
            CROSS((Point2-Point1),(Point3-Point1),nbox)
    if planehit == 3 :
         if Vecip1[3] == B1[3]:
            Point2=(B2[1], B1[2], B1[3])  
            Point3=(B1[1], B2[2], B1[3]) 
            CROSS((Point3-B1),(Point2-B1),nbox)
         elif Vecip1[3] == B2[3]:
            Point2=(B1[1],B2[2],B2[3])
            Point3=(B2[1],B1[2],B2[3])
            CROSS((Point2-B2),(Point3-B2),nbox)
    return PLANE