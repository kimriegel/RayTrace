#     BigBertha
#   Same as NASABOOM1EMParameterFile
#   Contains edits at the end 
#   These were made for compatibility and to eliminate bulk in main files 
import numpy as np

INPUTFILE = "inputNASABOOM1.txt"
Fs=24000.0
xinitial=145.0
yinitial=35.0
zinitial=0.0
radius=.15
soundspeed=348.537
ps=1.0
Temp=302.182778
time=.01
hr=20.0
theta=1.6863372
phi=3.44458181
boomspacing= 1 #2**-7 #.035#1   #.6
xmin=-1
ymin=30.0
zmin=0.0
xmax=-1
ymax=100.0
zmax=25.0
IMAX=75
h=10.0
absorbplanes=1
#allocate(tempalphabuilding(absorbplanes,8))
#Find way to rephrase
OUTPUTFILE='PythonTest.txt'
#Turn Radiosity on or off.  This will include diffuse reflections
radiosity=0
#Turn on complex absorption
complexabsorption=0

tempalphabuilding=np.zeros([absorbplanes,8])
if(complexabsorption == 1):
    
    tempalphabuilding[1]=[0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04]
else:
    tempalphabuilding=np.zeros([1,8])
#Enter an array for absorption of alpha ground octave bands between
#63 and 8000
tempalphaground=[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]
#Enter an array for absorption of Alpha Building octave bands between
#63 and 8000
tempalphabuilding[0]=[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]
#what percentage of the energy is reflected diffusely between 0,1
percentdiffuse=0.0

# #Broken all down to:
# complexabsorption = 1
# if complexabsorption == 1:
#     tempalphaground=np.array([[0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04],
#     [0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03],[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]])
# print(tempalphaground)


# Edits beyond here
#import Functions as fun
import math as m 
Vinitial=np.array([xinitial,yinitial,zinitial])
#xiinitial=m.cos(PF.phi)*m.sin(PF.theta)
#ninitial=m.sin(PF.phi)*m.sin(PF.theta)
#zetainitial=m.cos(PF.theta)
#length=m.sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*zetainitial)
#Finitial=np.array([xiinitial,ninitial,zetainitial])
                  #Show everything 
Finitial =np.array ([(m.cos(phi)*m.sin(theta)),(m.sin(phi)*m.sin(theta)),(m.cos(theta))])
tmp=(Finitial[0]*Vinitial[0]+Finitial[1]*Vinitial[1]+Finitial[2]*Vinitial[2])
#PLANEABC=np.array([Finitial[0],Finitial[1],Finitial[2],tmp])
PLANEABC = np.append(Finitial,tmp)

    #       Create initial boom array
yspace=boomspacing*abs(m.cos(phi))
zspace=boomspacing*abs(m.sin(theta))
if (xmin == xmax):
   RAYMAX=int((ymax-ymin)/yspace)*int((zmax-zmin)/zspace)   #This is the one we use. The others wil crash us
elif(ymin == ymax):
   RAYMAX=int((xmax-xmin)/xspace)*int((zmax-zmin)/zspace)
elif(zmin == zmax):
   RAYMAX=int((ymax-ymin)/yspace)*int((xmax-xmin)/xspace)
boomarray = np.zeros((RAYMAX,2))
print(RAYMAX , ' is the RAYMAX')
    #boomarray,sizex,sizey,sizez=fun.InitialGrid(boomspacing,PLANEABC[0],PLANEABC[1],PLANEABC[2],PLANEABC[3],theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,RAYMAX)
    # How about :
#boomarray,sizex,sizey,sizez=fun.InitialGrid(boomspacing,PLANEABC,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,RAYMAX)



def InitialGrid(radius,plane,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,arraysize):
    '''
    This function creates an equally spaced grid of size step apart

    Radius, the plane, its dimensions, and distances between parts 
     are defined in the parameterfile.    
    '''

    receiverarray=np.zeros((arraysize,3))
        #Both of these are defined the same as outsie the function.
    yspace=radius*abs(m.cos(phi))   
    zspace=radius*abs(m.sin(theta))
    xspace = 0  # pls define me 
    count = 0
    #i=0
    #j=0
        # For plane math
    A = plane[0]
    B = plane[1]
    C = plane[2]
    D = plane[3]    #D is tmp. I don't know what tmp is, and it would help a lot with understanding this math
        #These names are easier to search
    #Aleph = plane[0]
    #Bet   = plane[1]
    #Vet   = plane[2]
    #Gimel = plane[3]
    
    #boomarray = ((arraysize,3))    # For reference
    #zs = np.arange(1,int((zmax-zmin)//zspace)+1 )
    #ys = np.arange(1,int((ymax-ymin)//yspace)+1 )
    #    # For X
    #j = ys
    #i = zs
    #barray1 = ((Gimel-Bet) * (ymin+ys[:,None]*yspace) - Vet * (zmin+np.rot90(zs[:,None])*zspace))/Aleph
    #barray2 = ymin + j * yspace
    #barray3 = zmin + i * zspace
    #    #Has to be in reverse order to match 
    #        #I'll come back to this later
    #boomarray =np.array([[barray3],[barray2],[barray1]])
    ##boomarray = np.rotate(boomarray,3)
    #sizex=int((ymax-ymin)//(yspace))
    #sizey=int((zmax-zmin)//zspace)
    #sizez=1

    if xmin == xmax:    # This is the only working one, and therefore the one we use for our tests.
        for i in range(0,int((zmax-zmin)//zspace)):
            for j in range(0,int((ymax-ymin)//yspace)):
                receiverarray[count,0]=(D-B*(ymin+(j+1)*yspace)-C*(zmin+(i+1)*zspace))/A
                receiverarray[count,1]=ymin+(j+1)*yspace
                receiverarray[count,2]=zmin+(i+1)*zspace
                count=count+1
        sizex=int((ymax-ymin)//(yspace))
        sizey=int((zmax-zmin)//zspace)
        sizez=1
        #print(receiverarray[:,2])
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

boomarray,sizex,sizey,sizez=InitialGrid(boomspacing,PLANEABC,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,RAYMAX)
