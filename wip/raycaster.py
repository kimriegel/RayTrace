# Initialize all the important variables. 

import initializer as i

# import math libraries. 
import numpy as np 

#import the building Geometry. 
import BuildingGeometry as bg

# import relevant functions for Ray Trace
import Functions as fun
BG = bg
print('began rays')
# These Ray counts are where they hit the receivers when we need to check that things are working. 
#ray = 606                     # @ PF.boomspacing = 1
#ray = 455174                 # @ PF.boomspacing = 0.06
hitcount=0
doublehit=0
amplitude = i.frecuencias[:,1]/i.normalization
phase=i.frecuencias[:,2]
#if PF.defaultground:    # dummy code for now
#    ground = True 

def updateFreq(dx,alpha,diffusion):
    """
    Update ray phase and amplitude
    """
    nonlocal phase,amplitude        # works directly
    twopidx = twopi * dx
    tempphase = phase[:] - (twopidx/lamb)
    tempphase %= twopi
    phase = np.where( (tempphase > PI),      tempphase-twopi,      tempphase)
    amplitude *= ((1.0-alpha) * (1.0-diffusion) * np.exp(airabsorb*dx))

def cast():
    if (PF.h < (2*PF.radius)): 
        print('h is less than 2r')
        # break
        return
    F = np.array(Finitial)         # If not defined this way it will make them the same object. This will break the entire program. Do not change
    veci = boomarray[ray,:]
    for I in range(PF.IMAX):      # Making small steps along the ray path.  For each step we should return, location, phase and amplitude
        dxreceiver=HUGE
        # Find the closest sphere and store that as the distance
        for R in ears:
            # The way that tempreceiver works now, it's only used here and only should be used here. It's not defined inside the receiver because it's ray dependant.
            tempreceiver = R.SphereCheck(radius2,F,veci)    #distrance to receiver
            if (receiverhit >= 1):  #if you hit a receiver last time, don't hit it again
                if np.all(R.position ==lastreceiver):
                    tempreceiver=HUGE
                if np.all(F == checkdirection):
                    OC = R.position - veci
                    OCLength = np.dot(OC,OC)
                    if(OCLength < radius2):
                        tempreceiver=HUGE
            if(receiverhit >= 2):
                if np.all(R.position == lastreceiver):
                    tempreceiver=HUGE
                if (tempreceiver < dxreceiver):   
                    R.dxreceiver=tempreceiver
                    dxreceiver=tempreceiver
                    receiverpoint= R.position
                elif (tempreceiver == dxreceiver and tempreceiver != HUGE):
                    receivercheck=tempreceiver          
                    if np.all(R.position==receiverpoint):
                        doublehit=0
                    else:
                        #receiverpoint2 = R.position
                        R2 = R
                        doublehit=1
                        print('double hit')

            #     Check Intersection with ground plane
        if ground:  #if we are using our own ground or it is defined in the obj file
            GROUNDN=GROUNDABC
            GROUNDVD=GROUNDN[0]*F[0]+GROUNDN[1]*F[1]+GROUNDN[2]*F[2]
            if (groundhit==1):
                dxground=HUGE
            elif (GROUNDVD!=0.0):
                GROUNDVO=((GROUNDN[0]*veci[0]+GROUNDN[1]*veci[1]+GROUNDN[2]*veci[2])+GROUNDD)
                dxground1=(-1.000)*GROUNDVO*(1.000)/GROUNDVD
                dxground=dxground1
                Vecip1=veci+dxground*F
                tmp=(GROUNDABC[0]*Vecip1[0]+GROUNDABC[1]*Vecip1[1]+GROUNDABC[2]*Vecip1[2]+GROUNDD)                  
                if (dxground < 0.0):
                    dxground=HUGE
            else:
                dxground=HUGE

        #     Check intersection with building
        dxbuilding=HUGE
        hit=0
        planehit=0
        #     Check intersection with Boxes
        for Q in range(0, BG.BoxNumber):
            dxnear, dxfar, hit, planehit=fun.box(BG.BoxArrayNear[Q], BG.BoxArrayFar[Q], veci, F)
            if (dxnear < dxbuilding):
                dxbuilding=dxnear
                Vecip1=veci+np.multiply(dxbuilding,F)
                whichbox=Q
                nbox=fun.plane(Vecip1, BG.BoxArrayNear[whichbox], BG.BoxArrayFar[whichbox], planehit)
        #     Check intersection with Triangles
        if(BG.TriangleNumber > 0):
            for Q in range(0, BG.TriangleNumber):
                dxnear, behind = fun.Polygon(veci,F,Q,3,TriangleNumber,PointNumbers,Trianglearray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                if (dxnear < dxbuilding):
                    dxbuilding=dxnear
                    nbox=normal
                    whichbox=Q
        #     Check intersection with Squares
        if(BG.SquareNumber>0):
            for Q in range(0,BG.SquareNumber):
                    dxnear, behind=Polygon(veci,F,Q,4,SquareNumber,PointNumbers,SquareArray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                    if (dxnear < dxbuilding):
                        dxbuilding=dxnear
                        nbox=normal
                        whichbox=Q
        buildinghit=0
        receiverhit=0
        groundhit=0

        #     Check to see if ray hits within step size
        if (dxreceiver < PF.h or dxground < PF.h or dxbuilding < PF.h):
            dx=min(dxreceiver,dxground,dxbuilding)
            #tmpsum = tmpsum + dx
            #     if the ray hits a receiver, store in an array.  If the ray hits two, create two arrays to store in.
            for R in ears:
                if dx == R.dxreceiver:
                    print('Ray ',ray +1,' hit receiver ',R.recNumber,' at step ',I)
                    veci += (dx*F)
                    receiverhit=1
                    checkdirection=F
                    if(doublehit==1):
                          receiverhit=2
                    hitcount=hitcount+1
                    updateFreq(dx,alphanothing,0)
                    lastreceiver = receiverpoint
                    outputarray1[:,0] = frecuencias[:,0]
                    outputarray1[:,1:4] = receiverpoint[:]
                    outputarray1[:,5] = phase[:]
                    if doublehit == 1 :
                          #R2 = R      #Supposed to be other R, but just a placeholder for now
                          R.on_Hit(amplitude/2,phase)
                          R2.on_Hit(amplitude/2,phase)
                    else:
                          R.on_Hit(amplitude,phase)
            #     If the ray hits the ground then bounce off the ground and continue
            if (abs(dx-dxground)< 10.0**(-13.0)):                  
                #Vecip1=veci+np.multiply(dxground,F)
                veci += (dxground*F)
                tmp = np.dot(GROUNDABC,veci)

                if(tmp != GROUNDD): 
                    veci[2] = 0
                print('hit ground at ',I)
                dot1 = np.dot(F,nground)
                n2 = np.dot(nground,nground)
                F -= (2.0*(dot1/n2 *nground))
                length = np.sqrt(np.dot(F,F))
                groundhit=1
                twopidx=twopi*dxground
                #     Loop through all the frequencies
                updateFreq(dxground,alphaground,diffusionground)
                if(PF.radiosity==1 and (diffusionground!=0.0)):
                    for Q in range (0,PatchNo):
                        if (formfactors[0,Q,1]==1):
                            if(veci[0]<=(patcharray[Q,W,0]+0.5*patcharray[Q,W,3]) and veci[0]>=(patcharray[Q,W,0]-0.5*patcharray[Q,W,3])):
                                if(veci[1]<=(patcharray[Q,W,1]+0.5*patcharray[Q,W,4]) and veci[1]>=(patcharray[Q,W,1]-0.5*patcharray[Q,W,4])):
                                        if(veci[2]<=(patcharray[Q,W,2]+0.5*patcharray[Q,W,5]) and veci[2]>=(patcharray[Q,W,2]-0.5*patcharray[Q,W,5])):
                                            temp2=complex(abs(patcharray[Q,W,6])*np.exp(XJ*patcharray[Q,W,7]))
                                            temp3=complex(abs(amplitude[W]*(1.0-alphaground[W])*diffusionground*exp(-m*dxground))*exp(1j*phasefinal))
                                            temp4=temp2+temp3
                                            patcharray[Q,W,6]=abs(temp4)
                                            patcharray[Q,W,7]=np.arctan(temp4.imag,temp4.real)
            #     if the ray hits the building then change the direction and continue
            if (dx==dxbuilding):                  
                veci += (dx*F)
                print('hit building at step ',I)
                n2 = np.dot(nbox,nbox)
                nbuilding=nbox/np.sqrt(n2)
                dot1= np.dot(F,nbuilding)
                F -= (2.0*(dot1/n2 *nbuilding))
                length = np.sqrt(np.dot(F,F))
                buildinghit=1
                if PF.complexabsorption:
                    if PF.absorbplanes==2:
                        if (veci[2]>0.0) and (veci[2]<height1):
                            alpha = alphabuilding[0,:]
                        elif(veci[2]>height1 and veci[2]<=height2):
                                alpha=alphabuilding[1,:]
                    if(PF.absorbplanes==3):
                        if(veci[2]>height2 and veci[2] <=height3):
                                alpha=alphabuilding[2,:]
                    if(PF.absorbplanes==4):
                        if(veci[2]>height3):
                            alpha=alphabuilding[4,:]
                else:
                    alpha=alphabuilding[0,:]
                updateFreq(dx,alpha,diffusion)
                        
        #     If there was no interaction with buildings then proceed with one step. 
        else:
            veci += (PF.h*F)
            updateFreq(PF.h,alphanothing,0)
    print('finished ray', ray + 1)