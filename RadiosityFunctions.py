# $Revision: 1.0 $ : $Date:$ : $Author: kim

#***------------------------------------------------------------***
#***-------------------------- PATCHES -------------------------***
#***------------------------------------------------------------***

#     This function creates patch lengths using a geometric sum.

#use For loop instead of Do loop
#REPLACE l WITH L
def PATCHESSHORT(min,max,N,q,dd):
    W=max-min
    k = []
#     For box 1 we mesh the x direction
    for m in range (1,N):
        k[m] = W/2.0*(1-q)/(1-q**(N/2))
        if m <= (N/2):
            dd[m] = k[m]*q**(m-1)
        elif (m <= N) and (m >= (N/2)):
            dd[m] = k[m]*q**(N-m)
    return k
    #No errors returned. Does it work?
    #Find good values to plug


#Unfinished
def CREATEPATCHARRAY(ddm,ddl,Nm,Nl,x1,x2,y1,y2,z1,z2,patcharray,slope,b,slope1,b1,count,normal):
#This function creates an array of patches for each plane
#Raws
    double precision x1,y1,x2,y2,z1,z2,normal(3)
    #what is normal? a function or array? If array, why specificially index 3
    real temp1(3),temp2(3),temp3(3),vec1(3),vec2(3)
    #From Fortran, this is an array. Why is it referring to index 3.
    real vec3(3),slope,b,ddz,ddy,ddx,zcenter,xcenter,ycenter
    real slope1,b1
    double precision patcharray(Nm*ddL,6),x,y,z,ddm(Nm),ddL(ddL),d
    d = -x1*normal(1)-y1*normal(2)-z1*normal(3)
    #normal looks like a string with 3 inputs

#find how to do that
    if z1 == z2:
        count = 1
        print('ZPLANE',ddl, Nm, x1, y1, z1)
        for L in range (1,ddl):
            for m in range(1,Nm):
                x=x1-0.5*ddm[m]+sum(ddm[1:m])
                y=y1-0.5*ddl[L]+sum(ddl[1:L])
                #Translate SUM
                z = (-d-normal[0]*x-normal[1]*y)/normal[2]
                if m == 1:
                    zcenter = z
                if L == 1:
                    ddz = z-z1
                    #If loop may kill function early
                    #if problems come try While loop?
                else:
                    ddz = 0.5*(z-zcenter)
                if (y > slope*x+b) or (y < slope1*x+b1):
#                    GOTO 3
#Meep
#                patcharray(count,1)=x
#                patcharray(count,2)=y
#                patcharray(count,3)=z
#                patcharray(count,4)=ddm(m)
#                patcharray(count,5)=ddL(L)
#                patcharray(count,6)=ddz
                patcharray = [x,y,z,ddm[m],ddl[L],ddz]
                count=count+1

    elif y1 == y2:
        count=1
        print('YPLANE',Nl, Nm, x1, y1, z1)
        for L in range(1,Nl):
            for m in range(1,Nm):
                x = x1-0.5*ddm[m] + sum(ddm[1:m])
                z = z1-.5*ddl[L]+sum(ddl[1:L])
               #sum here
                y = (-d-normal[0]*x-normal[2]*z)/normal[1]
                if(m == 1):
                    ycenter = y
                if(L == 1):
                    ddy = y-y1
                else:
                    ddy = 0.5*(y-ycenter)
                if(z > slope*x+b) or (z < slope1*x+b1):
#                  GOTO 5
#               patcharray(count,1)=x
#               patcharray(count,2)=y
#               patcharray(count,3)=z
#               patcharray(count,4)=ddm(m)
#               patcharray(count,5)=ddy
#                patcharray(count,6)=ddl(l)
                patcharray = [x,y,z,ddm[m],ddy,ddL[L]]
                count=count+1
                #count is used to move thru patcharray string
                #I think
    elif x1 == x2:
        count = 1
        print('XPLANE', Nl, Nm, x1, y1, z1)
        for L in range (1, Nl):
            for m in range(1, Nm):
                y = y1- .5*ddm[m] + sum(ddm[1:m])
                z = z1- .5*ddl[l] + sum(dd1[1:l])
                x = (-d-normal[2]*z-normal[1]*y)/normal[0]

                if m == 1:
                    xcenter=x
                if(L == 1):
                    ddx= x-x1
                else:
                    ddx= 0.5*(x-xcenter)
                if(z > slope*y+b) or (z < slope1*y+b1):
                    patcharray = [x,y,z,ddx,ddm[m],ddl[l]]
                    count= count+1

# 5          CONTINUE
# 4       CONTINUE
#      elseif(x1.eq.x2)then
#         count=1
#         print*, 'XPLANE',Nl, Nm, x1, y1, z1
#         DO 6 l=1, Nl
#            DO 7 m=1, Nm
#
#               y=y1-.5*ddm(m)+SUM(ddm(1:m))
#               z=z1-.5*ddl(l)+SUM(ddl(1:l))
#               x=(-d-normal(3)*z-normal(2)*y)/normal(1)
#               if(m.eq.1)then
#                  xcenter=x
#               endif
#               if(l.eq.1)then
#                  ddx=x-x1
#               else
#                  ddx=0.5*(x-xcenter)
#               endif
#               if(z.gt.slope*y+b.or.z.lt.slope1*y+b1)then
#                  GOTO 7
#               endif
##               patcharray(count,1)=x
##               patcharray(count,2)=y
##               patcharray(count,3)=z
##               patcharray(count,4)=ddx
##               patcharray(count,5)=ddm(m)
##               patcharray(count,6)=ddl(l)
#                patcharray = [x,y,z,ddx,ddm(m),ddL(L)]
#               count=count+1
# 7          CONTINUE
# 6       CONTINUE
#      endif
#      END
#
