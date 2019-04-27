#     Create a point reciever
import numpy as np
planenum=1
planename1='Single Point'
arraysize=5
arraysize1=arraysize
sizex=2
sizey=2
sizez=1
sizex1=sizex
sizey1=sizey
sizez1=sizez
receiverarray=np.zeros((arraysize,3))
receiverarray1=np.zeros((arraysize1,3))
receiverarray[0,0]=93.4428213
receiverarray[0,1]=28.8397178
receiverarray[0,2]=.151
receiverarray[1,0]=64.5832
receiverarray[1,1]=28.5998
receiverarray[1,2]=.151
receiverarray[2,0]=64.5832
receiverarray[2,1]=28.5998
receiverarray[2,2]=7.9423
receiverarray[3,0]=-2.40793
receiverarray[3,1]=31.5003401
receiverarray[3,2]=.151
receiverarray[4,0]=75.11005
receiverarray[4,1]=28.4945787
receiverarray[4,2]=.151
receiverarray1=receiverarray
print('created receiver array')
