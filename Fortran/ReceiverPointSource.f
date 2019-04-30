C     Create a point reciever
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
      allocate(receiverarray(arraysize,3))
      allocate(receiverarray1(arraysize1,3))
      receiverarray(1,1)=93.4428213
      receiverarray(1,2)=28.8397178
      receiverarray(1,3)=.151
      receiverarray(2,1)=64.5832
      receiverarray(2,2)=28.5998
      receiverarray(2,3)=.151
      receiverarray(3,1)=64.5832
      receiverarray(3,2)=28.5998
      receiverarray(3,3)=7.9423
      receiverarray(4,1)=-2.40793
      receiverarray(4,2)=31.5003401
      receiverarray(4,3)=.151
      receiverarray(5,1)=75.11005
      receiverarray(5,2)=28.4945787
      receiverarray(5,3)=.151
      receiverarray1=recieverarray
      print*,'created receiver array'
