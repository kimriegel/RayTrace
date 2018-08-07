      FaceNormalNo=5
      allocate(FaceNormals(FaceNormalNo,3))
      FaceNormals(1,1:3)=(/-1,0,0/)
      FaceNormals(2,1:3)=(/0,1,0/)
      FaceNormals(3,1:3)=(/1,0,0/)
      FaceNormals(4,1:3)=(/0,-1,0/)
      FaceNormals(5,1:3)=(/0,0,1/)

      Boxnumber=1
      allocate(Boxarraynear(boxnumber,3))
      allocate(Boxarrayfar(boxnumber,3))
      Boxarraynear(1,1:3)=(/10,10,0/)
      boxarrayfar(1,1:3)=(/64.4322,46.9316,8.2423/)
      TriangleNumber=0
      SquareNumber=0
      PolyBuilding=0
