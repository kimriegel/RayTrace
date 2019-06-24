C     BigBertha
      INPUTFILE="inputNASABOOM1.txt"
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
      boomspacing=1
      xmin=-1
      ymin=30.0
      zmin=0.0
      xmax=-1
      ymax=100.0
      zmax=25.0
      IMAX=75
      h=10.0
      absorbplanes=1
      allocate(tempalphabuilding(absorbplanes,8))
      OUTPUTFILE='FortranTest1.dat'
C     Turn Radiosity on or off.  This will include diffuse reflections
      radiosity=0
C     Turn on complex absorption
      complexabsorption=0
      if(complexabsorption.eq.1)then
         tempalphabuilding(2,1:8)=(/0.55,0.55,0.25,0.18,0.12,0.07,0.04
     *        ,0.04/)
      endif
C     Enter an array for absorption of alpha ground octave bands between
C     63 and 8000
      tempalphaground=(/0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03/)
C     Enter an array for absorption of Alpha Building octave bands between
C     63 and 8000
      tempalphabuilding(1,1:8)=(/0.01,0.01,0.01,0.02,0.02,0.02,0.03,
     *     0.03/)
C     what percentage of the energy is reflected diffusely between 0,1
      percentdiffuse=0.1


