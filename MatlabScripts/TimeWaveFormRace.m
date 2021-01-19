%TimeWaveForm2.m creates the Matrix A (5x4x18000 Matrix) and writes the time/pressure data for each
%receiver to a .dat file with the correct format for use in
%DataGraphingTest.m

clear all
format short
% Read in the data from the Ray Trace Code
fid='C:\Users\Kory George\Documents\RayTrace\MatlabScripts\PythonTest.txt';
A=ReadTecPlotWorks(fid,'Single Point');
fPOutput_file = fopen('C:\Users\Kory George\Documents\RayTraceCode\MatlabScripts\ReceiverInput3.dat','w');
%^change here for receiver
tconst = 4.16666662E-05;


fprintf(fPOutput_file,'TITLE     = "Time Series Plot" \n');
fprintf(fPOutput_file,'VARIABLES = "Solution Time"\n');
fprintf(fPOutput_file,'"P[psf]"\n');
fprintf(fPOutput_file,'ZONE T="Time Series Plot Zone"\n');
fprintf(fPOutput_file,' STRANDID=0, SOLUTIONTIME=0\n');
fprintf(fPOutput_file,' I=18000, J=1, K=1, ZONETYPE=Ordered\n');
fprintf(fPOutput_file,' DATAPACKING=POINT\n');
fprintf(fPOutput_file,' DT=(SINGLE SINGLE )\n');
for tick = 1:1:18000
    fprintf(fPOutput_file,'%d',(tick*tconst));
    fprintf(fPOutput_file,' ');
    fprintf(fPOutput_file,'%d\t\n', A(3,4,tick)) ;
    %^Also here
%It works?
end
%tick is equal to I, change to make more concise. Or don't 
