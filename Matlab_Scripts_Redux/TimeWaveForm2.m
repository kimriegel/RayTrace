%TimeWaveForm2.m creates the Matrix A (5x4x18000 Matrix) and writes the time/pressure data for each
%receiver to a .dat file with the correct format for use in
%DataGraphingTest.m

clear all
format short
% Read in the data from the Ray Trace Code
%fid='C:\Users\Student\Documents\MATLAB\MatLabScripts\FortranTest3.dat';
fid='FortranTest3.dat';

textfile = fopen('FortranTest3Trial2.txt', 'a+');
A=ReadTecPlotWorks(fid,'Single Point')

fprintf (textfile, '%s \n', 'TITLE     = "Time Series Plot"')
fprintf (textfile, '%s \n', 'VARIABLES = "Solution Time"')
fprintf (textfile, '%s \n', '"P[psf]"')
fprintf (textfile, '%s \n', 'ZONE T="Time Series Plot Zone"')
fprintf (textfile, '%s \n', 'STRANDID=0, SOLUTIONTIME=0')
fprintf (textfile, '%s \n', 'I=18000, J=1, K=1, ZONETYPE=Ordered')
fprintf (textfile, '%s \n', 'DATAPACKING=POINT')
fprintf (textfile, '%s \n', 'DT=(SINGLE SINGLE )')
fprintf (textfile, '%d, \n', A(3,4,:))


