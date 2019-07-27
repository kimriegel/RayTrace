%TimeWaveForm2.m creates the Matrix A (5x4x18000 Matrix) and writes the time/pressure data for each
%receiver to a .dat file with the correct format for use in
%DataGraphingTest.m

clear all
format short
% Read in the data from the Ray Trace Code
fid='C:\Users\Kory George\Documents\RayTrace\MatlabScripts\PythonTest.txt';
A=ReadTecPlot(fid,'Single Point');

%A(3,4,:)
%fprintf(fid,'%d\n', A(3,4,:))
%fprintf('%d\t\n', A(3,4,:)) %Personal prefference
%for tick = 0:3000:18000
%    fprintf(tick, A(3:4:tick))
%end
%tick=1;
tconst = 4.16666662E-05;

%sample
%fprintf('Tertiary Receiver\n')
%fprintf('Pressure versus Time of Sonic Boom Data \n')
%fprintf('Time (Seconds) \t Pressure(Units) \n')
%for tick = 1:1:18 %Only for tests and showing others
%for tick = 1:1:18000 %the actual useful code
%    fprintf('%d',(tick*tconst))
%    fprintf('\t')
%    fprintf('%d\t\n', A(3,4,tick)) 
%It works?
%end


%All In Graph
%fprintf('All Receivers\n')
%fprintf('Pressure Versus Time of Each Receiver\n')
%fprintf('Time (Seconds) \t Receiver 1 \t Receiver 2 \t Receiver 3 \t Receiver 4 \t Receiver 5\n')
%for tick = 1:1000:18000
%       fprintf('%d', (tick*tconst))
%       fprintf('\t \t')
%       fprintf('%d\t', A(1,4,tick))
%       fprintf('\t \t\t')
%       fprintf('%d\t', A(2,4,tick))
%       fprintf('\t')
%       fprintf('%d\t', A(3,4,tick))
%       fprintf('\t')
%       fprintf('%d\t', A(4,4,tick))
%       fprintf('\t')
%       fprintf('%d\t\n', A(5,4,tick))
%end

%All In Graph
%fprintf('All Receivers\n')
%fprintf('Pressure Versus Time of Each Receiver\n')
%fprintf('Time (Seconds) \t Receiver 1 \t Receiver 2 \t Receiver 3 \t \t \t Receiver 4 \t \t Receiver 5\n')
%for tick = 1:1000:18000
%       fprintf('%d\t', (tick*tconst))
%       for rec = 1:1:5
%           fprintf('\t \t')
%           fprintf('%d\t', A(rec,4,tick))
%       end
%       fprintf('\n')
%end


fprintf('TITLE     = "Time Series Plot" \n')
fprintf('VARIABLES = "Solution Time"\n')
fprintf('"P[psf]"\n')
fprintf('ZONE T="Time Series Plot Zone"\n')
fprintf(' STRANDID=0, SOLUTIONTIME=0\n')
fprintf(' I=18000, J=1, K=1, ZONETYPE=Ordered\n')
fprintf(' DATAPACKING=POINT\n')
fprintf(' DT=(SINGLE SINGLE )\n')
for tick = 1:6000:18000 
    fprintf('%d',(tick*tconst))
    fprintf('\t')
    fprintf('%d\t\n', A(3,4,tick)) 
%It works?
end