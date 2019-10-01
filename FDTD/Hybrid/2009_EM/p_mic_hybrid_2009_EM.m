%% p_mic_hybrid_2009_EM.m
%
%   written by:     Sang Ik Cho
%   last modifed:   9 / 21 / 2011
%   e-mail:         stc142@psu.edu
%
%       This script uses GenerateHybrid.m & VerifyHybrid.m to generate the 
%   hybrid result of the pressure time histories recorded during the FDTD 
%   and ray tracing simulations of low-boom propagation and diffraction 
%   around the EM building during 2009 SonicBOBS experiment and compare it 
%   with experimental data for verification.  The hybrid result is saved to
%   a file called p_mics_hybrid.dat.  The first column of the file stores
%   the time vector and the subsequent columns store the pressure time 
%   histories at the microphones.

clc, clear all, close all;


%% USER-SPECIFIED INPUTS
% Load FDTD pressure time history result and set up a matrix with
% microphone locations where hybrid results are to be generated.  Also
% specify the sampling frequency for the time history data.
% *** NOTE: The first column of the matrix should correspond to the basis 
% microphone location where the shock arrival times of the FDTD and ray
% tracing results are aligned to make the hybrid result.  The time line for
% the other microphones then follow this basis microphone. ***
p_mics_FDTD = load('.\FDTD_Result\p_mics.dat');
p_mic_FDTD_19 = p_mics_FDTD(:,2);
p_mic_FDTD_29 = p_mics_FDTD(:,3);
p_FDTD = [p_mic_FDTD_19, p_mic_FDTD_29];  % Mic 19 selected as the basis microphone
tt_FDTD = p_mics_FDTD(:,1);
dt_FDTD = tt_FDTD(2)-tt_FDTD(1);
fs_FDTD = 1/dt_FDTD;

% Load ray tracing pressure time history result and set up a matrix with
% microphone locations where hybrid results are to be generated.  Also
% specify the sampling frequency for the time history data.
% *** NOTE: The Tecplot headers in the ray tracing results should be
% removed before loading them into Matlab.  The matrix must be set up with 
% the same order of microphones as the FDTD result. ***
p_mic_Rays_19 = load('.\Rays_Result\p_mic_19.dat');
p_mic_Rays_29 = load('.\Rays_Result\p_mic_29.dat');
p_Rays = [p_mic_Rays_19(:,2), p_mic_Rays_29(:,2)];
tt_Rays = p_mic_Rays_19(:,1);
dt_Rays = tt_Rays(2)-tt_Rays(1);
fs_Rays = 1/dt_Rays;

% Load experimentally recorded pressure time history and set up a matrix 
% with the same microphone locations specified for the hybrid result.  Also
% specify the sampling frequency of the recording.
% *** NOTE: The matrix must be set up with the same order of microphones as
% the FDTD result. ***
p_mic_NASA_19 = load('.\Experimental_Result\EMBuildingReciever2.txt');
p_mic_NASA_29 = load('.\Experimental_Result\EMBuildingReciever3.txt');
p_NASA = [p_mic_NASA_19, p_mic_NASA_29];
fs_NASA = 24000;

% Specify the filter parameters to generate the Parks-McClellan FIR filter
% used in hybrid modeling.
% *** NOTE: f_crossover should be set such that only appropriate frequency
% contents remain in the low-pass filtered FDTD result and high-pass 
% filtered ray tracing result.  Also, a large value of N_filt can cause 
% instability in the design algorithm and cause it to fail. ***
f_crossover = 70;  % The crossover frequency of the low-pass and high-pass filter pair.
f_margin = 10;      % The frequency roll off margin on both sides of the crossover frequency
N_filt = 1000;      % The FIR filter will be of order (2N_filt-1) with a delay of N_filt samples.

% Select the method to use in generating the hybrid result.
%   Method 1: Find the midpoint of initial shock by setting a pressure threshold at the half of the peak
%   Method 2: Match the time indices at the peak of derivatives of each numerical recordings
%   Method 3: Find the midpoint of initial shock by averaging the times where pressure reaches 10% and 90% of the peak
method_hybrid = 3;

% Select the method to use in comparing the hybrid result against the measurement.
%   Method 1: Find the midpoint of initial shock by setting a pressure threshold at the half of the peak
%   Method 2: Match the time indices at the peak of derivatives of each numerical recordings
%   Method 3: Find the midpoint of initial shock by averaging the times where pressure reaches 10% and 90% of the peak
method_compare = 3;


%% Function Calling
display 'Running GenerateHybrid...'
[p_hybrid,fs_hybrid] = GenerateHybrid(method_hybrid,p_FDTD,fs_FDTD,p_Rays,fs_Rays,f_crossover,f_margin,N_filt);

% Export the hybrid result to a file
N_hybrid = length(p_hybrid);
dt_hybrid = 1/fs_hybrid;
tt_hybrid = dt_hybrid*(0:N_hybrid-1);
f_p_hybrid = fopen('p_mics_hybrid.dat','w');
for i=1:N_hybrid
    fprintf(f_p_hybrid,'%10f %10f %10f\n',[tt_hybrid(i) p_hybrid(i,:)]);
end
fclose(f_p_hybrid);

% Visualize the hybrid result and verify it against the experimental data
display 'Running VerifyHybrid...'
VerifyHybrid(method_compare,p_hybrid,fs_hybrid,p_NASA,fs_NASA);

display '*******************************************************'
display '     Hybrid results are successfully generated'
display '*******************************************************'