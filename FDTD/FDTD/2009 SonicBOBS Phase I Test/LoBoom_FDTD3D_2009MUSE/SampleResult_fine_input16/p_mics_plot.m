%% p_mics_plot.m
%
%   written by Sang Ik Cho
%   last modified on 7 / 11 / 11
%   
%       This program plots the pressure time series of the low-boom
%   at microphone locations around the Museum building from 2009 SonicBobs
%   test recorded during an FDTD simulation.  The numerical data is read 
%   from an output file of the FDTD program called "p_mics.dat".  The plot 
%   is saved to a .tif image file in the same working directory.

clear all, close all;
scrsz = get(0,'ScreenSize');

p_mics = load('p_mics_fixed.dat');
num_mics = size(p_mics,2)-1;
tt = p_mics(:,1);
mics = p_mics(:,2:num_mics+1);
p_max = max(max(mics));
p_min = -max(max(-mics));

LineColors = [[0.0 0.7 0.0];...
              [1.0 0.0 0.0];...
              [1.0 0.0 1.0];...
              [0.0 0.0 1.0];...
              [0.7 0.9 0.2];...
              [0.9 0.6 0.0];...
              [0.7 0.0 0.7];...
              [0.8 0.8 0.2];...
              [0.0 0.7 0.9]];

figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2 scrsz(4)/1.8]);
h = axes('FontSize',16);
for l = 1:num_mics
    plot(tt,mics(:,l),'Linewidth',1.5,'Color',LineColors(mod(l-1,9)+1,:),'LineStyle','-'), hold on;
end
axis([0 tt(end) 1.2*p_min 1.2*p_max]);
% axis([0 tt(end) -0.8 1.5]);
xlabel('Time [s]','fontsize', 18); 
ylabel('Pressure [psf]','fontsize', 18);
title('Time trace of pressure at mics around the building');
set(gcf, 'Color', 'w');
print -dtiff -r300 p_mics.tif

