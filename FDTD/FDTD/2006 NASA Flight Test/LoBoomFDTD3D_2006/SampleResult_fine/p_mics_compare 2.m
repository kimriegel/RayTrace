%% p_mics_compare.m
%
%   written by Sang Ik Cho
%   last modified on 5 / 22 / 12
%   

clear all, close all;
scrsz = get(0,'ScreenSize');

%% Load Data
data_directory = 'C:\Documents and Settings\Sang Cho\My Documents\LoBoom\Field Test Data\2006 NASA Flight Test Data\';
day_num = 21;    % To select Precal on 6/21
session_num = 7;
session_time = 'Run 6_21_2006 9_26_11 AM_Boom';
boom_num = 1;
session_name = [session_time num2str(boom_num)];

data_name = [data_directory 'Mics Only Shortname\' session_name '.mat'];
load (data_name);
p_mics = load('p_mics.dat');

% Load the calibration data
Calib = xlsread([data_directory 'Calibration_more.xls']);
Shock_arrival = xlsread([data_directory 'Shock_arrival.xls']);

%% Prepare Waveforms
% Numerical
num_mics = size(p_mics,2)-1;
num_tt = p_mics(:,1);
num_p_mics = p_mics(:,2:num_mics+1);
num_179 = p_mics(:,2);
num_180 = p_mics(:,3);
num_181 = p_mics(:,4);
num_182 = p_mics(:,5);
num_183 = p_mics(:,6);
num_184 = p_mics(:,7);
num_185 = p_mics(:,8);
num_186 = p_mics(:,9);
num_187 = p_mics(:,10);
num_188 = p_mics(:,11);
num_189 = p_mics(:,12);
num_191 = p_mics(:,13);
num_192 = p_mics(:,14);
num_193 = p_mics(:,15);
num_194 = p_mics(:,16);
num_195 = p_mics(:,17);
num_196 = p_mics(:,18);
num_197 = p_mics(:,19);
num_198 = p_mics(:,20);
num_199 = p_mics(:,21);
num_200 = p_mics(:,22);
num_202 = p_mics(:,23);
num_207 = p_mics(:,24);
num_288 = p_mics(:,25);

p_max = max(max(num_p_mics));
p_min = -max(max(-num_p_mics));
num_dt = num_tt(2)-num_tt(1);
num_fs = 1/num_dt;
tmax = 0.4;

% Experimental
exp_dt = 1/SampFreq;
decimate_factor = 8;
fs = SampFreq/decimate_factor;
dt = 1/fs;
input_raw = Channel_190_Voltage/Calib(190,day_num)/47.88;
exp_179_raw = Channel_179_Voltage/Calib(179,day_num)/47.88;
exp_180_raw = Channel_180_Voltage/Calib(180,day_num)/47.88;
exp_181_raw = Channel_181_Voltage/Calib(181,day_num)/47.88;
exp_182_raw = Channel_182_Voltage/Calib(182,day_num)/47.88;
exp_183_raw = Channel_183_Voltage/Calib(183,day_num)/47.88;
exp_184_raw = Channel_184_Voltage/Calib(184,day_num)/47.88;
exp_185_raw = Channel_185_Voltage/Calib(185,day_num)/47.88;
exp_186_raw = Channel_186_Voltage/Calib(186,day_num)/47.88;
exp_187_raw = Channel_187_Voltage/Calib(187,day_num)/47.88;
exp_188_raw = Channel_188_Voltage/Calib(188,day_num)/47.88;
exp_189_raw = Channel_189_Voltage/Calib(189,day_num)/47.88;
exp_191_raw = Channel_191_Voltage/Calib(191,day_num)/47.88;
exp_192_raw = Channel_192_Voltage/Calib(192,day_num)/47.88;
exp_193_raw = Channel_193_Voltage/Calib(193,day_num)/47.88;
exp_194_raw = Channel_194_Voltage/Calib(194,day_num)/47.88;
exp_195_raw = Channel_195_Voltage/Calib(195,day_num)/47.88;
exp_196_raw = Channel_196_Voltage/Calib(196,day_num)/47.88;
exp_197_raw = Channel_197_Voltage/Calib(197,day_num)/47.88;
exp_198_raw = Channel_198_Voltage/Calib(198,day_num)/47.88;
exp_199_raw = Channel_199_Voltage/Calib(199,day_num)/47.88;
exp_200_raw = Channel_200_Voltage/Calib(200,day_num)/47.88;
exp_202_raw = Channel_202_Voltage/Calib(202,day_num)/47.88;
exp_207_raw = Channel_207_Voltage/Calib(207,day_num)/47.88;
exp_288_raw = Channel_288_Voltage/Calib(288,day_num)/47.88;

exp_arrival_time_191 = find(exp_191_raw>0.1);
exp_arrival_time_input = find(input_raw>0.1);
time_diff = (exp_arrival_time_191(1)-exp_arrival_time_input(1))*exp_dt;
num_arrival_time_191 = find(num_191>0.1);
pb = exp_arrival_time_191(1)-floor(num_arrival_time_191(1)*num_dt/exp_dt);
pe = pb+tmax/exp_dt;
exp_input = decimate(double(input_raw((pb:pe)-(time_diff+.01)/exp_dt)),decimate_factor);
exp_179 = decimate(double(exp_179_raw(pb:pe)),decimate_factor);
exp_180 = decimate(double(exp_180_raw(pb:pe)),decimate_factor);
exp_181 = decimate(double(exp_181_raw(pb:pe)),decimate_factor);
exp_182 = decimate(double(exp_182_raw(pb:pe)),decimate_factor);
exp_183 = decimate(double(exp_183_raw(pb:pe)),decimate_factor);
exp_184 = decimate(double(exp_184_raw(pb:pe)),decimate_factor);
exp_185 = decimate(double(exp_185_raw(pb:pe)),decimate_factor);
exp_186 = decimate(double(exp_186_raw(pb:pe)),decimate_factor);
exp_187 = decimate(double(exp_187_raw(pb:pe)),decimate_factor);
exp_188 = decimate(double(exp_188_raw(pb:pe)),decimate_factor);
exp_189 = decimate(double(exp_189_raw(pb:pe)),decimate_factor);
exp_191 = decimate(double(exp_191_raw(pb:pe)),decimate_factor);
exp_192 = decimate(double(exp_192_raw(pb:pe)),decimate_factor);
exp_193 = decimate(double(exp_193_raw(pb:pe)),decimate_factor);
exp_194 = decimate(double(exp_194_raw(pb:pe)),decimate_factor);
exp_195 = decimate(double(exp_195_raw(pb:pe)),decimate_factor);
exp_196 = decimate(double(exp_196_raw(pb:pe)),decimate_factor);
exp_197 = decimate(double(exp_197_raw(pb:pe)),decimate_factor);
exp_198 = decimate(double(exp_198_raw(pb:pe)),decimate_factor);
exp_199 = decimate(double(exp_199_raw(pb:pe)),decimate_factor);
exp_200 = decimate(double(exp_200_raw(pb:pe)),decimate_factor);
exp_202 = decimate(double(exp_202_raw(pb:pe)),decimate_factor);
exp_207 = decimate(double(exp_207_raw(pb:pe)),decimate_factor);
exp_288 = decimate(double(exp_288_raw(pb:pe)),decimate_factor);
N = length(exp_179);
tt = (0:N-1)*dt;
num_179 = interp1(num_tt,num_179,tt,'spline');
num_180 = interp1(num_tt,num_180,tt,'spline');
num_181 = interp1(num_tt,num_181,tt,'spline');
num_182 = interp1(num_tt,num_182,tt,'spline');
num_183 = interp1(num_tt,num_183,tt,'spline');
num_184 = interp1(num_tt,num_184,tt,'spline');
num_185 = interp1(num_tt,num_185,tt,'spline');
num_186 = interp1(num_tt,num_186,tt,'spline');
num_187 = interp1(num_tt,num_187,tt,'spline');
num_188 = interp1(num_tt,num_188,tt,'spline');
num_189 = interp1(num_tt,num_189,tt,'spline');
num_191 = interp1(num_tt,num_191,tt,'spline');
num_192 = interp1(num_tt,num_192,tt,'spline');
num_193 = interp1(num_tt,num_193,tt,'spline');
num_194 = interp1(num_tt,num_194,tt,'spline');
num_195 = interp1(num_tt,num_195,tt,'spline');
num_196 = interp1(num_tt,num_196,tt,'spline');
num_197 = interp1(num_tt,num_197,tt,'spline');
num_198 = interp1(num_tt,num_198,tt,'spline');
num_199 = interp1(num_tt,num_199,tt,'spline');
num_200 = interp1(num_tt,num_200,tt,'spline');
num_202 = interp1(num_tt,num_202,tt,'spline');
num_207 = interp1(num_tt,num_207,tt,'spline');
num_288 = interp1(num_tt,num_288,tt,'spline');

% Mic 179 (R1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_179,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_179,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic179','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_R1(179).tif -a1 -r600 -painters

% Mic 180 (R2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_180,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_180,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic180','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_R2(180).tif -a1 -r600 -painters

% Mic 181 (R3)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_181,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_181,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic181','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_R3(181).tif -a1 -r600 -painters

% Mic 182 (R4)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_182,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_182,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic182','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_R4(182).tif -a1 -r600 -painters

% Mic 183 (R5)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_183,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_183,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic183','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_R5(183).tif -a1 -r600 -painters

% Mic 184 (B2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_184,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_184,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic184','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_B2(184).tif -a1 -r600 -painters

% % Mic 185
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
% h = axes('FontSize',16);
% plot(tt,exp_185,'LineStyle','-','LineWidth',2,'color','b');
% hold on, plot(tt,num_185,'LineStyle','--','LineWidth',2,'color','r');
% xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
% legend('Experimental','Numerical');
% axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
% text(0.2,1.3,'Mic185','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
% set(gcf, 'Color', 'w');
% export_fig p_mic_185.tif -a1 -r600 -painters

% Mic 186 (G5)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_186,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_186,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic186','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G5(186).tif -a1 -r600 -painters

% Mic 187 (G6)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_187,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_187,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic187','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G6(187).tif -a1 -r600 -painters

% Mic 188 (Y1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_188,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_188,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic188','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_Y1(188).tif -a1 -r600 -painters

% Mic 189 (Y2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_189,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_189,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic189','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_Y2(189).tif -a1 -r600 -painters

% Mic 191 (G2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_191,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_191,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic191','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G2(191).tif -a1 -r600 -painters

% Mic 192 (G3)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_192,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_192,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic192','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G3(192).tif -a1 -r600 -painters

% % Mic 193
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
% h = axes('FontSize',16);
% plot(tt,exp_193,'LineStyle','-','LineWidth',2,'color','b');
% hold on, plot(tt,num_193,'LineStyle','--','LineWidth',2,'color','r');
% xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
% legend('Experimental','Numerical');
% axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
% text(0.2,1.3,'Mic193','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
% set(gcf, 'Color', 'w');
% export_fig p_mic_193.tif -a1 -r600 -painters

% Mic 194 (G1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_194,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_194,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic194','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G1(194).tif -a1 -r600 -painters

% % Mic 195
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
% h = axes('FontSize',16);
% plot(tt,exp_195,'LineStyle','-','LineWidth',2,'color','b');
% hold on, plot(tt,num_195,'LineStyle','--','LineWidth',2,'color','r');
% xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
% legend('Experimental','Numerical');
% axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
% text(0.2,1.3,'Mic195','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
% set(gcf, 'Color', 'w');
% export_fig p_mic_195.tif -a1 -r600 -painters

% Mic 196 (F1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_196,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_196,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic196','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_F1(196).tif -a1 -r600 -painters

% Mic 197 (F2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_197,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_197,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic197','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_F2(197).tif -a1 -r600 -painters

% Mic 198 (B1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_198,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_198,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic198','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_B1(198).tif -a1 -r600 -painters

% Mic 199 (G4)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_199,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_199,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic199','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_G4(199).tif -a1 -r600 -painters

% % Mic 200
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
% h = axes('FontSize',16);
% plot(tt,exp_200,'LineStyle','-','LineWidth',2,'color','b');
% hold on, plot(tt,num_200,'LineStyle','--','LineWidth',2,'color','r');
% xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
% legend('Experimental','Numerical');
% axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
% text(0.2,1.3,'Mic200','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
% set(gcf, 'Color', 'w');
% export_fig p_mic_200.tif -a1 -r600 -painters

% % Mic 202
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
% h = axes('FontSize',16);
% plot(tt,exp_202,'LineStyle','-','LineWidth',2,'color','b');
% hold on, plot(tt,num_202,'LineStyle','--','LineWidth',2,'color','r');
% xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
% legend('Experimental','Numerical');
% axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
% text(0.2,1.3,'Mic202','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
% set(gcf, 'Color', 'w');
% export_fig p_mic_202.tif -a1 -r600 -painters

% Mic 207 (S1)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_207,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_207,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic207','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_S1(207).tif -a1 -r600 -painters

% Mic 288 (S2)
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2.5 scrsz(4)/3]);
h = axes('FontSize',16);
plot(tt,exp_288,'LineStyle','-','LineWidth',2,'color','b');
hold on, plot(tt,num_288,'LineStyle','--','LineWidth',2,'color','r');
xlabel('Time [s]','fontsize', 18), ylabel('Pressure [psf]','fontsize', 18);
legend('Experimental','Numerical');
axis([0 tt(end) 1.2*p_min 1.2*p_max]); 
text(0.2,1.3,'Mic288','FontName','Arial','FontSize',18,'FontWeight','demi','HorizontalAlignment','center');
set(gcf, 'Color', 'w');
export_fig p_mic_S2(288).tif -a1 -r600 -painters

close all;
