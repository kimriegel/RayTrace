%% p_mic_analyze.m
%
%   written by Sang Ik Cho
%   last modifed on 2 / 16 / 2011
%
%   This program analyzes and compares the time series of pressure from the
%   FDTD to the NASA measurement data from 2009 SonicBOBS I test.

clear all, close all;
scrsz = get(0,'ScreenSize');

% %% NASA Measurement
% % Load C microphones data
% load('C:\Documents and Settings\Sang Cho\My Documents\LoBoom\Field Test Data\2009 SonicBOBS1 NASA & RT Data\Museum\NASA Data\Data\12September2009/SonicBOBS_091209180736.mat')
% BoomDate_C = FlightPass{1};                 % Date of the test
% BoomTime_C = FlightPass{2};                 % UTC time that Boom arrives at CA
% Aircraft_C = cell2mat(FlightPass{3});       % Aircraft_mics type
% FlightNum_C = FlightPass{4};                % Flight number
% BoomPass_C = FlightPass{5};                 % Which boom pass
% fs_NASA_C = 24000;
% dt_NASA_C = 1/fs_NASA_C;
% tt_NASA_raw_C = TimeSecUTC-BoomTime_C;      % Time w.r.t. BoomTime_mics (Recording always started 5sec before BoomTime_mics and ended 20sec after)
% 
% CalibratedMikeData_psf = CalibratedMikeData*47.88;
% p_NASA_raw_CA = CalibratedMikeData_psf(:,1);    % Pressure recording at CA
% p_NASA_raw_LB = CalibratedMikeData_psf(:,2);    % Pressure recording at LB
% p_NASA_raw_CC = CalibratedMikeData_psf(:,3);    % Pressure recording at CC
% p_NASA_raw_CD = CalibratedMikeData_psf(:,4);    % Pressure recording at CD
% p_NASA_raw_CE = CalibratedMikeData_psf(:,5);    % Pressure recording at CE
% p_NASA_raw_CF = CalibratedMikeData_psf(:,6);    % Pressure recording at CF
% p_NASA_raw_CG = CalibratedMikeData_psf(:,7);    % Pressure recording at CG
% p_NASA_raw_CH = CalibratedMikeData_psf(:,8);    % Pressure recording at CH
% p_NASA_raw_CI = CalibratedMikeData_psf(:,9);    % Pressure recording at CI
% p_NASA_raw_CJ = CalibratedMikeData_psf(:,10);   % Pressure recording at CJ
% 
% % Load W microphones data
% load('C:\Documents and Settings\Sang Cho\My Documents\LoBoom\Field Test Data\2009 SonicBOBS1 NASA & RT Data\Museum\NASA Data\Data\12September2009/SonicBOBSBB_091209180629.mat')
% BoomDate_W = FlightPass{1};                 % Date of the test
% BoomTime_W = FlightPass{2};                 % UTC time that Boom arrives at the 1st microphone in the list that successfully records the boom (should be WA if things went right)
% Aircraft_W = cell2mat(FlightPass{3});       % Aircraft type
% FlightNum_W = FlightPass{4};                % Flight number
% BoomPass_W = FlightPass{5};                 % Which boom pass
% fs_NASA_W = 1/(TimeSecUTC(2,1)-TimeSecUTC(1,1));  % Although the time codes for all W mics are different, their fs are the same.
% dt_NASA_W = 1/fs_NASA_W;
% 
% if (~strcmp(Aircraft_W,Aircraft_C) || FlightNum_W~=FlightNum_C || BoomPass_W~=BoomPass_C)
%     display('***** WARNING: Microphone data & BADS/BASS data being compared are for different boom events! *****');
% end
% tt_NASA_raw_WA = TimeSecUTC(:,1)-BoomTime_C;    % Time code at WA w.r.t boom arrival at CA
% tt_NASA_raw_WB = TimeSecUTC(:,2)-BoomTime_C;    % Time code at WB w.r.t boom arrival at CA
% tt_NASA_raw_WC = TimeSecUTC(:,3)-BoomTime_C;    % Time code at WC w.r.t boom arrival at CA
% tt_NASA_raw_WD = TimeSecUTC(:,4)-BoomTime_C;    % Time code at WD w.r.t boom arrival at CA
% tt_NASA_raw_WE = TimeSecUTC(:,5)-BoomTime_C;    % Time code at WE w.r.t boom arrival at CA
% tt_NASA_raw_WF = TimeSecUTC(:,6)-BoomTime_C;    % Time code at WF w.r.t boom arrival at CA
% tt_NASA_raw_WG = TimeSecUTC(:,7)-BoomTime_C;    % Time code at WG w.r.t boom arrival at CA
% tt_NASA_raw_WH = TimeSecUTC(:,8)-BoomTime_C;    % Time code at WH w.r.t boom arrival at CA
% tt_NASA_raw_WY = TimeSecUTC(:,9)-BoomTime_C;    % Time code at WY w.r.t boom arrival at CA
% tt_NASA_raw_WZ = TimeSecUTC(:,10)-BoomTime_C;   % Time code at WZ w.r.t boom arrival at CA
% 
% CalibratedOverpressure_psf = CalibratedOverpressure*47.88;
% p_NASA_raw_WA = CalibratedOverpressure_psf(:,1);    % Pressure recording at WA
% p_NASA_raw_WB = CalibratedOverpressure_psf(:,2);    % Pressure recording at WB
% p_NASA_raw_WC = CalibratedOverpressure_psf(:,3);    % Pressure recording at WC
% p_NASA_raw_WD = CalibratedOverpressure_psf(:,4);    % Pressure recording at WD
% p_NASA_raw_WE = CalibratedOverpressure_psf(:,5);    % Pressure recording at WE
% p_NASA_raw_WF = CalibratedOverpressure_psf(:,6);    % Pressure recording at WF
% p_NASA_raw_WG = CalibratedOverpressure_psf(:,7);    % Pressure recording at WG
% p_NASA_raw_WH = CalibratedOverpressure_psf(:,8);    % Pressure recording at WH
% p_NASA_raw_WY = CalibratedOverpressure_psf(:,9);    % Pressure recording at WY
% p_NASA_raw_WZ = CalibratedOverpressure_psf(:,10);   % Pressure recording at WZ
% 

%% FDTD Result
p_FDTD_raw = load ('p_mics.dat');
tt_FDTD_raw = p_FDTD_raw(:,1);
dt_FDTD_raw = p_FDTD_raw(2)-p_FDTD_raw(1);
fs_FDTD_raw = 1/dt_FDTD_raw;
N_FDTD_raw = length(tt_FDTD_raw);
T_FDTD_raw = dt_FDTD_raw*N_FDTD_raw;

p_FDTD_raw_179 = p_FDTD_raw(:,2);    % Pressure recording at CC
p_FDTD_raw_180 = p_FDTD_raw(:,3);    % Pressure recording at CD
p_FDTD_raw_181 = p_FDTD_raw(:,4);    % Pressure recording at CE
p_FDTD_raw_182 = p_FDTD_raw(:,5);    % Pressure recording at CF
p_FDTD_raw_183 = p_FDTD_raw(:,6);    % Pressure recording at CG
p_FDTD_raw_184 = p_FDTD_raw(:,7);    % Pressure recording at CH
p_FDTD_raw_185 = p_FDTD_raw(:,8);    % Pressure recording at CI
p_FDTD_raw_186 = p_FDTD_raw(:,9);    % Pressure recording at CJ
p_FDTD_raw_187 = p_FDTD_raw(:,10);   % Pressure recording at WA
p_FDTD_raw_188 = p_FDTD_raw(:,11);   % Pressure recording at WG
p_FDTD_raw_189 = p_FDTD_raw(:,12);   % Pressure recording at WY
p_FDTD_raw_191 = p_FDTD_raw(:,13);   % Pressure recording at WA
p_FDTD_raw_192 = p_FDTD_raw(:,14);   % Pressure recording at WG
p_FDTD_raw_193 = p_FDTD_raw(:,15);   % Pressure recording at WY
p_FDTD_raw_194 = p_FDTD_raw(:,16);   % Pressure recording at WA
p_FDTD_raw_195 = p_FDTD_raw(:,17);   % Pressure recording at WG
p_FDTD_raw_196 = p_FDTD_raw(:,18);   % Pressure recording at WY
p_FDTD_raw_197 = p_FDTD_raw(:,19);   % Pressure recording at WA
p_FDTD_raw_198 = p_FDTD_raw(:,20);   % Pressure recording at WG
p_FDTD_raw_199 = p_FDTD_raw(:,21);   % Pressure recording at WY
p_FDTD_raw_200 = p_FDTD_raw(:,22);   % Pressure recording at WA
p_FDTD_raw_202 = p_FDTD_raw(:,23);   % Pressure recording at WG
p_FDTD_raw_207 = p_FDTD_raw(:,24);   % Pressure recording at WY
p_FDTD_raw_289 = p_FDTD_raw(:,25);   % Pressure recording at WY


% Plot
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2 scrsz(4)/1.8]);
h = axes('FontSize',10);
plot(tt_FDTD_raw,p_FDTD_raw_179,'LineStyle','-','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_180,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_181,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_182,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_183,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_184,'LineStyle','-','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_185,'LineStyle','-','LineWidth',1,'Color',[0.9 0.6 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_186,'LineStyle','-','LineWidth',1,'Color',[0.7 0.0 0.7]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_187,'LineStyle','-','LineWidth',1,'Color',[0.8 0.8 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_188,'LineStyle','-','LineWidth',1,'Color',[0.0 0.7 0.9]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_189,'LineStyle','-','LineWidth',1,'Color',[0.3 0.3 0.3]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_191,'LineStyle',':','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_192,'LineStyle',':','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_193,'LineStyle',':','LineWidth',1,'Color',[1.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_194,'LineStyle',':','LineWidth',1,'Color',[1.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_195,'LineStyle',':','LineWidth',1,'Color',[0.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_196,'LineStyle',':','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_197,'LineStyle',':','LineWidth',1,'Color',[0.9 0.6 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_198,'LineStyle',':','LineWidth',1,'Color',[0.7 0.0 0.7]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_199,'LineStyle',':','LineWidth',1,'Color',[0.8 0.8 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_200,'LineStyle',':','LineWidth',1,'Color',[0.0 0.7 0.9]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_202,'LineStyle',':','LineWidth',1,'Color',[0.3 0.3 0.3]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_207,'LineStyle','--','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_289,'LineStyle','--','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;

title('Time history of pressure from FDTD simulation', 'fontsize', 18);
xlabel('Time [s]','fontsize', 18); 
ylabel('Pressure [Pa]','fontsize', 18);
legend('179','180','181','182','183','184','185','186','187','188','189','191','192','193','194','195','196','197','198','199','200','202','207','289');
axis([0 tt_FDTD_raw(end) -0.8 1.4]);
set(gcf, 'Color', 'w');
export_fig p_mic.tif -a1 -r600 -painters

% Plot
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2 scrsz(4)/1.8]);
h = axes('FontSize',10);
plot(tt_FDTD_raw,p_FDTD_raw_179,'LineStyle','-','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_180,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_181,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_182,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_183,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 1.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_184,'LineStyle','-','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_196,'LineStyle',':','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_197,'LineStyle',':','LineWidth',1,'Color',[0.9 0.6 0.0]), hold on;
plot(tt_FDTD_raw,p_FDTD_raw_198,'LineStyle',':','LineWidth',1,'Color',[0.7 0.0 0.7]), hold on;
title('Time history of pressure from FDTD simulation', 'fontsize', 18);
xlabel('Time [s]','fontsize', 18); 
ylabel('Pressure [Pa]','fontsize', 18);
legend('179','180','181','182','183','184','196','197','198');
axis([0 tt_FDTD_raw(end) -0.8 1.4]);
set(gcf, 'Color', 'w');
export_fig p_mic_select.tif -a1 -r600 -painters

% %% Interpolation (To match sampling rate across all data)
% fs = fs_NASA_C;
% dt = 1/fs;
% T = 0.5;
% N = T/dt;
% tt = (0:N)*dt;
% lead_time = 0.05;    % Amount of leading time before front shock arrival
% 
% tt_NASA_interp_WA = (tt_NASA_raw_WA(1)/dt:tt_NASA_raw_WA(end)/dt)*dt;
% tt_NASA_interp_WB = (tt_NASA_raw_WB(1)/dt:tt_NASA_raw_WB(end)/dt)*dt;
% tt_NASA_interp_WC = (tt_NASA_raw_WC(1)/dt:tt_NASA_raw_WC(end)/dt)*dt;
% tt_NASA_interp_WD = (tt_NASA_raw_WD(1)/dt:tt_NASA_raw_WD(end)/dt)*dt;
% tt_NASA_interp_WE = (tt_NASA_raw_WE(1)/dt:tt_NASA_raw_WE(end)/dt)*dt;
% tt_NASA_interp_WF = (tt_NASA_raw_WF(1)/dt:tt_NASA_raw_WF(end)/dt)*dt;
% tt_NASA_interp_WG = (tt_NASA_raw_WG(1)/dt:tt_NASA_raw_WG(end)/dt)*dt;
% tt_NASA_interp_WH = (tt_NASA_raw_WH(1)/dt:tt_NASA_raw_WH(end)/dt)*dt;
% tt_NASA_interp_WY = (tt_NASA_raw_WY(1)/dt:tt_NASA_raw_WY(end)/dt)*dt;
% tt_NASA_interp_WZ = (tt_NASA_raw_WZ(1)/dt:tt_NASA_raw_WZ(end)/dt)*dt;
% % p_NASA_interp_WA = interp1(tt_NASA_raw_WA,p_NASA_raw_WA,tt_NASA_interp_WA);    % Pressure recording at WA
% p_NASA_interp_WB = interp1(tt_NASA_raw_WB,p_NASA_raw_WB,tt_NASA_interp_WB);    % Pressure recording at WB
% % p_NASA_interp_WC = interp1(tt_NASA_raw_WC,p_NASA_raw_WC,tt_NASA_interp_WC);    % Pressure recording at WC
% % p_NASA_interp_WD = interp1(tt_NASA_raw_WD,p_NASA_raw_WD,tt_NASA_interp_WD);    % Pressure recording at WD
% % p_NASA_interp_WE = interp1(tt_NASA_raw_WE,p_NASA_raw_WE,tt_NASA_interp_WE);    % Pressure recording at WE
% % p_NASA_interp_WF = interp1(tt_NASA_raw_WF,p_NASA_raw_WF,tt_NASA_interp_WF);    % Pressure recording at WF
% p_NASA_interp_WG = interp1(tt_NASA_raw_WG,p_NASA_raw_WG,tt_NASA_interp_WG);    % Pressure recording at WG
% % p_NASA_interp_WH = interp1(tt_NASA_raw_WH,p_NASA_raw_WG,tt_NASA_interp_WH);    % Pressure recording at WH
% p_NASA_interp_WY = interp1(tt_NASA_raw_WY,p_NASA_raw_WY,tt_NASA_interp_WY);    % Pressure recording at WY
% % p_NASA_interp_WZ = interp1(tt_NASA_raw_WZ,p_NASA_raw_WZ,tt_NASA_interp_WZ);    % Pressure recording at WZ
% 
% tt_FDTD_interp = (tt_FDTD_raw(1)/dt:tt_FDTD_raw(end)/dt)*dt;
% p_FDTD_interp_CC = interp1(tt_FDTD_raw,p_FDTD_raw_CC,tt_FDTD_interp);  % Pressure recording at CC
% p_FDTD_interp_CD = interp1(tt_FDTD_raw,p_FDTD_raw_CD,tt_FDTD_interp);  % Pressure recording at CD
% p_FDTD_interp_CE = interp1(tt_FDTD_raw,p_FDTD_raw_CE,tt_FDTD_interp);  % Pressure recording at CE
% p_FDTD_interp_CF = interp1(tt_FDTD_raw,p_FDTD_raw_CF,tt_FDTD_interp);  % Pressure recording at CF
% p_FDTD_interp_CG = interp1(tt_FDTD_raw,p_FDTD_raw_CG,tt_FDTD_interp);  % Pressure recording at CG
% p_FDTD_interp_CH = interp1(tt_FDTD_raw,p_FDTD_raw_CH,tt_FDTD_interp);  % Pressure recording at CH
% p_FDTD_interp_CI = interp1(tt_FDTD_raw,p_FDTD_raw_CI,tt_FDTD_interp);  % Pressure recording at CI
% p_FDTD_interp_CJ = interp1(tt_FDTD_raw,p_FDTD_raw_CJ,tt_FDTD_interp);  % Pressure recording at CJ
% p_FDTD_interp_WB = interp1(tt_FDTD_raw,p_FDTD_raw_WB,tt_FDTD_interp);  % Pressure recording at WB
% p_FDTD_interp_WG = interp1(tt_FDTD_raw,p_FDTD_raw_WG,tt_FDTD_interp);  % Pressure recording at WG
% p_FDTD_interp_WY = interp1(tt_FDTD_raw,p_FDTD_raw_WY,tt_FDTD_interp);  % Pressure recording at WY
% 
% 
% %% Time Delay Calculation
% p_threshold = 5;    % The pressure threshold for detecting front shock arrival for NASA and FDTD CC microphones
% inds_NASA = find(p_NASA_raw_CC>p_threshold);
% inds_FDTD = find(p_FDTD_interp_CC>p_threshold);
% arrv_ind_NASA = inds_NASA(1);
% arrv_ind_FDTD = inds_FDTD(1);
% arrv_time_NASA = arrv_ind_NASA*dt+tt_NASA_raw_C(1);
% arrv_time_FDTD = arrv_ind_FDTD*dt+tt_FDTD_interp(1);
% 
% 
% %% Time Alignment of NASA & FDTD Data
% i_start_NASA_C = arrv_ind_NASA-round(lead_time*fs);
% i_end_NASA_C = i_start_NASA_C+N;
% p_NASA_CA = p_NASA_raw_CA(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CA
% p_NASA_LB = p_NASA_raw_LB(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at LB
% p_NASA_CC = p_NASA_raw_CC(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CC
% p_NASA_CD = p_NASA_raw_CD(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CD
% p_NASA_CE = p_NASA_raw_CE(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CE
% p_NASA_CF = p_NASA_raw_CF(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CF
% p_NASA_CG = p_NASA_raw_CG(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CG
% p_NASA_CH = p_NASA_raw_CH(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CH
% p_NASA_CI = p_NASA_raw_CI(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CI
% p_NASA_CJ = p_NASA_raw_CJ(i_start_NASA_C:i_end_NASA_C);         % NASA measurement at CJ
% 
% i_start_NASA_WA = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WA(1))*fs);
% i_start_NASA_WB = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WB(1))*fs);
% i_start_NASA_WC = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WC(1))*fs);
% i_start_NASA_WD = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WD(1))*fs);
% i_start_NASA_WE = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WE(1))*fs);
% i_start_NASA_WF = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WF(1))*fs);
% i_start_NASA_WG = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WG(1))*fs);
% i_start_NASA_WH = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WH(1))*fs);
% i_start_NASA_WY = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WY(1))*fs);
% i_start_NASA_WZ = i_start_NASA_C+round((tt_NASA_raw_C(1)-tt_NASA_raw_WZ(1))*fs);
% i_end_NASA_WA = i_start_NASA_WA+N;
% i_end_NASA_WB = i_start_NASA_WB+N;
% i_end_NASA_WC = i_start_NASA_WC+N;
% i_end_NASA_WD = i_start_NASA_WD+N;
% i_end_NASA_WE = i_start_NASA_WE+N;
% i_end_NASA_WF = i_start_NASA_WF+N;
% i_end_NASA_WG = i_start_NASA_WG+N;
% i_end_NASA_WH = i_start_NASA_WH+N;
% i_end_NASA_WY = i_start_NASA_WY+N;
% i_end_NASA_WZ = i_start_NASA_WZ+N;
% % p_NASA_WA = p_NASA_interp_WA(i_start_NASA_WA:i_end_NASA_WA);    % NASA measurement at WA
% p_NASA_WB = p_NASA_interp_WB(i_start_NASA_WB:i_end_NASA_WB);    % NASA measurement at WB
% % p_NASA_WC = p_NASA_interp_WC(i_start_NASA_WC:i_end_NASA_WC);    % NASA measurement at WC
% % p_NASA_WD = p_NASA_interp_WD(i_start_NASA_WD:i_end_NASA_WD);    % NASA measurement at WD
% % p_NASA_WE = p_NASA_interp_WE(i_start_NASA_WE:i_end_NASA_WE);    % NASA measurement at WE
% % p_NASA_WF = p_NASA_interp_WF(i_start_NASA_WF:i_end_NASA_WF);    % NASA measurement at WF
% p_NASA_WG = p_NASA_interp_WG(i_start_NASA_WG:i_end_NASA_WG);    % NASA measurement at WG
% % p_NASA_WH = p_NASA_interp_WH(i_start_NASA_WH:i_end_NASA_WH);    % NASA measurement at WH
% p_NASA_WY = p_NASA_interp_WY(i_start_NASA_WY:i_end_NASA_WY);    % NASA measurement at WY
% % p_NASA_WZ = p_NASA_interp_WZ(i_start_NASA_WZ:i_end_NASA_WZ);    % NASA measurement at WZ
% 
% i_start_FDTD = arrv_ind_FDTD-round(lead_time*fs);
% i_end_FDTD = i_start_FDTD+N;
% p_FDTD_CC = p_FDTD_interp_CC(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CC
% p_FDTD_CD = p_FDTD_interp_CD(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CD
% p_FDTD_CE = p_FDTD_interp_CE(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CE
% p_FDTD_CF = p_FDTD_interp_CF(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CF
% p_FDTD_CG = p_FDTD_interp_CG(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CG
% p_FDTD_CH = p_FDTD_interp_CH(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CH
% p_FDTD_CI = p_FDTD_interp_CI(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CI
% p_FDTD_CJ = p_FDTD_interp_CJ(i_start_FDTD:i_end_FDTD);          % FDTD measurement at CJ
% p_FDTD_WB = p_FDTD_interp_WB(i_start_FDTD:i_end_FDTD);          % FDTD measurement at WB
% p_FDTD_WG = p_FDTD_interp_WG(i_start_FDTD:i_end_FDTD);          % FDTD measurement at WG
% p_FDTD_WY = p_FDTD_interp_WY(i_start_FDTD:i_end_FDTD);          % FDTD measurement at WY
% 
% % Plot
% figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2 scrsz(4)/1.8]);
% h = axes('FontSize',16);
% plot(tt,p_NASA_CC,'LineStyle','-','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
% plot(tt,p_NASA_CD,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;
% plot(tt,p_NASA_CE,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 0.0]), hold on;
% plot(tt,p_NASA_CF,'LineStyle','-','LineWidth',1,'Color',[1.0 0.0 1.0]), hold on;
% plot(tt,p_NASA_CG,'LineStyle','-','LineWidth',1,'Color',[0.0 0.0 1.0]), hold on;
% plot(tt,p_NASA_CH,'LineStyle','-','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
% plot(tt,p_NASA_CI,'LineStyle','-','LineWidth',1,'Color',[0.9 0.6 0.0]), hold on;
% plot(tt,p_NASA_CJ,'LineStyle','-','LineWidth',1,'Color',[0.7 0.0 0.7]), hold on;
% plot(tt,p_NASA_WB,'LineStyle','-','LineWidth',1,'Color',[0.8 0.8 0.2]), hold on;
% plot(tt,p_NASA_WG,'LineStyle','-','LineWidth',1,'Color',[0.0 0.7 0.9]), hold on;
% plot(tt,p_NASA_WY,'LineStyle','-','LineWidth',1,'Color',[0.3 0.3 0.3]), hold on;
% plot(tt,p_FDTD_CC,'LineStyle',':','LineWidth',1,'Color',[0.0 0.7 0.0]), hold on;
% plot(tt,p_FDTD_CD,'LineStyle',':','LineWidth',1,'Color',[0.0 0.0 0.0]), hold on;
% plot(tt,p_FDTD_CE,'LineStyle',':','LineWidth',1,'Color',[1.0 0.0 0.0]), hold on;
% plot(tt,p_FDTD_CF,'LineStyle',':','LineWidth',1,'Color',[1.0 0.0 1.0]), hold on;
% plot(tt,p_FDTD_CG,'LineStyle',':','LineWidth',1,'Color',[0.0 0.0 1.0]), hold on;
% plot(tt,p_FDTD_CH,'LineStyle',':','LineWidth',1,'Color',[0.7 0.9 0.2]), hold on;
% plot(tt,p_FDTD_CI,'LineStyle',':','LineWidth',1,'Color',[0.9 0.6 0.0]), hold on;
% plot(tt,p_FDTD_CJ,'LineStyle',':','LineWidth',1,'Color',[0.7 0.0 0.7]), hold on;
% plot(tt,p_FDTD_WB,'LineStyle',':','LineWidth',1,'Color',[0.8 0.8 0.2]), hold on;
% plot(tt,p_FDTD_WG,'LineStyle',':','LineWidth',1,'Color',[0.0 0.7 0.9]), hold on;
% plot(tt,p_FDTD_WY,'LineStyle',':','LineWidth',1,'Color',[0.5 0.5 0.5]), hold on;
% title('Comparison between NASA & FDTD data', 'fontsize', 18);
% xlabel('Time [s]','fontsize', 18); 
% ylabel('Pressure [Pa]','fontsize', 18);
% legend('CC','CD','CE','CF','CG','CH','CI','CJ','WB','WG','WY');
% set(gcf, 'Color', 'w');
% export_fig p_mic_comparison.tif -a1 -r600 -painters
% 
close all;
