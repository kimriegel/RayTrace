function [] = VerifyHybrid(method,p_hybrid,fs_hybrid,p_NASA,fs_NASA)
%% VerifyHybrid.m
%
%   written by:     Sang Ik Cho
%   last modifed:   7 / 22 / 2011
%   e-mail:         stc142@psu.edu
%
%       This script generates pressure time series plots comparing the
%   FDTD-ray tracing hybrid result with experimental data.  For this visual 
%   comparison, the shock arrival times of the simulation and measurement
%   data must be performed, in the similar manner used for generating the 
%   hybrid result.  A choice from three different methods is given:
%       Method 1: Find the midpoint of initial shock by setting a pressure 
%                   threshold at the half of the peak
%       Method 2: Match the time indices at the peak of derivatives of each
%                   numerical recordings
%       Method 3: Find the midpoint of initial shock by averaging the times
%                   where pressure reaches 10% and 90% of the peak
%   The peak amplitude dB offset of the hybrid result from the measured 
%   data at each microphone is displayed on the screen.


%% Interpolation (To match sampling rate across all data)
dt_hybrid = 1/fs_hybrid;
N_hybrid = length(p_hybrid);
tt_hybrid = (0:N_hybrid-1)*dt_hybrid;
dt_NASA = 1/fs_NASA;
N_NASA = length(p_NASA);
tt_NASA = (0:N_NASA-1)*dt_NASA;

fs = max(fs_hybrid,fs_NASA);

% "Resampling" time vectors with new fs: 
% Multiplying the beginning and end of the time series by the new fs puts 
% every point that needs to be sampled at an integer marking, which then
% can be put in to an integer series by MATLAB. Dividing the new series by
% fs results in a new time vector.
tt_hybrid_interp = (tt_hybrid(1)*fs:tt_hybrid(end)*fs)/fs;
p_hybrid_interp = interp1(tt_hybrid,p_hybrid,tt_hybrid_interp);
tt_NASA_interp = (tt_NASA(1)*fs:tt_NASA(end)*fs)/fs;
p_NASA_interp = interp1(tt_NASA,p_NASA,tt_NASA_interp);


%% Time Alignment of FDTD & Rays Result
p_hybrid_interp_base = p_hybrid_interp(:,1);  % Time history at base microphone location
p_NASA_interp_base = p_NASA_interp(:,1);      % Time history at base microphone location

switch method
    case 1  % Method 1: Find the midpoint of initial shock by setting a pressure threshold at the half of the peak
        threshold = 0.5;    % The threshold for detecting bow shock arrival
        inds_hybrid = find(p_hybrid_interp_base>threshold*max(p_hybrid_interp_base));
        arrv_ind_hybrid = inds_hybrid(1);
        inds_NASA = find(p_NASA_interp_base>threshold*max(p_NASA_interp_base));
        arrv_ind_NASA = inds_NASA(1);

        % The shorter of the bow shock arrival times is assigned to be the 
        % "shock-ind"
        if arrv_ind_hybrid < arrv_ind_NASA
            shock_ind = arrv_ind_hybrid;
        else
            shock_ind = arrv_ind_NASA;
        end

        % Determine the starting index
        i_start_hybrid = arrv_ind_hybrid-shock_ind+1;
        i_start_NASA = arrv_ind_NASA-shock_ind+1;

    case 2  % Method 2: Match the time indices at the peak of derivatives of each numerical recordings
        threshold = 0.1;    % The threshold for detecting bow shock arrival
        inds_hybrid = find(p_hybrid_interp_base>threshold*max(p_hybrid_interp_base));
        arrv_ind_hybrid = inds_hybrid(1);
        p_hybrid_deriv_base = diff(p_hybrid_interp_base)*fs;
        p_NASA_deriv_base = diff(p_NASA_interp_base)*fs;
        max_deriv_ind_hybrid = find(p_hybrid_deriv_base==max(p_hybrid_deriv_base));
        max_deriv_ind_NASA = find(p_NASA_deriv_base==max(p_NASA_deriv_base));
        arrv_ind_NASA = arrv_ind_hybrid+max_deriv_ind_NASA-max_deriv_ind_hybrid;

        % The shorter of the bow shock arrival times is assigned to be the 
        % "shock-ind"
        if arrv_ind_hybrid < arrv_ind_NASA
            shock_ind = arrv_ind_hybrid;
        else
            shock_ind = arrv_ind_NASA;
        end

        % Determine the starting index
        i_start_hybrid = arrv_ind_hybrid-shock_ind+1;
        i_start_NASA = arrv_ind_NASA-shock_ind+1;

    case 3  % Method 3: Find the midpoint of initial shock by averaging the times where pressure reaches 10% and 90% of the peak
        peak10_hybrid = find(p_hybrid_interp_base>0.1*max(p_hybrid_interp_base));
        peak90_hybrid = find(p_hybrid_interp_base>0.9*max(p_hybrid_interp_base));
        arrv_ind_hybrid = round((peak10_hybrid(1)+peak90_hybrid(1))/2);
        peak10_NASA = find(p_NASA_interp_base>0.1*max(p_NASA_interp_base));
        peak90_NASA = find(p_NASA_interp_base>0.9*max(p_NASA_interp_base));
        arrv_ind_NASA = round((peak10_NASA(1)+peak90_NASA(1))/2);

        % The shorter of the bow shock arrival times is assigned to be the 
        % "shock-ind"
        if arrv_ind_hybrid < arrv_ind_NASA
            shock_ind = arrv_ind_hybrid;
        else
            shock_ind = arrv_ind_NASA;
        end

        % Determine the starting index
        i_start_hybrid = arrv_ind_hybrid-shock_ind+1;
        i_start_NASA = arrv_ind_NASA-shock_ind+1;
end

% Choose the shorter of the two simulation outputs to match the lengths
if (length(p_hybrid_interp)-i_start_hybrid)<(length(p_NASA_interp)-i_start_NASA)
    N = length(p_hybrid_interp)-i_start_hybrid;
else
    N = length(p_NASA_interp)-i_start_NASA;
end

T = 1/N;
dt = 1/fs;
tt = (0:N)*dt;

% Determine the ending index
i_end_hybrid = i_start_hybrid+N;
i_end_NASA = i_start_NASA+N;


% Truncate Time-Aligned Portions
p_hybrid_aligned = p_hybrid_interp(i_start_hybrid:i_end_hybrid,:);
p_NASA_aligned = p_NASA_interp(i_start_NASA:i_end_NASA,:);


%% Comparison Plots
num_mics = size(p_hybrid,2);

LineColors = [[0.0 0.7 0.0];...
              [1.0 0.0 0.0];...
              [1.0 0.0 1.0];...
              [0.0 0.0 1.0];...
              [0.7 0.9 0.2];...
              [0.9 0.6 0.0];...
              [0.7 0.0 0.7];...
              [0.8 0.8 0.2];...
              [0.0 0.7 0.9]];

% Experimental data at all microphones
figure('Position',[200 200 800 500]);
h = axes('FontSize',16);
for l = 1:num_mics
    plot(tt,p_NASA_aligned(:,l),'Linewidth',1.5,'Color',LineColors(mod(l-1,9)+1,:),'LineStyle','-'), hold on;
end
p_max = max(max(p_NASA_aligned));
p_min = max(max(-p_NASA_aligned));
axis([0 max(tt) -1.2*p_min 1.2*p_max]);
xlabel('Time [s]','fontsize',16); 
ylabel('Pressure [unit]','fontsize',16);
title('Experimental measurement of pressure time series at microphone locations','fontsize',16);
set(gcf, 'Color', 'w');

% Hybrid result at all microphones
figure('Position',[200 200 800 500]);
h = axes('FontSize',16);
for n = 1:num_mics
    plot(tt,p_hybrid_aligned(:,n),'Linewidth',1.5,'Color',LineColors(mod(n-1,9)+1,:),'LineStyle','-'), hold on;
end
p_max = max(max(p_hybrid_aligned));
p_min = max(max(-p_hybrid_aligned));
axis([0 max(tt) -1.2*p_min 1.2*p_max]);
xlabel('Time [s]','fontsize',16); 
ylabel('Pressure [unit]','fontsize',16);
title('Hybrid simulation result of pressure time series at microphone locations','fontsize',16);
set(gcf, 'Color', 'w');

% Comparison at individual microphones
for n = 1:num_mics
    figure('Position',[200 200 800 500]);
    h = axes('FontSize',16);
    plot(tt,p_NASA_aligned(:,n),'Linewidth',1.5,'Color','r','LineStyle','-'), hold on;
    plot(tt,p_hybrid_aligned(:,n),'Linewidth',1.5,'Color','k','LineStyle','--'), hold on;
    p_max = max(max([p_NASA_aligned p_hybrid_aligned]));
    p_min = max(max(-[p_NASA_aligned p_hybrid_aligned]));
    axis([0 max(tt) -1.2*p_min 1.2*p_max]);
    xlabel('Time [s]','fontsize',16); 
    ylabel('Pressure [unit]','fontsize',16);
    title(['Comparison of measurement and hybrid simulation results at microphone #' num2str(n)],'fontsize',16);
    legend('Measurement','Hybrid Simulation');
    set(gcf, 'Color', 'w');
end


%% Peak Amplitude dB Offset Calculation
for n = 1:num_mics
    p_peak_NASA = max(p_NASA_aligned(:,n));
    p_peak_hybrid = max(p_hybrid_aligned(:,n));
    dB_offset = 20*log10(p_peak_hybrid/p_peak_NASA);
    display(['Peak amplitude dB offset of hybrid result at microphone #' num2str(n) ' = ' num2str(dB_offset) ' dB']);
end