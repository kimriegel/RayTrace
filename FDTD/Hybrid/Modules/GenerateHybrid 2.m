function [p_hybrid,fs_hybrid] = GenerateHybrid(method,p_FDTD_raw,fs_FDTD_raw,p_Rays_raw,fs_Rays_raw,f_crossover,f_margin,N_filt)
%% GenerateHybrid.m
%
%   written by:     Sang Ik Cho
%   last modifed:   8 / 16 / 2011
%   e-mail:         stc142@psu.edu
%
%       This script takes the pressure time histories recorded during FDTD 
%   and ray tracing simulations of low-boom propagation and diffraction 
%   around building structures and generates a hybrid result of the two 
%   simulation methods, based on calculation of the shock arrival times and
%   frequency domain filtering.  FilterDesign.m is used to generate FIR 
%   filter coefficients.  A choice from three different methods of aligning
%   the shock arrival times is given:
%       Method 1: Find the midpoint of initial shock by setting a pressure 
%                   threshold at the half of the peak
%       Method 2: Match the time indices at the peak of derivatives of each
%                   numerical recordings
%       Method 3: Find the midpoint of initial shock by averaging the times
%                   where pressure reaches 10% and 90% of the peak


%% Interpolation (To match sampling rate across all data)
dt_FDTD_raw = 1/fs_FDTD_raw;
N_FDTD_raw = length(p_FDTD_raw);
tt_FDTD_raw = (0:N_FDTD_raw-1)*dt_FDTD_raw;
dt_Rays_raw = 1/fs_Rays_raw;
N_Rays_raw = length(p_Rays_raw);
tt_Rays_raw = (0:N_Rays_raw-1)*dt_Rays_raw;

fs_hybrid = max(fs_FDTD_raw,fs_Rays_raw);

% "Resampling" time vectors with new fs: 
% Multiplying the beginning and end of the time series by the new fs puts 
% every point that needs to be sampled at an integer marking, which then
% can be put in to an integer series by MATLAB. Dividing the new series by
% fs results in a new time vector.
tt_FDTD_interp = (tt_FDTD_raw(1)*fs_hybrid:tt_FDTD_raw(end)*fs_hybrid)/fs_hybrid;
p_FDTD_interp = interp1(tt_FDTD_raw,p_FDTD_raw,tt_FDTD_interp);
tt_Rays_interp = (tt_Rays_raw(1)*fs_hybrid:tt_Rays_raw(end)*fs_hybrid)/fs_hybrid;
p_Rays_interp = interp1(tt_Rays_raw,p_Rays_raw,tt_Rays_interp);


%% Filters Design
[b_LP,b_HP] = FilterDesign(fs_hybrid,f_crossover,f_margin,N_filt);


%% Filtering
% LP filtering FDTD result
p_FDTD_LP = filter(b_LP,1,p_FDTD_interp);

% HP filtering ray tracing result
p_Rays_HP = filter(b_HP,1,p_Rays_interp);


%% Time Alignment of FDTD & Rays Result
p_FDTD_interp_base = p_FDTD_interp(:,1);  % Time history at base microphone location
p_Rays_interp_base = p_Rays_interp(:,1);  % Time history at base microphone location

switch method
    case 1  % Method 1: Find the midpoint of initial shock by setting a pressure threshold at the half of the peak
        threshold = 0.5;    % The threshold for detecting bow shock arrival
        inds_FDTD = find(p_FDTD_interp_base>threshold*max(p_FDTD_interp_base));
        arrv_ind_FDTD = inds_FDTD(1);
        inds_Rays = find(p_Rays_interp_base>threshold*max(p_Rays_interp_base));
        arrv_ind_Rays = inds_Rays(1);

        % The shorter of the bow shock arrival times for FDTD and Rays 
        % results is assigned to be the "shock_ind", which is the shock 
        % arrival time in the hybrid result at the base microphone location
        if arrv_ind_FDTD < arrv_ind_Rays
            shock_ind = arrv_ind_FDTD;
        else
            shock_ind = arrv_ind_Rays;
        end

        % Determine the starting index
        i_start_FDTD = arrv_ind_FDTD-shock_ind+N_filt;
        i_start_Rays = arrv_ind_Rays-shock_ind+N_filt;

    case 2  % Method 2: Match the time indices at the peak of derivatives of each numerical recordings
        threshold = 0.1;    % The threshold for detecting bow shock arrival
        inds_FDTD = find(p_FDTD_interp_base>threshold*max(p_FDTD_interp_base));
        arrv_ind_FDTD = inds_FDTD(1);
        p_FDTD_deriv_base = diff(p_FDTD_interp_base)*fs_hybrid;
        p_Rays_deriv_base = diff(p_Rays_interp_base)*fs_hybrid;
        max_deriv_ind_FDTD = find(p_FDTD_deriv_base==max(p_FDTD_deriv_base));
        max_deriv_ind_Rays = find(p_Rays_deriv_base==max(p_Rays_deriv_base));
        arrv_ind_Rays = arrv_ind_FDTD+max_deriv_ind_Rays-max_deriv_ind_FDTD;

        % The shorter of the bow shock arrival times for FDTD and Rays 
        % results is assigned to be the "shock_ind", which is the shock 
        % arrival time in the hybrid result at the base microphone location
        if arrv_ind_FDTD < arrv_ind_Rays
            shock_ind = arrv_ind_FDTD;
        else
            shock_ind = arrv_ind_Rays;
        end

        % Determine the starting index
        i_start_FDTD = arrv_ind_FDTD-shock_ind+N_filt;
        i_start_Rays = arrv_ind_Rays-shock_ind+N_filt;

    case 3  % Method 3: Find the midpoint of initial shock by averaging the times where pressure reaches 10% and 90% of the peak
        peak10_FDTD = find(p_FDTD_interp_base>0.1*max(p_FDTD_interp_base));
        peak90_FDTD = find(p_FDTD_interp_base>0.9*max(p_FDTD_interp_base));
        arrv_ind_FDTD = round((peak10_FDTD(1)+peak90_FDTD(1))/2);
        peak10_Rays = find(p_Rays_interp_base>0.1*max(p_Rays_interp_base));
        peak90_Rays = find(p_Rays_interp_base>0.9*max(p_Rays_interp_base));
        arrv_ind_Rays = round((peak10_Rays(1)+peak90_Rays(1))/2);

        % The shorter of the bow shock arrival times for FDTD and Rays 
        % results is assigned to be the "shock_ind", which is the shock 
        % arrival time in the hybrid result at the base microphone location
        if arrv_ind_FDTD < arrv_ind_Rays
            shock_ind = arrv_ind_FDTD;
        else
            shock_ind = arrv_ind_Rays;
        end

        % Determine the starting index
        i_start_FDTD = arrv_ind_FDTD-shock_ind+N_filt;
        i_start_Rays = arrv_ind_Rays-shock_ind+N_filt;
end

% Choose the shorter of the two simulation outputs to match the lengths
if (length(p_FDTD_LP)-i_start_FDTD)<(length(p_Rays_HP)-i_start_Rays)
    N_hybrid = length(p_FDTD_LP)-i_start_FDTD;
else
    N_hybrid = length(p_Rays_HP)-i_start_Rays;
end

% Determine the ending index
i_end_FDTD = i_start_FDTD+N_hybrid;
i_end_Rays = i_start_Rays+N_hybrid;


%% Generate Hybrid Result
% The hybrid result is generated by superposing the LP filtered FDTD result
% and the HP filtered Rays result after time-aligning them.
p_FDTD_LP_aligned = p_FDTD_LP(i_start_FDTD:i_end_FDTD,:);
p_Rays_HP_aligned = p_Rays_HP(i_start_Rays:i_end_Rays,:);
p_hybrid = p_FDTD_LP_aligned+p_Rays_HP_aligned;
