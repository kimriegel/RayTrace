%% p_mics_condition.m
%
%   written by Sang Ik Cho
%   last modified on 5 / 22 / 11
%   
%       This program "cleans up" the instability that appeared at the end
%   of LDDRK4 simulation by replacing the end with tapered random noise...
%   Yea...

clear all, close all;
scrsz = get(0,'ScreenSize');

load ('p_mics.dat');
num_mics = size(p_mics,2)-1;
tt = p_mics(:,1);
mics = p_mics(:,2:num_mics+1);
mics_psf = mics/47.88;

for j=1:num_mics
    ind_instability1 = find(abs(mics_psf(1250:1462,j)-mics_psf(1249:1461,j))>0.02)-50;
    if ~isempty(ind_instability1)
        ind_start_decay(j) = ind_instability1(1)+1250;
        mics_psf(ind_start_decay(j):1462,j) = mics_psf(ind_start_decay(j):-1:2*ind_start_decay(j)-1462,j);
    end
    ind_instability2 = find(~isfinite(mics_psf(1250:end,j)))-5;
    if ~isempty(ind_instability2)
        ind_start_decay(j) = ind_instability2(1)+1250;
        mics_psf(ind_start_decay(j):1462,j) = mics_psf(ind_start_decay(j):-1:2*ind_start_decay(j)-1462,j);
    end
end

p_mics_fixed = [tt mics_psf];
save p_mics_fixed.dat p_mics_fixed -ascii
