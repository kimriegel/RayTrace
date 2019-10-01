%% Analy_Diffr_HalfPlane_Nwave.m
%
%   written by Sang-Ik Terry Cho
%   last modified on 1 / 11 / 2012
%
%       This program plots an analytical closed-form solution of a 2-D 
%   half plane diffraction problem found in Hewett (2011).  The solution
%   is for an incoming plane wave of a perfect N-wave.
%       The example case presented in this program models a rigid half
%   plane, which corresponds to dp/dn = 0 on the boundary.  The half plane
%   lies on the (+)ve x-axis starting at the origin and the radial angle is
%   defined going counterclockwise from the (+)ve x-axis.
%       The total pressure and each of its components are visualized in
%   separate figures at a user-specified time.
%
%   *** NOTE: The integral expression fails to produce analytic solution at
%   locations for which the Theta0 value becomes multiples of pi/4. ***

clear all, close all;
bluered=[linspace(0,1,100) ones(1,100); linspace(0,1,100) linspace(1,0,100); ones(1,100) linspace(1,0,100)]';

fig1 = figure('Position',[0 500 550 500]);
fig2 = figure('Position',[550 500 550 500]);
fig3 = figure('Position',[1100 500 600 500]);
fig4 = figure('Position',[0 50 1400 500]);

%% Parameters
% Time parameters
t_snapshot = 0.2;   % Time the snapshot of the field is to be taken
c = 343;

% Domain parameters
dx = 0.5;
x_min = -100;
x_max = 100;
y_min = -100;
y_max = 100;

% Incoming signature parameters
N_duration = 0.3;
N_amplitude = 1;
N_slope = N_duration/2/N_amplitude;

% theta0 = 2*pi/3;
% theta0_deg = theta0/pi*180;
theta0_deg = 30;
theta0 = theta0_deg/180*pi;

% Simulatiom parameter
reltol = 1e-2;

%% Preparing the simulation
xx = (x_min/dx:x_max/dx)*dx;
yy = (y_min/dx:y_max/dx)*dx;

[X,Y] = ndgrid(xx,yy);
Theta = atan2_mod(Y,X);
Theta1 = Theta-theta0;
Theta2 = Theta+theta0;
R = sqrt(X.^2+Y.^2);

P_inc = zeros(size(R));
P_ref = zeros(size(R));
P_dif1 = zeros(size(R));
P_dif2 = zeros(size(R));

%% Propagation
P_inc = Heaviside(pi-Theta1).*Heaviside(t_snapshot+R.*cos(Theta1)/c).*Signature(t_snapshot+R.*cos(Theta1)/c,N_duration,N_amplitude,N_slope);
P_ref = Heaviside(pi-Theta2).*Heaviside(t_snapshot+R.*cos(Theta2)/c).*Signature(t_snapshot+R.*cos(Theta2)/c,N_duration,N_amplitude,N_slope);

for i = 1:length(xx)
    for j = 1:length(yy)
        if t_snapshot-R(i,j)/c >= 0  % replacing Heaviside(t-R(i,j)/c)
            Integrand1 = @(S)Signature(t_snapshot-R(i,j)/c-S,N_duration,N_amplitude,N_slope)./(sqrt(S).*(S+R(i,j)/c*(1+cos(Theta1(i,j)))));
            Integrand2 = @(S)Signature(t_snapshot-R(i,j)/c-S,N_duration,N_amplitude,N_slope)./(sqrt(S).*(S+R(i,j)/c*(1+cos(Theta2(i,j)))));
            P_dif1(i,j) = -sign(pi-Theta1(i,j))*sqrt(R(i,j)/c*(1+cos(Theta1(i,j))))/(2*pi)*Heaviside(t_snapshot-R(i,j)/c)*quadgk(Integrand1,0,t_snapshot-R(i,j)/c,'RelTol',reltol);
            P_dif2(i,j) = -sign(pi-Theta2(i,j))*sqrt(R(i,j)/c*(1+cos(Theta2(i,j))))/(2*pi)*Heaviside(t_snapshot-R(i,j)/c)*quadgk(Integrand2,0,t_snapshot-R(i,j)/c,'RelTol',reltol);
        end
    end
end

%     % Take an average value at the neighboring grid points if the integration produced 'inf'
%     I1 = ~isfinite(P_dif1);
%     I1(1) = 0;
%     I2 = ~isfinite(P_dif2);
%     I2(1) = 0;
%     P_dif1(I1)=(P_dif1(find(I1)+1)+P_dif1(find(I1)-1))/2;
%     P_dif2(I2)=(P_dif2(find(I2)+1)+P_dif2(find(I2)-1))/2;

P_tot = P_inc+P_ref+P_dif1+P_dif2;

mkdir('Snapshots');
cd('Snapshots');
figure(fig1);
h = axes('FontSize',16);
pcolor(X,Y,P_inc+P_ref);
% title('(a) Incident + Reflected','fontsize',18);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
% colorbar('FontSize',16);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2);
set(gcf, 'Color', 'w');
export_fig(sprintf('Snapshot_Analy_IncRef_%g.png',theta0_deg),'-a1','-r300','-painters');

figure(fig2);
h = axes('FontSize',16);
pcolor(X,Y,P_dif1+P_dif2);
% title('(b) Diffracted','fontsize',18);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
% colorbar('FontSize',16);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
set(gcf, 'Color', 'w');
export_fig(sprintf('Snapshot_Analy_Diffr_%g.png',theta0_deg),'-a1','-r300','-painters');

figure(fig3);
h = axes('FontSize',16);
pcolor(X,Y,P_tot);
% title('(c) Total','fontsize',18);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
colorbar('FontSize',16);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2);
set(gcf, 'Color', 'w');
export_fig(sprintf('Snapshot_Analy_Tot_%g.png',theta0_deg),'-a1','-r300','-painters');

figure(fig4);
subplot(1,3,1), pcolor(X,Y,P_inc+P_ref);
title('(a) Incident + Reflected','fontsize',14);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
colorbar('FontSize',12);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2);
set(gcf, 'Color', 'w');

subplot(1,3,2), pcolor(X,Y,P_dif1+P_dif2);
title('(b) Diffracted','fontsize',14);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
colorbar('FontSize',12);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2);
set(gcf, 'Color', 'w');

subplot(1,3,3), pcolor(X,Y,P_tot);
title('(c) Total','fontsize',14);
colormap(bluered);
shading flat;
caxis([-2.0 2.0]);
colorbar('FontSize',12);
axis image;
line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2);
set(gcf, 'Color', 'w');
export_fig(sprintf('Snapshot_Analy_%g.png',theta0_deg),'-a1','-r300','-painters');

cd('..');