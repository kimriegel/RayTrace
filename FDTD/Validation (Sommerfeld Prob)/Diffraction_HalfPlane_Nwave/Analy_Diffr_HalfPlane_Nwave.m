function Analy_Diffr_HalfPlane_Nwave(theta0_deg)
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
%   separate figures for every time step.
%
%   *** NOTE: The integral expression fails to produce analytic solution at
%   various locations for Theta0 values being multiples of pi/4. ***

clear all, close all;
bluered=[linspace(0,1,100) ones(1,100); linspace(0,1,100) linspace(1,0,100); ones(1,100) linspace(1,0,100)]';

fig1 = figure('Position',[0 500 550 500]);
fig2 = figure('Position',[550 500 550 500]);
fig3 = figure('Position',[1100 500 550 500]);
fig4 = figure('Position',[0 0 550 500]);
fig5 = figure('Position',[550 0 550 500]);
fig6 = figure('Position',[1100 0 550 500]);

%% Parameters
% Time parameters
dt = 0.02;
t_start = -0.2;
t_end = 0.8;
c = 343;

% Domain parameters
dx = 1.0;
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

N = round((t_end-t_start)/dt);
tt = (0:N-1)*dt+t_start;

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
for n = 1:N
    P_inc = Heaviside(pi-Theta1).*Heaviside(tt(n)+R.*cos(Theta1)/c).*Signature(tt(n)+R.*cos(Theta1)/c);
    P_ref = Heaviside(pi-Theta2).*Heaviside(tt(n)+R.*cos(Theta2)/c).*Signature(tt(n)+R.*cos(Theta2)/c);

    for i = 1:length(xx)
        for j = 1:length(yy)
            if tt(n)-R(i,j)/c >= 0  % replacing Heaviside(tt(n)-R(i,j)/c)
                Integrand1 = @(S)Signature(tt(n)-R(i,j)/c-S)./(sqrt(S).*(S+R(i,j)/c*(1+cos(Theta1(i,j)))));
                Integrand2 = @(S)Signature(tt(n)-R(i,j)/c-S)./(sqrt(S).*(S+R(i,j)/c*(1+cos(Theta2(i,j)))));
                P_dif1(i,j) = -sign(pi-Theta1(i,j))*sqrt(R(i,j)/c*(1+cos(Theta1(i,j))))/(2*pi)*Heaviside(tt(n)-R(i,j)/c)*quadgk(Integrand1,0,tt(n)-R(i,j)/c,'RelTol',reltol);
                P_dif2(i,j) = -sign(pi-Theta2(i,j))*sqrt(R(i,j)/c*(1+cos(Theta2(i,j))))/(2*pi)*Heaviside(tt(n)-R(i,j)/c)*quadgk(Integrand2,0,tt(n)-R(i,j)/c,'RelTol',reltol);
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

    figure(fig1);
    pcolor(X,Y,P_inc);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Incident, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_inc(n) = getframe;

    figure(fig2);
    pcolor(X,Y,P_ref);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Reflected, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_ref(n) = getframe;
    
    figure(fig3);
    pcolor(X,Y,P_inc+P_ref);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Incident + Reflected, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-2.0 2.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_incref(n) = getframe;
    
    figure(fig4);
    pcolor(X,Y,P_dif1);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Diffracted from incident, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_dif1(n) = getframe;

    figure(fig5);
    pcolor(X,Y,P_dif2);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Diffracted from reflected, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_dif2(n) = getframe;
    
    figure(fig6);
    pcolor(X,Y,P_tot);
    xlabel('x [m]'), ylabel('y [m]');
    title(['Total field, \theta = ' num2str(theta0_deg) '\circ, t = ' num2str(tt(n)) 's']);
    colormap(bluered);
    shading flat;
    caxis([-2.0 2.0]);
    colorbar;
    axis image;
    line([0 xx(end)],[0 0],[0 0],'Color','k','LineWidth',2)
    M_tot(n) = getframe;
end

%% Animation
mkdir('Animations');
movie2avi(M_inc,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_inc.avi'],'compression','Cinepak','fps',5);
movie2avi(M_ref,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_ref.avi'],'compression','Cinepak','fps',5);
movie2avi(M_incref,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_incref.avi'],'compression','Cinepak','fps',5);
movie2avi(M_dif1,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_dif1.avi'],'compression','Cinepak','fps',5);
movie2avi(M_dif2,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_dif2.avi'],'compression','Cinepak','fps',5);
movie2avi(M_tot,['.\Animations\Analy_Diffr_HalfPlane_Nwave' num2str(theta0_deg) '_tot.avi'],'compression','Cinepak','fps',5);

close all;
display(['Finished executing Analy_Diffr_HalfPlane_Nwave.m for \theta = ' num2str(theta0_deg)]);

%% Functions
function Theta = atan2_mod(X,Y)
    Theta = atan2(X,Y);
    I = find(Theta<0);
    Theta(I) = Theta(I)+2*pi;
end

function X = Heaviside(x)
    X = ones(size(x));
    I = find(x<0);
    X(I) = 0;
end
function res = Signature(S)
    res = Heaviside(S).*Heaviside(N_duration-S).*(N_amplitude-S/N_slope);
end

end