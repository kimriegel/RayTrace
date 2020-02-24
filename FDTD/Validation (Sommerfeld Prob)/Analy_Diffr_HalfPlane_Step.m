function Analy_Diffr_HalfPlane_Step()
%% Analy_Diffr_HalfPlane_Step.m
%
%   written by Sang-Ik Terry Cho
%   last modified on 12 / 20 / 2011
%
%       This program plots an analytical closed-form solution of a 2-D 
%   half plane diffraction problem found in Hewett (2011).  The solution
%   is for an incoming plane wave of Dirac delta function and consists 
%   only of elementary functions.
%       The example case presented in this program models a rigid half
%   plane, which corresponds to dp/dn = 0 on the boundary.  The half plane
%   lies on the (+)ve x-axis starting at the origin and the radial angle is
%   defined going counterclockwise from the (+)ve x-axis.
%       The total pressure and each of its components are visualized in
%   separate figures for every time step.

clear all, close all;
bluered=[linspace(0,1,25) ones(1,25); linspace(0,1,25) linspace(1,0,25); ones(1,25) linspace(1,0,25)]';

fig1 = figure('Position',[0 500 500 500]);
fig2 = figure('Position',[500 500 500 500]);
fig3 = figure('Position',[1000 500 500 500]);
fig4 = figure('Position',[0 0 500 500]);
fig5 = figure('Position',[500 0 500 500]);
fig6 = figure('Position',[1000 0 500 500]);

%% Domain
xx = (-150:150);
yy = (-150:150);

N = 100;
dt = 0.5;
tt = (0:N-1)*dt;

[X,Y] = ndgrid(xx,yy);
Theta = atan2_mod(Y,X);
R = sqrt(X.^2+Y.^2);

%% Parameter
c = 1;
theta0 = pi/8;
Theta1 = Theta-theta0;
Theta2 = Theta+theta0;

%% Propagation
for n = 1:N
    P_inc = Heaviside(pi-Theta1).*Heaviside(tt(n)+R.*cos(Theta1)/c);
    P_ref = Heaviside(pi-Theta2).*Heaviside(tt(n)+R.*cos(Theta2)/c);
    P_dif1 = -sign(pi-Theta1).*Heaviside(tt(n)-R/c)/pi.*atan(sqrt((c*tt(n)-R)./(R.*(1+cos(Theta1)))));
    P_dif2 = -sign(pi-Theta2).*Heaviside(tt(n)-R/c)/pi.*atan(sqrt((c*tt(n)-R)./(R.*(1+cos(Theta2)))));
    P_step = P_inc+P_ref+P_dif1+P_dif2;
    
    figure(fig1);
    pcolor(X,Y,P_inc);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
%     colorbar;
    axis image;
    
    figure(fig2);
    pcolor(X,Y,P_ref);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
%     colorbar;
    axis image;

    figure(fig3);
    pcolor(X,Y,P_inc+P_ref);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-2.0 2.0]);
%     colorbar;
    axis image;

    figure(fig4);
    pcolor(X,Y,P_dif1);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
%     colorbar;
    axis image;
    
    figure(fig5);
    pcolor(X,Y,P_dif2);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-1.0 1.0]);
%     colorbar;
    axis image;

    figure(fig6);
    surf(X,Y,10*P_step);
    xlabel('x'), ylabel('y');
    colormap(bluered);
    shading flat;
    caxis([-20.0 20.0]);
%     colorbar;
    axis image;
    axis off
end

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

end