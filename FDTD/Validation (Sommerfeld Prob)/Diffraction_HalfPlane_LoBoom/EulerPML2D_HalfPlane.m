function dU = EulerPML2D_HalfPlane(U,CFL,z0,origin_i,origin_j)
%%  EulerPML2D_HalfPlane.m
%
%   written by Sang Cho
%   last modified on 1 / 16 / 12
%
%       This function takes as an input a structure containing the 
%   pressure and velocity fields and computes the divergence of flux term 
%   in the Euler equations, which is to be used as an argument for Runge 
%   Kutta time integration scheme.  The domain boundaries are rigid.
%       A rigid half plane is included which extends in (+)ve x-direction
%   from the origin location passed to the function.
%   

dim = size(U.px);
I_MAX = dim(1);
J_MAX = dim(2);

dU.px = zeros(I_MAX,J_MAX);
dU.py = zeros(I_MAX,J_MAX);
dU.u = zeros(I_MAX,J_MAX);
dU.v = zeros(I_MAX,J_MAX);
U.p = U.px+U.py;

for j = 1:J_MAX
    for i=1:I_MAX
        % x-components
        if i==1
            dU.px(i,j) = -CFL/2*z0*(U.u(i+1,j)+U.u(i,j));
            dU.u(i,j) = -CFL/2/z0*(U.p(i+1,j)-U.p(i,j));
        elseif i==I_MAX
            dU.px(i,j) = -CFL/2*z0*(-U.u(i,j)-U.u(i-1,j));
            dU.u(i,j) = -CFL/2/z0*(U.p(i,j)-U.p(i-1,j));
        else
            dU.px(i,j) = -CFL/2*z0*(U.u(i+1,j)-U.u(i-1,j));
            dU.u(i,j) = -CFL/2/z0*(U.p(i+1,j)-U.p(i-1,j));
        end
        % y-components
        if (j==origin_j+1)&&(i>origin_i+1)
            dU.py(i,j) = -CFL/2*z0*(U.v(i,j+1)+U.v(i,j));
            dU.v(i,j) = -CFL/2/z0*(U.p(i,j+1)-U.p(i,j));
        elseif (j==origin_j)&&(i>origin_i+1)
            dU.py(i,j) = -CFL/2*z0*(-U.v(i,j)-U.v(i,j-1));
            dU.v(i,j) = -CFL/2/z0*(U.p(i,j)-U.p(i,j-1));            
        elseif j==1
            dU.py(i,j) = -CFL/2*z0*(U.v(i,j+1)+U.v(i,j));
            dU.v(i,j) = -CFL/2/z0*(U.p(i,j+1)-U.p(i,j));
        elseif j==J_MAX
            dU.py(i,j) = -CFL/2*z0*(-U.v(i,j)-U.v(i,j-1));
            dU.v(i,j) = -CFL/2/z0*(U.p(i,j)-U.p(i,j-1));
        else
            dU.py(i,j) = -CFL/2*z0*(U.v(i,j+1)-U.v(i,j-1));
            dU.v(i,j) = -CFL/2/z0*(U.p(i,j+1)-U.p(i,j-1));
        end
    end
end