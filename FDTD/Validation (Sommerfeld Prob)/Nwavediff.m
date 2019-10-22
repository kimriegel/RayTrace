function [X,Y,Phi,R] = Nwavediff()

% Diffraction of an arbitrary plane pulse by a rigid half-plane.
% In particular, Sang Ik Cho, Penn State University was interested in the
% N-wave case, where the incident wave is, say, H(t)*H(2-t)(1-t).

reltol=1e-6;
theta0=pi/6;T='PiO2';      %T=theta0
c=1;
x=0.1;
y=0.1;
r=norm([x,y]);
at=atan2(y,x);
if at<0
    theta=2*pi+at;
else
    theta=at;
end
% ---------Compute integral-----------------------------------------------
tstep = 0.01;
T=(-3:tstep:5);
Inc=zeros(length(T));
Ref=zeros(length(T));
Diff=zeros(length(T));
Total=zeros(length(T));

thetap=theta+theta0;
thetam=theta-theta0;

for k=1:length(T)
    t=T(k);

    if t>-r/c*cos(thetam) && pi>=thetam
        Inc(k) = Inc(k) + heaviside(pi-thetam)*g(t+r/c*cos(thetam));
    end
    if t>-r/c*cos(thetap) && pi>=thetap
        Ref(k) = Ref(k) + heaviside(pi-thetap)*g(t+r/c*cos(thetap));
    end
    if t>r/c
        Diff(k)=Diff(k) -...
                sign(pi-thetam)/(2*pi)*sqrt(r/c*(1+cos(thetam)))*quadgk(@(s)integrandm(s),0,t-r/c,'RelTol',reltol)-...
                sign(pi-thetap)/(2*pi)*sqrt(r/c*(1+cos(thetap)))*quadgk(@(s)integrandp(s),0,t-r/c,'RelTol',reltol);
    end
end
   
Total = Inc+Ref+Diff;

function res = heaviside(s)
    if s < 0
        res = 0;
    else
        res =1;
    end
end
function res = g(s)
    res=zeros(1,length(s));
    for l=1:length(s)
        if s(l)<2
        res(l)=1-s(l);        
        end
    end   
end
function res=integrandm(s)
    res = zeros(1,length(s));
    for l=1:length(s)
        res(l)=g(t-r/c -s(l))/(sqrt(s(l))*(s(l)+r/c*(1+cos(thetam))));
    end        
end
function res=integrandp(s)
    res = zeros(1,length(s));
    for l=1:length(s)
        res(l)=g(t-r/c -s(l))/(sqrt(s(l))*(s(l)+r/c*(1+cos(thetap))));
    end        
end
% ------------Plot results--------------------------------------------------
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2-100 scrsz(3) scrsz(4)/2])
subplot(1,3,1);plot(T,Total), axis([min(T) max(T) -2 2]);
title(['Total field, x=' num2str(x), ', y=' num2str(y)])
subplot(1,3,2);plot(T,Inc+Ref), axis([min(T) max(T) -2 2]);
title(['Inc+Ref fields, x=' num2str(x), ', y=' num2str(y)])
subplot(1,3,3);plot(T,Diff), axis([min(T) max(T) -2 2]);
title(['Diffracted field, x=' num2str(x), ', y=' num2str(y)])
end