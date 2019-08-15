%clear all
format short
% Read in the data from the Ray Trace Code
fid=fopen('~/RayTrace/PythonTest.txt');
%A=fscanf(fid,'%f',[7,inf]);
A=ReadTecPlotWorks(fid,'Single Point');
Aplus=A';
% Sort the array first by the x direction, then the Y direction, Then the
% time of arrival.  This allows us to group together arrays by point and
% then by time.
Asorted=sortrows(Aplus,[2,3,7]);

x=Asorted(1,2);
y=Asorted(1,3);
Grid(1,:)=[x y];
mod=0;
% loop through the array and create a new array with a dimension for each
% point in space.  
j=1;
m=2;
for n=1:length(Asorted)
    if ((Asorted(n,2)==x)&&(Asorted(n,3)==y))
         iter=n-mod;
         B(iter,:,j)=Asorted(n,:);
    else
        mod=n-1;
        x=Asorted(n,2);
        y=Asorted(n,3);
        Grid(m,:)=[x y];
        iter=n-mod;
        j=j+1;
        n=n+1;
        B(iter,:,j)=Asorted(n,:);
    end
end
[l,e,loc]=size(B);
mod=0;
iter=1;
k=1;
% loop through the new array and, for each point in space, create 
%a new dimension, for each time of arrival.
for j=1:loc
    mod=0;
    time=B(1,7,j);
    for m=1:l
        if (B(m,7,j)==time)
            iter=m-mod;
            C(iter,:,j,k)=B(iter,:,j);
        else
            mod=m-1;
            time=B(m,7,j);
            if(time~=0)               
                iter=m-mod;
                k=k+1;
                C(iter,:,j,k)=B(iter,:,j);
            end
        end
    end
end
% create an array that has just the fft for each point at each arrival
% time.
for j=1:loc
    for p=1:l+1
        m=2*l+1-p;
    if p==1
        new1(p,5,j)=0;
        new1(p,1,j)=0;
        new1(p,2:4,j)=C(p,2:4,j);
        new1(p,6,j)=C(p,7,j);
    else
        new1(p,5,j)=C(p-1,6,j).*exp(i.*C(p-1,5,j));
        new1(p,1,j)=C(p-1,1,j);
        new1(p,2:4,j)=C(p-1,2:4,j);
        new1(p,6,j)=C(p-1,7,j);
    end
    if(p<(l))
        new1(m,5,j)=conj(C(p,6,j).*exp(i.*C(p,5,j)));
        new1(m,1,j)=C(p,1,j);
        new1(m,2:4,j)=C(p,2:4,j);
        new1(m,6,j)=C(p,7,j);
    end
    end
end
% plot the inverse fft for a location of my chosing.  Displace the initial
% time by the arrival time.
T=.00005;
Fs=1/T;
L=1024;
t=((0:L-1)*T)+new1(1,6,1);
figure 
plot(t,(ifft(new1(1:1024,5,1)/2)*((2*l))))  
    