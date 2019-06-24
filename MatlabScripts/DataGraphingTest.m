inputfile = 'C:\Users\Kory George\Documents\RayTrace\inputNASABoom1.txt'; 
%FORTRAN_TEST_simplefile='C:\Users\Kory George\Documents\RayTrace\MatlabScripts\NASABoom1_0Diff_1Alpha_rec1.dat';
FORTRAN_TEST_simplefile = 'C:\Users\Kory George\Documents\RayTrace\PythonTest.txt';
%FORTRAN_TEST_simplefile='C:\Users\Kory George\Documents\RayTrace\MatlabScripts\FortranTestDiffusion_Receiver4.dat';
[header,FORTRAN_TEST_simple]=hdrload(FORTRAN_TEST_simplefile);
[header,input]=hdrload(inputfile);

Fs=24000;
T=1/Fs;
N=length(input);
NFFT = N;
dt=1/Fs;
tt=((0:N-1)*dt)';
t=(-NFFT+1:NFFT-1)*T;
f_crossover = 40;
f_margin = 5;
f_trans_lo = f_crossover-f_margin;
f_trans_hi = f_crossover+f_margin;
f = [0 f_trans_lo f_trans_hi (Fs/2)]/(Fs/2);

a_LP = [1 1 0 0];
a_HP = [0 0 1 1];
N_filt = 2000;   % The linear phase filter will be of order (2N_filt-1) with
                % a delay of N_filt samples.
freqz_l = 12000;

% Constructing the linear phase filter
b_LP = firpm(2*(N_filt-1),f,a_LP,[1 1]);
[h_LP,ff] = freqz(b_LP,1,freqz_l,Fs);
Gd_LP = grpdelay(b_LP,1,ff,Fs);
b_HP = firpm(2*(N_filt-1),f,a_HP,[1 1]);
[h_HP,ff] = freqz(b_HP,1,freqz_l,Fs);
Gd_HP = grpdelay(b_HP,1,ff,Fs);

%X=0:1/fs1:(1/fs1*(12000-1));
%xi=0:1/fs:(1/fs1*(12000-1));
%input_unfiltered2=interp1(X,input_unfiltered1,xi);

figure, h = plot(ff,20*log10(abs(h_LP)),...
                ff,20*log10(abs(h_HP)),...
                ff,20*log10(abs(h_LP+h_HP)));
%set(h,'LineWidth',2,{'Color'},{[0 0 1];[1 0 0];[0.7 0.9 0.2]});
axis([0 f_crossover*2 -50 5]);

title(['FIR filter Magnitude Response with N_{filt}=' num2str(N_filt)],'fontsize',16);
xlabel('Frequency [Hz]','fontsize',14), ylabel('Magnitude [dB]','fontsize',14);
grid on, legend('Low pass filter','High pass filter','Combined');
saveas(gcf,'Filter_MagResponse.tif');

%input_LP = filter(b_LP,1,input_unfiltered);
input_HP = filter(b_HP,1,input);
%NASABoom1_0Diff_1Alpha_rec1_complex_HP=filter(b_HP,1,NASABoom1_0Diff_1Alpha_rec1_complex);
FORTRAN_TEST_simple_HP=filter(b_HP,1,FORTRAN_TEST_simple);
% data3_HP=filter(b_HP,1,input_unfiltered3);
% data4_HP=filter(b_HP,1,input_unfiltered4);
print length(input_HP);
print length(FORTRAN_TEST_simple_HP);

corr1=xcorr(input_HP(1:18000),FORTRAN_TEST_simple_HP(:,2));
figure
plot(t,corr1)
timeshift_input_meas=-6.25*10^-4;
timeshift_input_rec1=.1467;
timeshift_sim_meas=timeshift_input_rec1+timeshift_input_meas;
FORTRAN_TEST_simple_HP = FORTRAN_TEST_simple_HP;
%^only for presentation, delete when done 

%T=N/fs;
figure, h = plot(tt-timeshift_sim_meas,FORTRAN_TEST_simple_HP(:,2));
figure, h = plot(tt-timeshift_sim_meas,NASABoom1_0Diff_1Alpha_rec1_simple_HP(:,2),
tt-timeshift_sim_meas,NASABoom1_0Diff_2Alpha_rec1_simple_HP(:,2),
tt-timeshift_sim_meas,NASABoom1_0Diff_3Alpha_rec1_simple_HP(:,2),
tt-timeshift_sim_meas,NASABoom1_0Diff_4Alpha_rec1_simple_HP(:,2),
tt,NASABoom1_EM17Measurement_HP(1:15000));

%set(h,'LineWidth',1,{'Color'},{[1 0 0];[1 1 0];[0 1 0];[0 0 1];[0 0 0]});
%axis([0 T 1.2*min(data_HP(:,2)) 1.2*max(data_HP(:,2))]);
title(('Simulated Vs. Measured Data for Various \alpha values at Receiver EM17'),'fontsize',20);
xlabel('Time [s]','fontsize',18), ylabel('Pressure [Pa]','fontsize',18);
grid on, legend('Simulated All Hard','Simulated Average','Simulated Complex','Simpulated All Soft','Measured Data');
%saveas(gcf,'/Volumes/Hermes/School/Results/EMBuildingSimpleGeometry/Figures/0Diff_Comp.tiff');

index=int16(abs(timeshift_input_rec1/T)+1);
index1=int16(abs(timeshift_input_meas/T)+1);
for j=1:NFFT

if (j<NFFT-index)
    FORTRAN_TEST_simple_shifted(j)=FORTRAN_TEST_simple_HP(j+index,2)-2*input_HP(j);
%    new4(j)=input_unfiltered4(j+index-1)-2*input(j);
elseif(j>=NFFT-index)
    FORTRAN_TEST_simple_shifted(j)=FORTRAN_TEST_simple_HP(j-NFFT+index+1,2);
%    new4(j)=2*input(j);
end
end

figure 
plot(tt,FORTRAN_TEST_simple_shifted);

corr_FORTRAN_TEST_simple=xcorr(FORTRAN_TEST_simple_shifted);


figure, h = plot(t,corr_FORTRAN_TEST_simple);
set(h,'LineWidth',1);
%axis([0 T 1.2*min(data_HP(:,2)) 1.2*max(data_HP(:,2))]);
title('corr-0Diff-1Alpha-rec1_simple','fontsize',16);
xlabel('Time [s]','fontsize',14), ylabel('Pressure [Pa]','fontsize',14);
grid on;
%saveas(gcf,'/Volumes/Hermes/School/Results/EMBuildingSimpleGeometry/Figure
%s/0Diff_Comp.tiff');


% figure, h = plot(tt,new,tt,new1,tt,new2,tt,new3,tt,new4);
% set(h,'LineWidth',1,{'Color'},{[1 0 0];[1 1 0];[0 1 0];[0 1 1];[0 0 1]});
% %axis([0 T 1.2*min(data_HP(:,2)) 1.2*max(data_HP(:,2))]);
% title(['Simulated Vs. Measured Data at Receiver EM17'],'fontsize',16);
% xlabel('Time [s]','fontsize',14), ylabel('Pressure [Pa]','fontsize',14);
% grid on, legend('Simulated 1 Alpha','Simulated 2 Alpha','Simulated 3 Alpha','Simulated 4 Alpha','Measured Data');
% %saveas(gcf,'/Volumes/Hermes/School/Results/EMBuildingSimpleGeometry/Figures/0Diff_Comp.tiff');
% %plot(tt,new,'r',tt,new1,'y',tt,new2,'g',tt,new3,'c',tt,new4(1:15000),'
% 
% corr=xcorr(new4,new);
% corr1=xcorr(new4,new1);
% figure
% plot(t,corr)
% figure
% plot(t,corr1)
%save(['input_LP_' num2str(f_crossover) '_' num2str(N_filt) '.txt'],'input_LP','-ascii');
%save(['input_HP_' num2str(f_crossover) '_' num2str(N_filt) '.txt'],'input_HP','-ascii');

% Fs=24000;
% f=Fs/2*linspace(0,1,length(data)/2);
% size=length(data(1044:18000))
% Y=fft(data(1044:18000,2),length(data(1044:18000)))/length(data(1044:18000));
% X=fft(data5(1:18000-1043),length(data(1044:18000)))/length(data(1044:18000));
% figure
% plot(data(1044:18000,1),data(1044:18000,2),'r',data(1044:18000,1),data5(1:18000-1043),'g')
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Data Vs. Simulation time domain')
% figure
% plot(data(1044:18000,1),((data(1044:18000,2).^2-data5(1:18000-1043).^2)))
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Comparison percent error time domain')
% figure
% plot(f(1:length(data(1044:18000))/2),((Y(1:length(data(1044:18000))/2)-X(1:length(data(1044:18000))/2))./X(1:length(data(1044:18000))/2))*100)
% xlabel('frequency')
% ylabel('pressure(Pa)')
% title('Comparison percent error frequency domain')
% figure
% plot(freqOctaves,(abs(octaves-octaves1))./octaves1*100)
% xlabel('frequency one third octaves')
% ylabel('pressure(Pa)')
% title('Comparison percent error octaves band frequency domain')
% time=data(1044:18000,1);
%  X=data(1044:18000,2);
%  Y=data5(1:18000-1043);
%  Z=data3(1044:18000,2);
%  A=data4(1044:18000,2);
% Xfreq=fft(X)/size;
% Yfreq=fft(Y)/size;
% Zfreq=fft(Z)/size;
% Afreq=fft(A)/size;
% PwSp=Yfreq.*conj(Yfreq);
% PwSp2=Yfreq.*conj(Zfreq);
% PwSp3=Yfreq.*conj(Afreq);
% figure
% plot(time,X,'g',time,Y,'r')
% figure
% plot(f(1:size/2),Xfreq(1:size/2),'g',f(1:size/2),Yfreq(1:size/2),'r')
% figure
% plot(f(1:size/2),abs(PwSp(1:size/2)),'g',f(1:size/2),abs(PwSp2(1:size/2)),'r',f(1:size/2),abs(PwSp3(1:size/2)),'black')
% figure
% plot(time,ifft(PwSp),'g',time,ifft(PwSp2),'r',time,ifft(PwSp3),'black')
% diff1=Y-X;
% diff2=Y-Z;
% diff3=Y-A;
% rms1=100*sqrt(sum(diff1.*diff1))/length(diff1)
% rms2=100*sqrt(sum(diff2.*diff2))/length(diff2)
% rms3=100*sqrt(sum(diff3.*diff3))/length(diff3)



% figure
% plot(data(:,1),data(:,2),'r')
% hold on
% plot(data1(:,1),data1(:,2),'b')
% plot(data2(:,1),data2(:,2),'g')
% plot(data3(:,1),data3(:,2),'c')
% plot(data4(:,1),data4(:,2),'m')
% plot(data1(:,1)+.0435,data5(:,1),'k')
% 
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Acoustically Hard Case')
% 
% file='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_0Diff_2Alpha_timeseries.dat';
% file1='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_10Diff_2Alpha_timeseries.dat';
% file2='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_20Diff_2Alpha_timeseries.dat';
% file3='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_50Diff_2Alpha_timeseries.dat';
% file4='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_100Diff_2Alpha_timeseries.dat';
% file5='/Volumes/Hermes/School/ResultsPart1/NASABoom1/NASABoom1_Reciever1.txt';
% [header,data]=hdrload(file);
% [header,data1]=hdrload(file1);
% [header,data2]=hdrload(file2);
% [header,data3]=hdrload(file3);
% [header,data4]=hdrload(file4);
% [header,data5]=hdrload(file5);
% 
% figure
% plot(data(:,1),data(:,2),'r')
% hold on
% plot(data1(:,1),data1(:,2),'b')
% plot(data2(:,1),data2(:,2),'g')
% plot(data3(:,1),data3(:,2),'c')
% plot(data4(:,1),data4(:,2),'m')
% plot(data1(:,1)+.0435,data5(:,1),'k')
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Average Absorption Case')
% 
% file='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_0Diff_3Alpha_timeseries.dat';
% file1='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_10Diff_3Alpha_timeseries.dat';
% file2='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_20Diff_3Alpha_timeseries.dat';
% file3='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_50Diff_3Alpha_timeseries.dat';
% file4='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_100Diff_3Alpha_timeseries.dat';
% file5='/Volumes/Hermes/School/ResultsPart1/NASABoom1/NASABoom1_Reciever1.txt';
% [header,data]=hdrload(file);
% [header,data1]=hdrload(file1);
% [header,data2]=hdrload(file2);
% [header,data3]=hdrload(file3);
% [header,data4]=hdrload(file4);
% [header,data5]=hdrload(file5);
% 
% figure
% plot(data(:,1),data(:,2),'r')
% hold on
% plot(data1(:,1),data1(:,2),'b')
% plot(data2(:,1),data2(:,2),'g')
% plot(data3(:,1),data3(:,2),'c')
% plot(data4(:,1),data4(:,2),'m')
% plot(data1(:,1)+.0435,data5(:,1),'k')
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Complicated Absorption Case')
% 
% file='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_0Diff_4Alpha_timeseries.dat';
% file1='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_10Diff_4Alpha_timeseries.dat';
% file2='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_20Diff_4Alpha_timeseries.dat';
% file3='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_50Diff_4Alpha_timeseries.dat';
% file4='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_100Diff_4Alpha_timeseries.dat';
% file5='/Volumes/Hermes/School/ResultsPart1/NASABoom1/NASABoom1_Reciever1.txt';
% [header,data]=hdrload(file);
% [header,data1]=hdrload(file1);
% [header,data2]=hdrload(file2);
% [header,data3]=hdrload(file3);
% [header,data4]=hdrload(file4);
% [header,data5]=hdrload(file5);
% 
% figure
% plot(data(:,1),data(:,2),'r')
% hold on
% plot(data1(:,1),data1(:,2),'b')
% plot(data2(:,1),data2(:,2),'g')
% plot(data3(:,1),data3(:,2),'c')
% plot(data4(:,1),data4(:,2),'m')
% plot(data1(:,1)+.0435,data5(:,1),'k')
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('Totally Absorptive Case')
% 
% file='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha1/Timeseries/Reciever1/NASABoom1_20Diff_1Alpha_timeseries.dat';
% file1='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha2/Timeseries/Reciever1/NASABoom1_20Diff_2Alpha_timeseries.dat';
% file2='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha3/Timeseries/Reciever1/NASABoom1_20Diff_3Alpha_timeseries.dat';
% file3='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_20Diff_4Alpha_timeseries.dat';
% %file4='/Volumes/Hermes/School/ResultsPart1/NASABoom1/Alpha4/Timeseries/Reciever1/NASABoom1_100Diff_4Alpha_timeseries.dat';
% file5='/Volumes/Hermes/School/ResultsPart1/NASABoom1/NASABoom1_Reciever1.txt';
% [header,data]=hdrload(file);
% [header,data1]=hdrload(file1);
% [header,data2]=hdrload(file2);
% [header,data3]=hdrload(file3);
% %[header,data4]=hdrload(file4);
% [header,data5]=hdrload(file5);
% 
% figure
% plot(data(:,1),data(:,2),'r')
% hold on
% plot(data1(:,1),data1(:,2),'b')
% plot(data2(:,1),data2(:,2),'g')
% plot(data3(:,1),data3(:,2),'c')
% %plot(data4(:,1),data4(:,2),'m')
% plot(data1(:,1)+.0435,data5(:,1),'k')
% xlabel('time(s)')
% ylabel('pressure(Pa)')
% title('20% Diffuse Case')
