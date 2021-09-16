
% Target Transmitter variables 
azimuth_arr_angle=10;
elevation_arr_angle=0;
snr=-10;
%----------------------------------------------------
%3-antenna array elements using Non-harmonic spacing
lambda=(3*10^8)/(6*10^9);       %6GHz
num_elements=3;
d1=lambda;                      %Antenna Spacing 1
d2=(3/2)*lambda;                %Antenna Spacing 2

%---------------------------------------------------
%Signal generation with with Pulse and SNR
%---------------------------------------------------
A=1;
fc=6*10^9;
w=2*pi*fc;
PW=0.1*10^-6 ;                  %Pulse width 0.05microseconds
PRI=0.2*10^-6 ;                 %Pulse repetition interval 0.2microseconds
fs = 50*10^9;                   %Sampling Frequency:40GHz
dt= 1/fs;
stopTime = 1*10^-6;             %seconds 1microSeconds
t =(0:dt:stopTime);             %Not sure about the stop time
c= 3*10^8;
%-----------------------------------------------------------------
%Azimuth Angle calculations 
%-----------------------------------------------------------------
%ref antenna x(t)= Asin(wt+phi)                
x1=A*sin(w*t);                                   %The recieved signal on reciever 1            

%Antenna 2
phi1=2*pi*d1*sin(azimuth_arr_angle)*(fc/c);      %The phase delay on receiver 2
x2=A*sin(w*t+phi1);                              %The recieved signal on reciever 2

%Antenna 3
phi2=2*pi*(d2)*sin(azimuth_arr_angle)*(fc/c);    %The phase delay on receiver 3
x3=A*sin(w*t+phi2);                                   %The recieved signal on reciever 3

%-------------------------------------------------------------------------------------
%------ADD White Gaussian Noise to signals 
%-----------------------------------------------------------------------------------
x1_noise=awgn(x1,snr,'measured');
x2_noise=awgn(x2,snr,'measured');
x3_noise=awgn(x3,snr,'measured');

%-----------------------------------------------------------------------------------
%-------DOWN CONVERSION
%-----------------------------------------------------------------------------------
f_local=0.9*fc;                  %Local0scillator must be fraction of Fc:0.9*fc
x_local=sin(2*pi*f_local*t);     %Locall oscillator signal

x1_local= x1_noise.*x_local;
x2_local= x2_noise.*x_local;
x3_local= x3_noise.*x_local;

%-----IFT and filtering ------------------
X1_local = fft(x1_local);
X2_local = fft(x2_local);
X3_local = fft(x3_local);

%-----Low pass Filter --------------------

L=length(t);
f2 = fs*(0:L-1)/L;

rect =rectangularPulse(0,fc,f2);
x1_IF =X1_local.*rect;
x2_IF =X2_local.*rect;
x3_IF =X3_local.*rect;

%P2 = abs(x3_IF/L);
%P1 = P2(1:L/2+1);
%P1(2:end-1) = 2*P1(2:end-1);
%f = fs*(0:(L/2))/L;
%plot(f,P1) 
%title('Single-Sided Amplitude Spectrum of X(t)')
%xlabel('f (Hz)')
%ylabel('|P1(f)|')

%-----------------------------------------------------------------------------------
%      Frequency and Phase calculation
%-----------------------------------------------------------------------------------

[M1,I1]=max(x1_IF);   %Get the dominant frequency and get the phase at that point
[M2,I2]=max(x2_IF);
[M3,I3]=max(x3_IF);

phase_shift1=angle(M1);
phase_shift2=angle(M2);
phase_shift3=angle(M3);

aoa1=angle(asind((lambda*(phase_shift2-phase_shift1)/2*pi+1)/d1));
aoa2=angle(asind((lambda*(phase_shift3-phase_shift2)/2*pi+2)/d2));
%figure()
%plot(t,x1_local);
%second antenna 
%phi=
%x2=sin(w*t+phi)

%third antenna

