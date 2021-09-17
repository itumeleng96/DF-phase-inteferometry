%Target Transmitter variables 
azimuth_arr_angle=30;
elevation_arr_angle=-30;
snr=0;

%----------------------------------------------------
%3-antenna array elements using Non-harmonic spacing
fc=6*10^9;                             %6GHz,12GHz,18GHz
lambda=(3*10^8)/(fc); 
lambda_spacing = (3*10^8)/(12*10^9);
num_elements=3;
d1=lambda_spacing/2;                    %Antenna Spacing 1
d2=(1)*lambda_spacing;                  %Antenna Spacing 2

%---------------------------------------------------
%Signal generation with with Pulse and SNR
%---------------------------------------------------
A=1;
w=2*pi*fc;
PW=0.1*10^-6 ;                  %Pulse width 0.05microseconds
PRI=0.2*10^-6 ;                 %Pulse repetition interval 0.2microseconds
fs = 40*10^9;                   %Sampling Frequency:40GHz
dt= 1/fs;
stopTime = 0.1*10^-6;           
t =(0:dt:stopTime);             %Not sure about the stop time
c= 3*10^8;

%-----------------------------------------------------------------
%Azimuth Angle calculations 
%-----------------------------------------------------------------
%ref antenna x(t)= Asin(wt+phi)                
x1=A*sin(w*t);                                   %The recieved signal on reciever 1            

%Antenna 2
phi2= 2*pi*d1*sind(azimuth_arr_angle)*(fc/c);      %The phase delay on receiver 2
x2=A*sin(w*t+phi2);                              %The recieved signal on reciever 2

%Antenna 3
phi3=2*pi*(d2)*sind(azimuth_arr_angle)*(fc/c);    %The phase delay on receiver 3
x3=A*sin(w*t+phi3);                              %The recieved signal on reciever 3

%--------------------------------------------------------------------------------------
%Elevation Angle calculations 
%---------------------------------------------------------------------------------------
%ref antenna x(t)= Asin(wt+phi)                
y1=A*sin(w*t);                                         %The recieved signal on reciever 1 -y             

%Antenna 2
y_phi2= 2*pi*d1*sind(elevation_arr_angle)*(fc/c);      %The phase delay on receiver 2 -y 
y2=A*sin(w*t+y_phi2);                                  %The recieved signal on reciever 2 - y

%Antenna 3
y_phi3=2*pi*(d2)*sind(elevation_arr_angle)*(fc/c);    %The phase delay on receiver 3 - y
y3=A*sin(w*t+y_phi3);                                 %The recieved signal on reciever 3 - y


rms_azimuth=zeros(1,10);
rms_elevation=zeros(1,10);
snr_array = zeros(1,10);

for j = 1:10
    aoa_azimuth=zeros(1,100);
    aoa_elevation = zeros(1,100);
    
    for i = 1:100
        %-------------------------------------------------------------------------------------
        %------ADD White Gaussian Noise to signals 
        %-----------------------------------------------------------------------------------
        x1_noise=awgn(x1,snr,'measured');
        x2_noise=awgn(x2,snr,'measured');
        x3_noise=awgn(x3,snr,'measured');

        y1_noise=awgn(y1,snr,'measured');
        y2_noise=awgn(y2,snr,'measured');
        y3_noise=awgn(y3,snr,'measured');

        %-----------------------------------------------------------------------------------
        %-------DOWN CONVERSION
        %-----------------------------------------------------------------------------------
        f_local=0.9*fc;                  %Local0scillator must be fraction of Fc:0.9*fc
        x_local=sin(2*pi*f_local*t);     %Locall oscillator signal

        x1_local= x1_noise.*x_local;
        x2_local= x2_noise.*x_local;
        x3_local= x3_noise.*x_local;


        y1_local= y1_noise.*x_local;
        y2_local= y2_noise.*x_local;
        y3_local= y3_noise.*x_local;

        %-----FFT and filtering ------------------
        X1_local = fft(x1_local);
        X2_local = fft(x2_local);
        X3_local = fft(x3_local);

        Y1_local = fft(y1_local);
        Y2_local = fft(y2_local);
        Y3_local = fft(y3_local);

        %-----Low pass Filter --------------------

        L=length(t);
        f2 = fs*(0:L-1)/L;

        rect =rectangularPulse(0,fc,f2);
        X1_IF =X1_local.*rect;
        X2_IF =X2_local.*rect;
        X3_IF =X3_local.*rect;

        Y1_IF =Y1_local.*rect;
        Y2_IF =Y2_local.*rect;
        Y3_IF =Y3_local.*rect;


        %-----------------------------------------------------------------------------------
        %      Frequency and Phase calculation
        %-----------------------------------------------------------------------------------

        % ----DEMODULATION 
        x1_demod = ifft(X1_IF).*x_local;
        x2_demod = ifft(X2_IF).*x_local;
        x3_demod = ifft(X3_IF).*x_local;

        y1_demod = ifft(Y1_IF).*x_local;
        y2_demod = ifft(Y2_IF).*x_local;
        y3_demod = ifft(Y3_IF).*x_local;

        %Low pass filter 
        X1_demod =fft(x1_demod);
        X2_demod =fft(x2_demod);
        X3_demod =fft(x3_demod);

        Y1_demod =fft(y1_demod);
        Y2_demod =fft(y2_demod);
        Y3_demod =fft(y3_demod);

        %Get the dominant frequency and get the phase at that point
        [M1,I1]=max(X1_demod);   
        [M2,I2]=max(X2_demod);
        [M3,I3]=max(X3_demod);

        [y_M1,y_I1]=max(Y1_demod);   
        [y_M2,y_I2]=max(Y2_demod);
        [y_M3,y_I3]=max(Y3_demod);

        phase_shift1=angle(M1);
        phase_shift2=angle(M2);
        phase_shift3=angle(M3);


        y_phase_shift1=angle(y_M1);
        y_phase_shift2=angle(y_M2);
        y_phase_shift3=angle(y_M3);

        aoa1=(asind(lambda*(((phase_shift2-phase_shift1)/(2*pi))+0)/d1));
        aoa2=[asind(lambda*(((phase_shift3-phase_shift1)/(2*pi))-1)/d2),asind(lambda*(((phase_shift3-phase_shift1)/(2*pi))+0)/d2),asind(lambda*(((phase_shift3-phase_shift1)/(2*pi))+1)/d2)];

        aoa1_y=(asind(lambda*(((y_phase_shift2-y_phase_shift1)/(2*pi))+0)/d1));    
        aoa2_y = [asind(lambda*(((y_phase_shift3-y_phase_shift1)/(2*pi))-1)/d2), asind(lambda*(((y_phase_shift3-y_phase_shift1)/(2*pi))+0)/d2), asind(lambda*(((y_phase_shift3-y_phase_shift1)/(2*pi))+1)/d2)];

        [val,idx]=min(abs(aoa2-aoa1));
        %angle_of_arrival_final_azimuth=aoa2(idx);
        [y_val,y_idx]=min(abs(aoa2_y-aoa1_y));
        %angle_of_arrival_final_elevation=aoa2_y(y_idx);
        aoa_azimuth(i)=aoa2(idx);
        aoa_elevation(i)=aoa2_y(y_idx);
    end
    rms_azimuth(j)= rms(aoa_azimuth-azimuth_arr_angle);
    rms_elevation(j) = rms(aoa_elevation-elevation_arr_angle);
    snr_array(j)=snr;
    snr=snr+2;
end
%Plot RMSE
figure();
plot(snr_array,rms_azimuth,'-s');
title('RMSE vs SNR ');
xlabel('SNR(dB)');
ylabel('RMSE');
legend('RMSE 30 degrees azimuth angle (6GHz)');

figure();
plot(snr_array,rms_elevation,'-s');
title('RMSE vs SNR ');
xlabel('SNR(dB)');
ylabel('RMSE');
legend('RMSE -30 degrees elevation angle (6GHz)');
