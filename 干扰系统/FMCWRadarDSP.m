%RCS仍然有一定的问题，不应该为定值
clear all;
close all;
clc;

%计算回波信号的功率
Pt_dBm = 12;    %雷达的发射功率，单位为dBm
Pt_W = 0.01585; %雷达的发射功率，单位为W
Gt_dB = 10.5;  %雷达的发射天线增益，单位为dBi（dB）
Gt = 11.22; %雷达的发射天线增益
% Gr_dB = 30; %雷达接收天线的增益,单位为dB
% Gr = 1000;  %雷达接收天线的增益
Gr_dB = 10.5; %雷达接收天线的增益,单位为dB
Gr = 11.22;  %雷达接收天线的增益
R = 7;    %雷达与目标物体的距离，单位为m
fRadar = 77*10^9;  %雷达的工作频率
c = 3*10^8; %光速
lambda = c/fRadar;  %波长，单位为m
RCS = 10;   %雷达散射截面面积，单位为m^2
RCS_dBm = 10*log10(RCS);
Aeff = Gr*lambda^2/(4*pi);
% Pr = Pt_W*Gt*RCS*Aeff/(4*pi*R^2*4*pi*R^2); %雷达接收到的回波信号的功率
Pr = Pt_W*Gt*RCS*Aeff/(4*pi*R^2*4*pi*R^2); %雷达接收到的回波信号的功率
% Pr_dBm = 10*log10(Pr);
Pr_dBm = 10*log10(Pr*1000);
% Pr_dBm_1 = Pt_dBm+Gt_dB+Gr_dB+RCS_dBm+10*log10(lambda^2/((4*pi)^3*R^4));


%计算回波信号的振幅
N = 353;    %采样点
fs = 3000*10^3;    %采样频率
tADC=(N-1)/fs;  %采样持续时间
tStart=6.4*10^(-6);  %采样开始时间
t=tStart:1/fs:tStart+tADC;


td=2*R/c;   %回波信号的时延
tChirp=125.07*10^(-6);   %chirp的持续时间
k=25.49*10^12; %斜率
B=k*tChirp; %带宽
fc=77*10^9; %载波频率
f=fc+k*t;

%雷达发射的信号
st = cos(2*pi*fc*t + pi*k*t.^2);


%雷达接收的信号
sr = cos(2*pi*fc*(t-td) + pi*k*(t-td).^2);
  plot(t,st);
  hold on;
  plot(t,sr);
  
RX_Gain_dB = 30;
RX_Gain = 10^(RX_Gain_dB / 10);
Pr_Amplified = Pr * RX_Gain;

syms A;
eqn = (A/sqrt(2))^2/50-Pr_Amplified == 0;
ACal=double(solve(eqn)); %根据Pr_Amplified计算回波信号的振幅
A_Amplified = ACal(2,1);


s_I_Amplified = A_Amplified*cos((2*pi*k*td)*t+2*pi*fc*td-pi*k*td^2); %回波信号
s_Q_Amplified = A_Amplified*sin((2*pi*k*td)*t+2*pi*fc*td-pi*k*td^2);



VIN = 0.5;
s_I_Amplified = s_I_Amplified*(2^15)/VIN;
s_Q_Amplified = s_Q_Amplified*(2^15)/VIN;



time=0:1/fs:tADC;
figure
plot(time,s_I_Amplified,'b');
hold on;
plot(time,s_Q_Amplified,'r');
hold on;


din = s_I_Amplified + 1j* s_Q_Amplified; % complex value
figure
plot(abs(din))
title('complex value')



range_win = hanning(N);
range_win = range_win';
wSum=sum(range_win(1,:));
temp = din.*range_win;
temp_fft =  abs(fftshift(fft(temp, N)));

index = -N/2 + 0.5: 1: N/2 - 0.5;
freq_bin = (index - 1) * fs / N;
range_bin = freq_bin * c / 2 / k;





Pout_fft = 20*log10(temp_fft);
b_adc = 16;
cf = 20*log10(2^(b_adc-1))+20*log10(wSum)-20*log10(sqrt(2));
Pout_dBFs = Pout_fft-cf;

figure
plot(range_bin, Pout_dBFs)
% yticks([-140 -120 -100 -80 -60 -40 -20 0])
title('fft(after window)')
xlabel('range(m)')



