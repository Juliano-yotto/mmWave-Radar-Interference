%RCS��Ȼ��һ�������⣬��Ӧ��Ϊ��ֵ
clear all;
close all;
clc;

%����ز��źŵĹ���
Pt_dBm = 12;    %�״�ķ��书�ʣ���λΪdBm
Pt_W = 0.01585; %�״�ķ��书�ʣ���λΪW
Gt_dB = 10.5;  %�״�ķ����������棬��λΪdBi��dB��
Gt = 11.22; %�״�ķ�����������
% Gr_dB = 30; %�״�������ߵ�����,��λΪdB
% Gr = 1000;  %�״�������ߵ�����
Gr_dB = 10.5; %�״�������ߵ�����,��λΪdB
Gr = 11.22;  %�״�������ߵ�����
R = 7;    %�״���Ŀ������ľ��룬��λΪm
fRadar = 77*10^9;  %�״�Ĺ���Ƶ��
c = 3*10^8; %����
lambda = c/fRadar;  %��������λΪm
RCS = 10;   %�״�ɢ������������λΪm^2
RCS_dBm = 10*log10(RCS);
Aeff = Gr*lambda^2/(4*pi);
% Pr = Pt_W*Gt*RCS*Aeff/(4*pi*R^2*4*pi*R^2); %�״���յ��Ļز��źŵĹ���
Pr = Pt_W*Gt*RCS*Aeff/(4*pi*R^2*4*pi*R^2); %�״���յ��Ļز��źŵĹ���
% Pr_dBm = 10*log10(Pr);
Pr_dBm = 10*log10(Pr*1000);
% Pr_dBm_1 = Pt_dBm+Gt_dB+Gr_dB+RCS_dBm+10*log10(lambda^2/((4*pi)^3*R^4));


%����ز��źŵ����
N = 353;    %������
fs = 3000*10^3;    %����Ƶ��
tADC=(N-1)/fs;  %��������ʱ��
tStart=6.4*10^(-6);  %������ʼʱ��
t=tStart:1/fs:tStart+tADC;


td=2*R/c;   %�ز��źŵ�ʱ��
tChirp=125.07*10^(-6);   %chirp�ĳ���ʱ��
k=25.49*10^12; %б��
B=k*tChirp; %����
fc=77*10^9; %�ز�Ƶ��
f=fc+k*t;

%�״﷢����ź�
st = cos(2*pi*fc*t + pi*k*t.^2);


%�״���յ��ź�
sr = cos(2*pi*fc*(t-td) + pi*k*(t-td).^2);
  plot(t,st);
  hold on;
  plot(t,sr);
  
RX_Gain_dB = 30;
RX_Gain = 10^(RX_Gain_dB / 10);
Pr_Amplified = Pr * RX_Gain;

syms A;
eqn = (A/sqrt(2))^2/50-Pr_Amplified == 0;
ACal=double(solve(eqn)); %����Pr_Amplified����ز��źŵ����
A_Amplified = ACal(2,1);


s_I_Amplified = A_Amplified*cos((2*pi*k*td)*t+2*pi*fc*td-pi*k*td^2); %�ز��ź�
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



