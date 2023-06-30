% 大论文
% 一个目标
% 人为正弦干扰影响分析

%% 最大探测距离约51.89m（c*fs/2s），距离分辨率约0.09m(c/2B)
%% 最大测量速度   无
%% 静止单目标
%% 139行 为什么要加窗？
% fs=200
% t=cell(2,1);
% t(1,1)={1:1/fs:2};
% t(2,1)={1:1/fs:3};

clear all;
close all;
clc;

c = 3*10^8; %光速
R1 = 1.3; %雷达与目标物体1的距离，单位为m


td_1 = 2*R1/c;   %目标物体1对应的回波信号时延

tChirp=62.48*10^(-6);   %chirp的持续时间
k=26.023*10^12; %斜率
B=k*tChirp; %带宽  约1.6GHz
fc=77*10^9; %载波频率 

N = 512;    %中频信号采样点
fs = 9000*10^3;    %采样频率9MHz
tADC=(N-1)/fs;  %采样持续时间 56.7us
t_start=4*10^(-6);  %采样开始时间
t=t_start:1/fs:(t_start+tADC);

% 计算回波信号的功率
Pt_dBm = 12;    %雷达的发射功率，单位为dBm
Pt_W = 0.01585; %雷达的发射功率，单位为W
% Pt_dBm = 18;    %雷达的发射功率，单位为dBm
% Pt_W = 0.06309573444801932; %雷达的发射功率，单位为W
Pt_dB = 10*log10(Pt_W);
Gt_dB = 10.5;  %雷达的发射天线增益，单位为dBi（dB）
Gt = 11.22; %雷达的发射天线增益
Gr_dB = 10.5; %雷达接收天线的增益,单位为dB
Gr = 11.22;  %雷达接收天线的增益

lambda = c/fc;  %波长，单位为m
RCS_1 = 0.2;   %目标物体1对应的雷达散射截面面积，单位为m^2


Pr_W_1 = Pt_W*Gt*Gr*lambda^2*RCS_1/((4*pi)^3*R1^4);


R_r = 50;
Ar_1 = sqrt(Pr_W_1 * R_r * 2);


Ar_1_Amplified = 31.62*Ar_1;


VRef = 1.8;
ADCCode_1 = Ar_1_Amplified*(2^16)/VRef; %目标物体1对应的ADC采样结果振幅


%人为正弦干扰信号参数
% PI_W = 5*10^(-7); %接收到的人为正弦干扰信号功率

% 根据相关信息计算接收到的干扰信号功率
P_itx = 10; %人为正弦干扰发射功率，单位为dBm
txAntGain = 20; %喇叭天线增益，单位为dB
rxAntGain = 10.5; %雷达接收天线增益，单位为dB
Ri = 1.3; %喇叭天线与雷达之间的距离
PI_Calculation_dBm = P_itx + txAntGain + rxAntGain - 10*log10((4*pi*Ri/lambda)^2); %单位为dBm
PI_W = 10^(PI_Calculation_dBm / 10)*10^(-3);


% PI_W = 3.16*10^(-7); %接收到的人为正弦干扰信号功率
A_I = sqrt(PI_W * R_r * 2);
A_I_Amplified = 31.62*A_I;
ADCCode_I = A_I_Amplified*(2^15)/VRef; %人为正弦干扰对应的中频信号采样结果振幅
fi = 78*10^9; %人为正弦干扰信号的频率
% fi = 78*10^9; %人为正弦干扰信号的频率
% rand返回一个在区间 (0,1) 内均匀分布的随机数
phi = 2*180*rand;


BIF = 15*10^6; %中频带宽为15MHz 低通滤波器带宽15MHz
tGlitch = BIF/k; %受到正弦干扰影响的持续时间
t_start_i = (fi - fc)/k; %开始受到正弦干扰影响的时间

%未受到干扰影响的I、Q两路中频信号的时域表达式
s_I = ADCCode_1*cos((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);
s_Q = ADCCode_1*sin((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);

% s_I = ADCCode_1*cos((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);
% s_Q = ADCCode_1*sin((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);

figure
plot(t  - t_start,s_I,'Color','blue','LineWidth',4);
xlabel('时间','FontSize',36,'FontWeight','bold');
ylabel('ADC值','FontSize',36,'FontWeight','bold');
ylim([-500 500])
hold on;
plot(t  - t_start,s_Q,'Color','red','LineWidth',4);
hold on;
legend('I路中频信号','Q路中频信号','FontSize',28,'FontWeight','bold')



%受到干扰影响的I、Q两路中频信号的时域表达式
t1 = t_start:1/fs:t_start_i;
t2 = t_start_i:1/fs:(t_start_i+tGlitch);
t3 = (t_start_i+tGlitch):1/fs:(t_start+tADC); %时间分为干扰前、干扰信号段、干扰后

s_I_1 = ADCCode_1*cos((2*pi*k*td_1)*t1+2*pi*fc*td_1-pi*k*td_1^2);
s_I_2 = ADCCode_1*cos((2*pi*k*td_1)*t2+2*pi*fc*td_1-pi*k*td_1^2) + ADCCode_I*sin(2*pi*(fi-fc)*t2-pi*k*t2.^2 + phi); %加了随机相位？
s_I_3 = ADCCode_1*cos((2*pi*k*td_1)*t3+2*pi*fc*td_1-pi*k*td_1^2);

s_Q_1 = ADCCode_1*sin((2*pi*k*td_1)*t1+2*pi*fc*td_1-pi*k*td_1^2);
s_Q_2 = ADCCode_1*sin((2*pi*k*td_1)*t2+2*pi*fc*td_1-pi*k*td_1^2) + ADCCode_I*cos(2*pi*(fi-fc)*t2-pi*k*t2.^2 + phi); %受到干扰的信号段
s_Q_3 = ADCCode_1*sin((2*pi*k*td_1)*t3+2*pi*fc*td_1-pi*k*td_1^2);

t_sum = [t1 t2 t3];
s_I_sum = [s_I_1 s_I_2 s_I_3];
s_Q_sum = [s_Q_1 s_Q_2 s_Q_3];


figure
plot(t_sum - t_start,s_I_sum,'Color','blue','LineWidth',4);
xlabel('时间','FontSize',36,'FontWeight','bold');
ylabel('ADC值','FontSize',36,'FontWeight','bold');
hold on
plot(t_sum - t_start,s_Q_sum,'Color','red','LineWidth',4);
hold on
legend('I路中频信号','Q路中频信号','FontSize',28,'FontWeight','bold')

din = s_I + 1j* s_Q; % 未受到干扰的中频IQ信号
din_i = s_I_sum + 1j* s_Q_sum; % 受到人为正弦干扰的中频采样信号

range_win = hanning(N); %加长为N=512的汉宁窗 
range_win = range_win';
wSum=sum(range_win(1,:));

temp = din.*range_win;  %无干扰信号加窗fft
temp_fft =  abs(fftshift(fft(temp, N))); 

temp_i = din_i.*range_win;  %%干扰信号加窗fft
temp_i_fft =  abs(fftshift(fft(temp_i, N)));

index = -N/2 + 0.5: 1: N/2 - 0.5;
freq_bin = (index - 1) * fs / N;
range_bin = freq_bin * c / 2 / k;

Pout_fft = 20*log10(temp_fft);
Pout_fft_i = 20*log10(temp_i_fft);

b_adc = 16;
cf = 20*log10(2^(b_adc-1))+20*log10(wSum)-20*log10(sqrt(2));

Pout_dBFs = Pout_fft - cf;
Pout_dBFs_i = Pout_fft_i - cf;

% Pout_threshold = 20*log10(XT);
% Pout_threshold_dBFs = Pout_threshold - cf;
% 
% 
% i1 = -i_length/2 + 0.5: 1: i_length/2 - 0.5;
% freq_bin_CFAR = (i1 - 1) * fs / N;
% range_bin_CFAR = freq_bin_CFAR * c / 2 / k;
% Pout_CFAR = 20*log10(XT) - cf;

% M=12;
% pro_N=3;
% Q = N-M-pro_N;
% 
% PAD=10^(-4);
% [ind, XT] = cfar_whyu(temp_fft, M, pro_N, PAD);
% Pout_XT = 20*log10(XT);
% Pout_XT_dBFs = Pout_XT-cf;
% Q_ind = -Q/2+0.5:1:Q/2-0.5;
% freq_Q=(Q_ind-1)*fs /Q;
% range_Q = freq_Q*c/2/k;

% %正常工作时的FFT结果
% figure
% plot(range_bin, Pout_dBFs,'b','LineWidth',4);
% % title('fft(after window)')
% xlabel('距离（米）','FontSize',36,'FontWeight','bold')
% ylabel('FFT结果（dBFS）','FontSize',36,'FontWeight','bold')
% axis([-30 30, -250 0])
% hold on

%受到干扰时的FFT结果
figure
plot(range_bin, Pout_dBFs_i,'b','LineWidth',4);
% title('fft(after window)')
xlabel('距离（米）','FontSize',36,'FontWeight','bold')
ylabel('FFT结果（dBFS）','FontSize',36,'FontWeight','bold')
axis([-30 30, -150 0])
hold on


% plot(range_Q,Pout_XT_dBFs,'r');
% hold on

% result = [range_bin.' Pout_dBFs' Pout_dBFs_i.'];
% result_s = xlswrite('D:\matlabprojects\DifferentTransmittedPower\18dBm-DifferentTransmittedPower-fft.xls', result);    % 将result写入到wind.xls文件中
