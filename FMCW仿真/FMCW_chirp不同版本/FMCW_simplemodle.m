%--2022/10/2--%
%--杨帆--%
%本程序对简单FMCW波雷达测速、测距进行仿真。不考虑雷达信号发射间隔endle_time 保证得到中频为单频信号进行处理。为后续进行抗干扰仿真铺垫

%% 问题：
%% 1.matlab写程序时，要想完整表达发射信号与接收信号，那么时间轴应该加以区分  2022/10/6 已解决√
%% 2.时间延迟是否要算成变量？此程序写成变量td(i)，但是师姐的程序里面是定值td.TI文档中也近似忽略IF信号频率与物体速度依赖，看成定值 2022/10/6 已解决√
%% 3.试一下main.m里面的数据 采样率也要换！！！采样点数与采样率匹配问题
clear all;
close all;
clc;

%% 雷达参数设置
% maxR = 200;           % 雷达最大探测目标的距离
% rangeRes = 1;         % 雷达的距离分率
% maxV = 70;            % 雷达最大检测目标的速度
fc= 77e9;             % 雷达工作频率 载频
c = 3e8;              % 光速

%% 目标参数设置
r0 = 150; % 目标距离设置 (max = 200m)
v0 = 50; % 目标速度设置 (min =-70m/s, max=70m/s)


%% 雷达发射FMCW参数设置
T_chirp=7.33e-6;       %chirp周期
B_chirp=150e6;         %发射信号调频带宽
slope=B_chirp/T_chirp; %线性调频斜率68MHz/us


%  slope = 29.306e12;
%  T_chirp = 29.56e-6; %chirp周期
%  B_chirp = 750e6;    %发射信号调频带宽
%  fs=20e6 
N_chirp=128;                          %chirp数量 
N_ADC=1024;                           %ADC采样点数
fs=N_ADC/T_chirp;                     %中频模拟信号采样率 

t=linspace(0,N_chirp*T_chirp,N_ADC*N_chirp); %发射信号和接收信号的采样时间，在MATLAB中的模拟信号是通过数字信号无限采样生成
Tx=zeros(1,length(t));                       %发射信号
Rx=zeros(1,length(t));                       %接收信号
Mix = zeros(1,length(t));                    %差频、差拍、拍频、中频信号

r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% 动目标信号生成
fint = fc;
Tint = N_chirp*T_chirp;
Bint = B_chirp;
slope_int = Bint/Tint;

for i=1:length(t) 
    
    r_t(i) = r0 + v0*t(i); % 距离更新
    td(i) = 2*r_t(i)/c;    % 延迟时间
    
    Tx(i) = cos(2*pi*(fc*t(i) + (slope*t(i)^2)/2));                 %发射信号 实数信号
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2/2))); %接收信号 实数信号
    int(i) = cos(2*pi*(fint*t(i) + slope_int*t(i)^2/2));            %FMCW干扰信号
    Rxint(i) = Rx(i) + int(i);                                      %加干扰
    
    if i<=1024
         freq(i)=fc+slope*i;       %发射信号时频图 只取第一个chirp
         freq_echo(i)=fc+slope*i;  %回波信号频谱延迟
    end

    Mix(i) = Tx(i).*Rx(i);        %差频、差拍、拍频、中频信号
    Mixint(i) = Tx(i).* Rxint(i); %含干扰信号的中频差拍信号
end
Tx;
Rx;
Rxint;

%% 时频分析 存在交叉项
sigint = reshape(Mixint,N_ADC,N_chirp);
wlen = 64; %窗函数大小
hop = wlen/8;
nfft =1024;
win = blackman(wlen,'periodic');
% [S_if,F_if,T_if] = spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
figure;
spectrogram(sigint(:,128), win, wlen - hop, nfft, fs);
title('受干扰中频信号短时傅里叶时频图');


%发射信号时域图
figure;
plot(Tx(1:1024));
xlabel('点数');
ylabel('幅度');
title('TX发射信号时域图');
hold on;
%接收信号时域图
plot(Rx(1:1024));
xlabel('点数');
ylabel('幅度');
title('RX接收信号时域图');

% %发射信号时频图
% figure;
% plot(t(1:1024),freq);
% xlabel('时间');
% ylabel('频率');
% title('TX发射信号时频图');


% %接收信号时域图
% figure;
% plot(Rx(1:1024));
% xlabel('点数');
% ylabel('幅度');
% title('RX接收信号时域图');

%接收信号与发射信号的时频图
% figure;
% plot(t(1:1024),freq);
% hold on;
% plot(td(1:1024)+t(1:1024),freq);
% xlabel('时间');
% ylabel('频率');
% title('接收信号与发射信号时频图');
% legend ('TX','RX');

%中频信号频谱 和频信号观察
%figure;
% plot(db(abs(fft(Mix(1:1024*256)))));%查看宽带的和频信号 将chirp的点数改为1024*256即可看到有一个门信号，但注意计算机内存。
% xlabel('频率');
% ylabel('幅度');
% title('中频信号频谱');

% figure;
% plot(db(abs(fft(Mix(1:1024)))));%查看宽带的和频信号 将chirp的点数改为1024*256即可看到有一个门信号，但注意计算机内存。只针对第一个chirp对应的中频1024点
% xlabel('频率');
% ylabel('幅度');
% title('中频信号频谱');

%% 低通滤波 截止频率30MHz  采样频率120MHz
% Mix=lowpass(Mix(1:1024*256),30e6,120e6);
% plot(db(abs(fft(Mix(1:1024*256)))));
% xlabel('频率');
% ylabel('幅度');
% title('中频信号低通滤波器');

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
signal = reshape(Mix,N_ADC,N_chirp); %*对混频信号矩阵进行重构，变成1024×128的矩阵，方便后续进行距离-速度维FFT

figure;
mesh(signal);
xlabel('脉冲数')
ylabel('距离门数');
title('中频信号时域');

%% 距离维FFT
sig_fft = fft(signal,N_ADC)./N_ADC;
sig_fft = abs(sig_fft);
sig_fft = sig_fft(1:(N_ADC/2),:);

figure;
plot(sig_fft(:,1));
xlabel('距离（频率）');
ylabel('幅度')
title('单个chirp的FTF结果')


%% 距离FFT结果谱矩阵
figure;
mesh(sig_fft);
xlabel('距离（频率）');
ylabel('chirp脉冲数')
zlabel('幅度')
title('距离维FTF结果')

%% 速度维FFT

Mix=reshape(Mix,[N_ADC,N_chirp]); %*MIX信号矩阵进行重构
sig_fft2 = fft2(Mix,N_ADC,N_chirp);

sig_fft2 = sig_fft2(1:N_ADC/2,1:N_chirp);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
%RDM = 10*log10(RDM) ;
doppler_axis = linspace(-100,100,N_chirp);
range_axis = linspace(-200,200,N_ADC/2)*((N_ADC/2)/400);

figure;
mesh(doppler_axis,range_axis,RDM);
xlabel('多普勒通道'); ylabel('距离通道'); zlabel('幅度（dB）');
title('速度维FFT 距离多普勒谱');