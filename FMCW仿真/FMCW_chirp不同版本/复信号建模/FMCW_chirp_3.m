% --2022/10/21-- %
% --杨帆-- %
% FMCW仿真3.0.ver
%% 可调整目标个数
%% 采用复信号形式
%% 添加正弦干扰 速度单元出现干扰条带 各速度单元上的检测造成影响

clear all;
close all;
clc;

%% 目标参数设置
R_target = [45];
V_target = [3.5];
A_target = [0.8 0.1 0.5]; %目标对应的接收增益
N_target = length(R_target);

%% FMCW雷达仿真参数设置
c = 3e8;
fc = 77e9;                %chirp起始频率
fc_noise = 77.6e9;        %正弦干扰信号频率
B_chirp = 1.6e9;          %chirp带宽 来源师姐的程序数据
T_chirp = 62.48e-6;       %chirp持续时间
N_ADC = 1024;             %chirp周期内采样点数
N_chirp = 128;            %chirp个数
Fs = N_ADC/T_chirp;       %采样率 
t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
slope = B_chirp/T_chirp;  %chirp调频斜率
lambda = c / fc ;         %波长
IF_store = zeros(N_chirp,N_ADC);

%% 估算FMCW雷达性能参数
distance_max = c*Fs/(2*slope);             %雷达探测最大距离
distance_res = c/(2*B_chirp);              %距离分辨率
velocity_max = lambda/(4*T_chirp);         %雷达最大测量速度
velocity_res = lambda/(2*N_chirp*T_chirp); %速度分辨率

%% 发射信号参数
AT = 1;                  %发射信号增益
t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
%AR = 0.8;                 %回波信号衰减的比例值
Anoise = 10;              %干扰信号增益
fnoise = 78e9;            %干扰信号功率
%t = t - Tr/2; %将fc作为中心频率

%% 发射、接收数据生成
for i = 0:1:N_chirp-1  %chirp周期循环
    %% 目标位置多普勒信息cell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %距离更新
        t_d{k,1} = 2*distant{k,1}/c;                             %时延更新
    end
    
    St = AT*exp((1i*2*pi)*(fc*(t-i*T_chirp)+slope/2*t.^2));    %发射信号
    Sn = AT*exp(1i*2*pi*fc_noise*(t-i*T_chirp));               %干扰信号
    slope_N = 1.09e13;
    SN = AT*exp((1i*2*pi)*(fc_noise*(t-i*T_chirp)+slope_N/2*t.^2)); 
    Snconj = conj(Sn);                                         %干扰信号共轭
    Sr_sum = zeros(1,N_ADC);                                   %回波信号叠加空矩阵
    
    %% 生成回波信号
    for m = 1:1:N_target
        Sr{m,1} = A_target(m)*AT*exp((1i*2*pi)*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2));
        Srconj{m,1} = conj(Sr{m,1});
        Sr_sum = Sr_sum + Srconj{m,1};
    end
    
    %% 混频得到中频
    IF = St .* Sr_sum ;
%     IF = St .* (Sr_sum+Sn);    %中频信号有干扰时
%     SNR = 10;                      %信噪比
%     IF = awgn(IF,SNR,'measured');  %给中频信号加高斯白噪声，在添加噪声的时候，要进行能量的测量
    IF_mat_1(i+1,:) = IF;            %将带有噪声的中频信号保存
end
save('Ego_vehicle_复信号.mat', 'IF_mat_1');

%% 汉宁窗函数
range_win = hanning(N_ADC);      %生成range窗
doppler_win = hanning(N_chirp);  %生成doppler窗

%% range fft
for i = 1:1:N_chirp
    temp = IF_mat_1(i,:) .* range_win';
    temp_fft = fft(temp,N_ADC)./N_ADC;
    IF_mat_1(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N_ADC
    temp = IF_mat_1(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,N_chirp));
    IF_mat_1(:,j) = temp_fft;
end

IF_mat_1 = abs(IF_mat_1);
RDM = 10*log10(IF_mat_1);

%% 画图
figure;
distance_temp = (0:N_ADC - 1) * Fs * c / N_ADC / 2 / slope;
speed_temp = (-N_chirp / 2:N_chirp / 2 - 1) * lambda / T_chirp / N_chirp / 2;
[X,Y] = meshgrid(distance_temp,speed_temp);
mesh(X,Y,RDM);
xlabel('distance(m)');
ylabel('velocity(m/s)');
zlabel('Amplitude');
title('range-doppler-3D fft');

figure;
imagesc(distance_temp,speed_temp,RDM);
title('range-doppler');
xlabel('distance(m)');
ylabel('velocity(m/s)');
