% --2022/10/25-- %
% --杨帆-- %
% FMCW仿真5.0.ver

%% 可调整目标个数
%% 雷达参数设置中 带宽1.6G、chirp周期62.48us是师姐的程序  √
%% 用实信号建模 再用正交基带采样（拟用希尔伯特变换实现）    √

clear all;
close all;
clc;

%% 目标参数设置
R_target = [45];
V_target = [3];
N_target = length(R_target);

%% FMCW雷达仿真参数设置
c = 3e8;
fc = 77e9;                %chirp起始频率
B_chirp = 1.6e9;          %chirp带宽 来源师姐的程序数据
T_chirp = 62.48*10^(-6);  %chirp持续时间
N_ADC = 1024;             %chirp周期内采样点数
N_chirp = 128;            %chirp个数
Fs = N_ADC/T_chirp;       %采样率 
% t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
t = 0:1/Fs:(N_ADC-1)/Fs;  %**改了
slope = B_chirp/T_chirp;  %chirp调频斜率
lambda = c / fc ;         %波长
IF_store = zeros(N_chirp,N_ADC);

%% 估算FMCW雷达性能参数
distance_max = c*Fs/(2*slope);             %雷达探测最大距离
distance_res = c/(2*B_chirp);              %距离分辨率
velocity_max = lambda/(4*T_chirp);         %雷达最大测量速度
velocity_res = lambda/(2*N_chirp*T_chirp); %速度分辨率

%% 发射信号参数
AT = 10;                  %发射信号增益
AR = 0.8;                 %回波信号衰减的比例值
Anoise = 10;              %干扰信号增益
fnoise = 78e9;            %干扰信号功率

%% 干扰参数设置
fint = [77.2e9 77.3e9 77.5e9 77.55e9 77.6e9  78e9 78.5e9]; %宽频带持续干扰正弦波 频率
R_int = [1.5 20 10 7 8 45 50 ];                            %干扰的位置
V_int = [1.5 8 2 10 12 4 5 ];                              %干扰速度
N_int = length(fint);                                      %干扰个数

%% 发射、接收数据生成
for i = 0:1:N_chirp-1  %chirp周期循环
    %% 目标位置多普勒信息cell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %距离更新
        t_d{k,1} = 2*distant{k,1}/c;                             %时延更新
    end
    
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2));   %发射信号
    St = zeros(1,N_ADC);
    Sr_sum = zeros(1,N_ADC);                 %回波信号叠加空矩阵
    St = AT*cos(2*pi*(fc*t+slope/2*t.^2));   %发射信号
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2));
    
    %% 干扰多普勒信息
    for n = 1:1:N_int
        distant_int{n,1} = R_int(n) +V_int(n).*(t+i*T_chirp);    %干扰距离更新
        td_int{n,1} = 2*distant_int{n,1}/c;                      %干扰时延更新
    end
    
    %% 生成干扰信号
    S_int_sum = zeros(1,N_ADC);                           %干扰信号空矩阵 1024点
    
    for int = 1:1:N_int %干扰个数循环
        
        S_int{int,1} = AT*exp(1i*2*pi*fint(int)*(t-i*T_chirp-td_int{int,1})); %正弦干扰信号 
%         S_int{int,1} = 2*AT*exp(1i*2*pi*(fc*(t-i*T_chirp-td_int{int,1})+slope/2*(t-td_int{int,1}).^2)); %正弦干扰信号
        S_int_sum = S_int_sum + S_int{int,1};
    
    end
   
    %% 生成回波信号
    for m = 1:1:N_target
%         Sr{m,1} = AR*AT*cos(2*pi*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2));
        Sr{m,1} = AR*AT*cos(2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));
%         Sr{m,1} = AR*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));
        Sr_sum = Sr_sum + Sr{m,1};
    end
    
%     Sr_itf = Anoise*exp(1i*2*pi*fnoise*(t+i*T_chirp));                   
%     Sr_itf_conj = conj(Snoise);                                           %正弦干扰信号共轭
    
    %% 混频得到中频
%     IF_I = St .* (Sr_sum + S_int_sum);                      %中频实信号
    IF_I = St .* Sr_sum ;                                   %无干扰中频实信号
%     IF = IF_I + 1i*IF_I;

    %% 低通滤波
%     IF_I = filter(Hd,IF_I);
    
    %% 变成复信号
    IF_Q = imag(hilbert(IF_I));               %希尔伯特变换得到中频信号的Q分量      
    IF = IF_I + 1i*IF_Q;
%     IF = St .* (Sr_conj+Snoise_conj);         %中频信号有干扰时
%     SNR = 10;                                 %信噪比
%     IF = awgn(IF,SNR,'measured');             %给中频信号加高斯白噪声，在添加噪声的时候，要进行能量的测量
    IF_mat(i+1,:) = IF_I;                       %将带有噪声的中频信号保存
    IF_I_mat(i+1,:) = IF_I;
    IF_Q_mat(i+1,:) = IF_Q;     
end

save('Ego_vehicle_复信号.mat', 'IF_mat');

%% 汉宁窗函数
range_win = hanning(N_ADC);      %生成range窗
doppler_win = hanning(N_chirp);  %生成doppler窗

%% range fft
for i = 1:1:N_chirp
    temp = IF_mat(i,:) .* range_win';
    temp_fft = fft(temp,N_ADC);
    IF_mat(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N_ADC
    temp = IF_mat(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,N_chirp));
    IF_mat(:,j) = temp_fft;
end


IF_mat = abs(IF_mat(:,1:N_ADC/2));

RDM = 10*log10(IF_mat);

%% 画图
% figure;
% mesh(IF_mat);
% xlabel('fast-time');
% ylabel('slow_time');
% zlabel('Amplitude');
% title('中频采样信号时域');
% 
% figure;
% plot(IF_mat(1,:));
% hold on
% plot(IF_Q_mat(1,:));
% xlabel('samples');
% ylabel('Amplitude');
% title('I-Q路信号');
% legend('I路','Q路');

figure;
distance_temp = (0:N_ADC/2 - 1) * Fs * c / N_ADC / 2 / slope;
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