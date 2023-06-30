% --2022/10/21-- %
% --杨帆-- %
% FMCW仿真3.0.ver

%% 此程序主要讨论idel_time对速度测量的影响
%% 可调整目标个数
%% 采用复信号形式
%% 加入空闲时间idel_time

clear all;
close all;
clc;

%% 目标参数设置
R_target = [1.3 7 60];
V_target = [0.5 3 10];
N_target = length(R_target);

%% FMCW雷达仿真参数设置
c = 3e8;
fc = 77e9;                %chirp起始频率
B_chirp = 300e6%1.6e9;          %chirp带宽 来源师姐的程序数据
T_chirp = 30e-6%62.48*10^(-6);  %chirp持续时间
%N_ADC = 512%1024;             %chirp周期内采样点数
N_chirp = 256%128;            %chirp个数
Fs = 16.67e6;%N_ADC/T_chirp;       %采样率 
N_ADC = floor(T_chirp * Fs);
t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
slope = B_chirp/T_chirp;  %chirp调频斜率
lambda = c / fc ;         %波长
IF_store = zeros(N_chirp,N_ADC);

%% IDEL TIME设置(用于时域信号建模)
idel_time = 4e-6;                   %idel_time持续时间
N_idel = floor(idel_time * Fs);     %补零点数
t_zeros = 0:1/Fs:idel_time-1/Fs;    %确定采样零点在idel time中的时间位置
idel_zeros = 0 * t_zeros;           %补零矩阵

%% 估算FMCW雷达性能参数
distance_max = c*Fs/(2*slope);             %雷达探测最大距离
distance_res = c/(2*B_chirp);              %距离分辨率
velocity_max = lambda/(4*T_chirp);         %雷达最大测量速度
velocity_res = lambda/(2*N_chirp*T_chirp); %速度分辨率

%% 发射信号参数
AT = 10;                  %发射信号增益
t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
AR = 0.8;                 %回波信号衰减的比例值
Anoise = 10;              %干扰信号增益
fnoise = 78e9;            %干扰信号功率
%t = t - Tr/2; %将fc作为中心频率

%% 发射、接收数据生成
for i = 0:1:N_chirp-1  %chirp周期循环
    %% 目标位置多普勒信息cell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp+i*idel_time); %距离更新
        t_d{k,1} = 2*distant{k,1}/c;                                         %时延更新
    end
    
    St = AT*exp((1i*2*pi)*(fc*(t-i*T_chirp-i*idel_time)+slope/2*t.^2));   %发射信号
    St_full = [St,idel_zeros];                                            %空闲时间内补零
    Sr_sum = zeros(1,N_ADC);                                              %回波信号叠加空矩阵
    Srfull_sum = zeros(1,N_ADC + N_idel);
    
    %% 生成回波信号
    for m = 1:1:N_target
        Sr{m,1} = AR*AT*exp((1i*2*pi)*(fc*(t-t_d{m,1}-i*T_chirp-i*idel_time)+slope/2*(t-t_d{m,1}).^2));
        Srfull{m,1} = [Sr{m,1},idel_zeros];
        Srfull_conj{m,1} = conj(Srfull{m,1});
        Srconj{m,1} = conj(Sr{m,1});
        Sr_sum = Sr_sum + Srconj{m,1};
        Srfull_sum = Srfull_sum +Srfull_conj{m,1};
    end
%     Sr_itf = Anoise*exp(1i*2*pi*fnoise*(t+i*T_chirp));                   
%     Sr_itf_conj = conj(Snoise);                                           %正弦干扰信号共轭
    
    %% 混频得到中频
    IF = St .* Sr_sum;
    IF_full = St_full .* Srfull_sum;
%     IF = St .* (Sr_conj+Snoise_conj);         %中频信号有干扰时
    SNR = 10;                      %信噪比
    IF = awgn(IF,SNR,'measured');  %给中频信号加高斯白噪声，在添加噪声的时候，要进行能量的测量
    IF_mat(i+1,:) = IF;          %将带有噪声的中频信号保存
    IFfull_mat(i+1,:) = IF_full;
end


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

% %% 汉宁窗函数
% range_win = hanning(N_ADC);      %生成range窗1124点
% doppler_win = hanning(N_chirp);  %生成doppler窗
% 
% IFfull_mat1 = zeros(N_chirp,512); %用于存放二维fft值

% %% range fft
% for i = 1:1:N_chirp
%     temp = IF_mat(i,:) .* range_win';
%     temp_fft = fft(temp,512); %还是进行1024点变换
%     IFfull_mat1(i,:) = temp_fft;
% end
% %% doppler fft
% for j = 1:1:512
%     temp = IFfull_mat1(:,j) .* doppler_win;
%     temp_fft = fftshift(fft(temp,N_chirp));
%     IFfull_mat1(:,j) = temp_fft;
% end

%% 画图
figure;
distance_temp = (0:N_ADC - 1) * Fs * c / N_ADC / 2 / slope;
speed_temp = (-N_chirp / 2:N_chirp / 2 - 1) * lambda / T_chirp / N_chirp / 2;
[X,Y] = meshgrid(distance_temp,speed_temp);
mesh(X,Y,(abs(IF_mat)));
xlabel('distance(m)');
ylabel('velocity(m/s)');
zlabel('Amplitude');
title('range-doppler-3D fft');

figure;
imagesc(distance_temp,speed_temp,abs(IF_mat));
title('range-doppler');
xlabel('distance(m)');
ylabel('velocity(m/s)');