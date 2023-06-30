%% 超参数
% c = 3e8; %光速
% fc = 76.5e9;  %发射信号载频 中心频率
% bw = 500e6;  %发射信号带宽
% Tr = 10e-6;  %扫频时间 也就是周期
% N = 256;  %采样点
% Fs = 25.6e6; %采样率
% M = 256;  %chirp的数目
% k = bw/Tr; %chirp斜率
% index = 1:1:N;  %产生点向量
% IF_mat = zeros(M,N);  %存储带有噪声的中频信号

%% 参数设置
c = 3e8;
fc = 77e9; %chirp起始频率
bw = 1.6e9; %chirp带宽
Tr = 62.48*10^(-6) %chirp持续时间
N = 1024; %chirp周期内采样点数
M = 128; %chirp个数
Fs = N/Tr; %采样率 
k = bw/Tr; %chirp调频斜率
index = 1:1:N;
IF_mat = zeros(M,N);

%% 发射信号参数
AT = 1;  %发射信号增益
t = 0:1/Fs:Tr-1/Fs;  %时间向量 确定256个点在一个Tr中的每个时刻
%t = t - Tr/2; %将fc作为中心频率

%% 回波信号参数
distance = 80;  %目标距离雷达50m的距离
t_d = 2 * distance / c;  %目标距离雷达的延迟
velocity = 10;  %目标距雷达的相对速度为30m/s
f_d = 2 * (fc - bw/2) * velocity / c;  %多普勒频移
AR = 0.8;  %回波信号衰减的比例值

%% 生成数据
for i = 1:1:M  %chirp的循环
%     St = AT*exp((1i*2*pi)*(fc*(t-i*Tr)+k/2*t.^2)); %发射信号
    St = AT*cos(2*pi*(fc*(t-i*Tr)+k/2*t.^2));
%     Sr = AR*AT*exp((1i*2*pi)*((fc-f_d)*(t-t_d-i*Tr)+k/2*(t-t_d).^2));  %回波信号
    Sr = AR*AT*cos(2*pi*((fc-f_d)*(t-t_d-i*Tr)+k/2*(t-t_d).^2));  %回波信号
    %% 求回波信号的共轭
%     Sr = conj(Sr);  %求回波信号的共轭
    
    %% 求中频信号
    IF = St .* Sr + 1i*Sr.*sin(2*pi*(fc*(t-i*Tr)+k/2*t.^2));  %求中频信号
    SNR = 1e15;  %信噪比
    IF_with_Noise = awgn(IF,SNR,'measured');  %给中频信号加高斯白噪声，在添加噪声的时候，要进行能量的测量
    IF_mat(i,:) = IF_with_Noise;  %将带有噪声的中频信号保存
    IF_mat_1(i,:) = IF_with_Noise;
    IF_mat1(i,:) = IF;%用于频率分析 不用管
end
save('Ego_vehicle.mat', 'IF_mat');  %进行数据的保存

%% 加载数据
IF_mat = cell2mat(struct2cell(load('Ego_vehicle.mat','IF_mat')));
%% 超参数
% c = 3e8; %光速
% fc = 76.5e9;  %发射信号载频 中心频率
% bw = 500e6;  %发射信号带宽
% Tr = 10e-6;  %扫频时间 也就是周期
% N = 256;  %采样点
% Fs = 25.6e6; %采样率
% M = 256;  %chirp的数目
% k = bw/Tr; %chirp斜率

lambda = c / (fc - bw/2);
point = 1:1:N;  %产生点向量

%% 中频信号时频分析
tfr_IF =  tfrstft(IF_mat_1(1,1:1024)');
figure;
imagesc(abs(tfr_IF));
title('短时傅里叶时频图');

% %% 频率分析
% figure;
% IF_fft = fft(IF_mat1(1,:),1024);
% plot(IF_fft);

%% 生成窗
range_win = hamming(N);  %生成range窗
doppler_win = hamming(M);  %生成doppler窗
%% range fft
for i = 1:1:M
    temp = IF_mat(i,:) .* range_win';
    temp_fft = fft(temp,N);
    IF_mat(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N
    temp = IF_mat(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,M));
    IF_mat(:,j) = temp_fft;
end
%% 画图
figure;
distance_temp = (0:N - 1) * Fs * c / N / 2 / k;
speed_temp = (-M / 2:M / 2 - 1) * lambda / Tr / M / 2;
[X,Y] = meshgrid(distance_temp,speed_temp);
mesh(X,Y,(abs(IF_mat)));
xlabel('距离(m)');
ylabel('速度(m/s)');
zlabel('信号幅值');
title('2维FFT处理三维视图');

figure;
speed_temp = -speed_temp;
imagesc(distance_temp,speed_temp,abs(IF_mat));
title('距离-多普勒视图');
xlabel('距离(m)');
ylabel('速度(m/s)');