% --2022/11/10-- %
% --杨帆-- %
%% FMCW仿真
%% 可调整目标个数
%% 信号形式可选 *记得改成复信号形式进行下面的仿真
%% 根据论文《A Peer-to-Peer Interference Analysis for Automotive Chirp Sequence Radars》
%% 仿交叉干扰FMCW（非正弦） 调N_intchirp参数（一个干扰FMCW chirp的周期与发射chirp周期比值） 等于1时为平行干扰，产生假目标。 

clear all;
close all;
clc;

%% 目标参数设置
R_target = [45 50 80];
V_target = [10 6 14];
A_target = [1 0.1 0.5];      %不同目标对应的接收衰减
N_target = length(R_target); %目标个数

%% FMCW雷达仿真参数设置
c = 3e8;
fc = 77e9;                %chirp起始频率
B_chirp = 150e6;          %chirp带宽 来源师姐的程序数据
T_chirp = 26e-6;          %chirp持续时间
slope = B_chirp/T_chirp;  %chirp调频斜率
Fs = 16e6;                %采样率 
N_ADC = T_chirp *Fs;      %一个周期的采样点数 26*16=416点
N_chirp = 128;             %chirp周期数
% t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定采样点在Tchirp中的位置
t = 0:1/Fs:(N_ADC-1)*1/Fs;  %时间向量 确定采样点在Tchirp中的位置
lambda = c / fc ;         %波长

IF_store = zeros(N_chirp,N_ADC);

%% 估算FMCW雷达性能参数
distance_max = c*Fs/(4*slope);             %雷达探测最大距离 *注意论文中是/4，但是原理里面写的/2 
f_IF_max = slope * 2 * distance_max/c /2   %中频带宽
distance_res = c/(2*B_chirp);              %距离分辨率
velocity_max = lambda/(4*T_chirp);         %雷达最大测量速度
velocity_res = lambda/(2*N_chirp*T_chirp); %速度分辨率
f_d = 2 * fc * V_target / c;

%% 发射信号参数

AT = 1;  %发射信号增益

%% CW干扰参数设置
fint_CW = [77.2e9 77.3e9 77.5e9 77.55e9 77.6e9  78e9 78.5e9]; %CW波宽频带持续干扰正弦波 频率
R_int_CW = [1.5 20 10 7 8 45 50 ];                            %干扰的位置
V_int_CW = [1.5 8 2 10 12 4 5 ];                              %干扰速度
N_int_CW = length(fint_CW);                                     %CW干扰个数

%% FMCW干扰信号参数设置
R_int_FMCW = [80];
V_int_FMCW = [25];
N_int_FMCW = length(R_int_FMCW);

N_intchirp = 128;                             % 一个FMCW干扰跨越的发射chirp周期数
T_int_FMCW = N_intchirp*T_chirp;            % 干扰信号周期
N_int_division = ceil(N_chirp/N_intchirp);  %向上取整 干扰FMCW的周期数

fint_FMCW = fc;                         % FMCW干扰起始频率
B_int = B_chirp;                        % FMCW干扰带宽与发射chirp相同

slope_int = B_int/T_int_FMCW;        % FMCW干扰的调频斜率
N_FMCWint = floor(T_int_FMCW * Fs);  % 一个周期干扰信号点数 416* N_int_chirp
% t_int = 0:1/Fs:T_int-1/Fs;      % 干扰信号采样点位置
t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;   % 干扰信号采样点位置

%% 生成FMCW干扰信号
% S_int_FMCW = zeros(1,N_FMCWint);
% S_int_FMCW_targetsum = zeros(1,N_FMCWint);
S_int_FMCW_sum = zeros(1,N_FMCWint *N_int_division);

for int_FMCW = 0:1:N_int_division-1 %fmcw干扰周期数循环
    
    S_int_FMCW_targetsum = zeros(1,N_FMCWint); % !**注意**！定义的位置 
    
    %% FMCW干扰位置多普勒信息cell
    for n_fmcw = 1:1:N_int_FMCW
        distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw).*(t_int+int_FMCW *T_int_FMCW);  %CW干扰距离更新
        td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                       %干扰时延更新   *注意这里时延不是两倍了 干扰是直接到达接收端的
%         S_int_FMCW{n_fmcw,1} = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW干扰信号
        S_int_FMCW{n_fmcw,1} = AT*cos(2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2));      %FMCW干扰信号 *实信号
        S_int_FMCW_targetsum = S_int_FMCW_targetsum + S_int_FMCW{n_fmcw,1};
    end
    
%     S_int_FMCW = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW干扰信号
        S_int_FMCW_sum(int_FMCW*N_FMCWint+1:(int_FMCW+1)*N_FMCWint) = S_int_FMCW_targetsum;
end



%% 发射、接收数据生成 添加干扰信号
for i = 0:1:N_chirp-1  %chirp周期循环
    
    %% 目标位置多普勒信息cell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %距离更新
        t_d{k,1} = 2*distant{k,1}/c;                             %时延更新
    end
    
    %% CW干扰位置多普勒信息cell
    for n_cw = 1:1:N_int_CW
        distant_int_CW{n_cw,1} = R_int_CW(n_cw) +V_int_CW(n_cw).*(t+i*T_chirp);    %CW干扰距离更新
        td_int_CW{n_cw,1} = distant_int_CW{n_cw,1}/c;                              %干扰时延更新   *注意这里时延不是两倍了 干扰是直接到达接收端 没有经过目标反射
    end
    
%     %% FMCW干扰位置多普勒信息cell
%     for n_fmcw = 1:1:N_int_FMCW
%         distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw) .*(t+i*T_chirp); %FMCW干扰距离更新
%         td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                %FMCW干扰时延更新  *注意这里时延不是两倍了
%     end

  %% 发射信号
%     St = AT*exp(1i*2*pi*(fc*(t-i*T_chirp)+slope/2*t.^2)); %发射信号
    St = AT*cos(2*pi*(fc*t+slope/2*t.^2)); %发射信号
    
    Sr_sum = zeros(1,N_ADC);                              %一个chirp对应回波信号叠加空矩阵 1024点
    S_int_sum = zeros(1,N_ADC);                           %干扰信号空矩阵 1024点
    
    %% CW干扰信号
    for int_CW = 1:1:1 %干扰个数循环
        
        S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %正弦干扰信号(对应不同位置的干扰信号) 
        S_int_sum = S_int_sum + S_int_CW{int_CW,1};
    
    end
    
    
    %% 生成没有干扰的回波信号
    for m = 1:1:1
        
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %目标回波信号
        Sr{m,1} = A_target(m)*AT*cos(2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %目标回波信号
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2)) + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %加入FMCW干扰
        Sr_sum = Sr_sum + Sr{m,1};
        
    end
    
    %% ----添加噪声模块-----
    %% CW噪声
    if i <= 50
        
        Sr_sum_withCW = Sr_sum + S_int_sum;
        
    else
        Sr_sum_withCW = Sr_sum; 
    end
    
    %% FMCW噪声

    Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %加入干扰
    
    %% 混频得到中频
%     IF = St .*conj(Sr_sum);            %无干扰
%     IF = St .*conj(Sr_sum_withCW);     %CW噪声
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW噪声

    IF_I = St .* Sr_sum_withFMCW;   %FMCW噪声
    
    IF = IF_I + 1i*Sr_sum_withFMCW .*sin(2*pi*(fc*t+slope/2*t.^2));             %变成解析信号 90°相移器 除去负频率部分
    
    %% 设置低通滤波
    %% 设计巴特沃斯低通滤波器 模拟or数字？
%     omg_p = 8e6;
%     omg_s = 8.15e6;
%     w_p = 2*pi*omg_p;
%     w_s = 2*pi*omg_s;
%     Rp = 1;
%     As = 30
%     [N,wc] = buttord(w_p,w_s,Rp,As,'s'); %'s'代表模拟滤波器
%     [B,A] = butter(N,wc);
    
%     %% 滤波器幅频曲线
%     f_filter = 0:10:1.65e7;
%     w = 2*pi * f_filter;
%     Hk = freqs(B,A,w);
%     figure;
%     plot(f_filter/1e6,abs(Hk));
    
    %% 滤波
%     IF = filter(B,A,IF);
   

    IF_mat(i+1,:) = IF;              %用于时频分析
    IF_mat_withnoise(i+1,:) = IF;    %将带有噪声的中频信号保存
    IF_I_mat(i+1,:) = IF_I;          %中频I路信号
    
    St_mat(i+1,:) = St;
    Sr_mat(i+1,:) = Sr_sum;
    
end

N_fastfft = 512;
N_slowfft =128;


%% 画中频信号频谱图
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %取中频矩阵一行1024点
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure;
plot(f_IF,IF_realfft);
title('中频信号fft')
xlabel('frequency(Hz)');


%% 信号时频分析 短时傅里叶变换(还有些问题)
wlen = 416; %窗函数大小
hop = wlen/4;
nfft =416;
win = blackman(wlen,'periodic');
% [S_if,F_if,T_if] = spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
figure;
spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
title('受干扰中频信号短时傅里叶时频图');



%% 加汉宁窗函数
range_win = hanning(N_ADC);      %生成range窗
doppler_win = hanning(N_chirp);  %生成doppler窗

%% range fft
for i = 1:1:N_chirp
    
    temp = IF_mat_withnoise(i,:); %.* range_win'; %对每一行进行fft
    temp_fft = fft(temp,N_ADC)./N_ADC;

    IF_mat_withnoise(i,:) = temp_fft;           %存入IF_mat_withnoise
    
end

%% doppler fft
for j = 1:1:N_ADC
    
    temp = IF_mat_withnoise(:,j); %.* doppler_win; %再对每一列进行fft
    temp_fft = fftshift(fft(temp,N_chirp))./N_chirp;      
    IF_mat_withnoise(:,j) = temp_fft;            %存入IF_mat_withnoise
    
end

IF_mat_withnoise = abs(IF_mat_withnoise);
RDM = 10*log10(IF_mat_withnoise); %换算dB

Gs_I = 10*log10(N_ADC)%干扰增益比

%% 画RD图
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
