%% 2022/11/18
%% 本程序用于抗交叉干扰算法研究
%% 干扰1：宽频长时间正弦信号 算干扰开始时间与结束时间 在每个chirp加入干扰
%% 针对长时间且频率不单一的正弦干扰信号（跨越的带宽很大） 程序中雷达带宽150M 
%% 设置正弦干扰[77.04G, 77.08G, 77.12G]跨越带宽0.08G=80M

%% 收发都是实信号建模


clear all;
close all;
clc;

%% 目标参数设置
R_target = [40 50 80];
V_target = [1.5 6 14];
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
N_chirp = 75;             %chirp周期数
% t = 0:1/Fs:T_chirp-1/Fs;  %时间向量 确定采样点在Tchirp中的位置
t = 0:1/Fs:(N_ADC-1)*1/Fs;%时间向量 确定采样点在Tchirp中的位置
lambda = c / fc ;         %波长

BAAF = 8e6;               %中频滤波器带宽8M
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
fint_CW = [ 77.040e9 77.08e9]; %CW波宽频带持续干扰正弦波 频率 
R_int_CW = [ 0 0 0];                            %干扰的位置
V_int_CW = [ 0 0 0];                              %干扰速度
N_int_CW = length(fint_CW);                                   %CW干扰个数


t_int_start = (fint_CW - fc - BAAF)/slope;   %干扰开始影响时间                     
n_int_start = floor(t_int_start*Fs);       %开始对应的点数
t_int_end = (fint_CW - fc + BAAF)/slope;     %干扰结束影响时间
n_int_end = floor(t_int_end*Fs);           %结束对应的点数
t_int_consist = 2*BAAF/slope;              %干扰影响总时长
n_int_consist = floor(t_int_consist*Fs);       %干扰影响总点数


% %% FMCW干扰信号参数设置
% R_int_FMCW = [80];
% V_int_FMCW = [25];
% N_int_FMCW = length(R_int_FMCW);
% 
% N_intchirp = 5;                             % 一个FMCW干扰跨越的发射chirp周期数
% T_int_FMCW = N_intchirp*T_chirp;            % 干扰信号周期
% N_int_division = ceil(N_chirp/N_intchirp);  %向上取整 干扰FMCW的周期数
% 
% fint_FMCW = fc;                         % FMCW干扰起始频率
% B_int = B_chirp;                        % FMCW干扰带宽与发射chirp相同
% 
% slope_int = B_int/T_int_FMCW;        % FMCW干扰的调频斜率
% N_FMCWint = floor(T_int_FMCW * Fs);  % 一个周期干扰信号点数 416* N_int_chirp
% % t_int = 0:1/Fs:T_int-1/Fs;      % 干扰信号采样点位置
% t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;   % 干扰信号采样点位置
% 
% %% 生成FMCW干扰信号
% % S_int_FMCW = zeros(1,N_FMCWint);
% % S_int_FMCW_targetsum = zeros(1,N_FMCWint);
% S_int_FMCW_sum = zeros(1,N_FMCWint *N_int_division);

% for int_FMCW = 0:1:N_int_division-1 %fmcw干扰周期数循环
%     
%     S_int_FMCW_targetsum = zeros(1,N_FMCWint); % !**注意**！定义的位置 
%     
%     %% FMCW干扰位置多普勒信息cell
%     for n_fmcw = 1:1:N_int_FMCW
%         distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw).*(t_int+int_FMCW *T_int_FMCW);  %CW干扰距离更新
%         td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                       %干扰时延更新   *注意这里时延不是两倍了 干扰是直接到达接收端的
% %         S_int_FMCW{n_fmcw,1} = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW干扰信号
%         S_int_FMCW{n_fmcw,1} = AT*cos(2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2));      %FMCW干扰信号 *实信号
%         S_int_FMCW_targetsum = S_int_FMCW_targetsum + S_int_FMCW{n_fmcw,1};
%     end
%     
% %     S_int_FMCW = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW干扰信号
%         S_int_FMCW_sum(int_FMCW*N_FMCWint+1:(int_FMCW+1)*N_FMCWint) = S_int_FMCW_targetsum;
% end



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
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2)); %发射信号
    St = AT*cos(2*pi*(fc*(t-i*T_chirp)+slope/2*t.^2)); %发射信号
    
    Sr_sum = zeros(1,N_ADC);                              %一个chirp对应回波信号叠加空矩阵 1024点
    S_int_sum = zeros(1,N_ADC);                           %干扰信号空矩阵 1024点
    
%     %% CW干扰信号
%     for int_CW = 1:1:3 %干扰个数循环
%         
% %         S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %正弦干扰信号(对应不同位置的干扰信号)
%         S_int_CW{int_CW,1} = 10*cos(2*pi*fint_CW(int_CW)*(t-td_int_CW{int_CW,1})); %CW干扰 *实信号
%         S_int_sum = S_int_sum + S_int_CW{int_CW,1};
%         
%     end
    
    
    %% 生成没有干扰的回波信号
    for m = 1:1:1
        
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %目标回波信号
        Sr{m,1} = A_target(m)*AT*cos(2*pi*(fc*(t-i*T_chirp-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %目标回波信号
        Sr_sum = Sr_sum + Sr{m,1};
        
    end
    
    Sr_sum_withCW = zeros(1,N_ADC);
    Sr_sum_addCW = zeros(1,N_ADC);
    
     %% CW干扰信号
    for int_CW = 1:1:N_int_CW %干扰个数循环
        
%         S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %正弦干扰信号(对应不同位置的干扰信号)
        S_int_CW{int_CW,1} = 20*cos(2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %CW干扰 *实信号
%         S_int_sum = S_int_sum + S_int_CW{int_CW,1};
        Sr_sum_addCW = [Sr_sum(1:n_int_start(int_CW)-1),Sr_sum(n_int_start(int_CW):n_int_end(int_CW))+S_int_CW{int_CW,1}(n_int_start(int_CW):n_int_end(int_CW)),Sr_sum(n_int_end(int_CW)+1:N_ADC)];
        Sr_sum_withCW = Sr_sum_withCW + Sr_sum_addCW;
        
    end
    
    
    %% ----添加噪声模块-----
    %% CW噪声
%     if i <= 80
%         
% %         Sr_sum_withCW = Sr_sum + S_int_sum;
%         Sr_sum_withCW = [Sr_sum(1:n_int_start-1),Sr_sum(n_int_start:n_int_end)+S_int_sum(n_int_start:n_int_end),Sr_sum(n_int_end+1:N_ADC)];
%         
%     else
%         Sr_sum_withCW = Sr_sum; 
%     end
    
    %% FMCW噪声

%     Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %加入干扰
    
    %% 混频得到中频
%     IF = St .*conj(Sr_sum);            %无干扰
%     IF = St .*conj(Sr_sum_withCW);     %CW噪声
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW噪声
    
    %% 混频得到中频
%     IF = St .* Sr_sum;
%     IF = IF + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));           %无干扰
    
%     IF = St .* Sr_sum_withFMCW; 
%     IF = IF + 1i*Sr_sum_withFMCW .*sin(2*pi*(fc*t+slope/2*t.^2));  %FMCW干扰
    
    IF_I = St .*Sr_sum_withCW;                                       %接收机混频得到I分量
    IF = IF_I + 1i*Sr_sum_withCW .*sin(2*pi*(fc*t+slope/2*t.^2));    %宽带正弦干扰
%     IF = hilbert(IF_I);
    IF_I_mat(i+1,:) = IF_I;                                          %* 接收机中频信号I分量 用于时频分析
    
%     IF = St .*conj(Sr_sum);                                        %复信号建模情况  *无干扰                                
    
    IF_wonoise = St .*Sr_sum + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));  %* 无噪声中频信号进行时频分析用
    IF_wonoise_mat(i+1,:) = IF_wonoise;                                    
    IF_wonoise_I_mat(i+1,:) = St .*Sr_sum;                                 %* 无噪声中频信号I分量 进行时频分析用

   

    IF_mat(i+1,:) = IF;              %用于时频分析
    IF_mat_withnoise(i+1,:) = IF;    %将带有噪声的中频信号保存
    St_mat(i+1,:) = St;
    Sr_mat(i+1,:) = Sr_sum;
    
end


%% 画中频信号时域图
figure(1);
IF_realtime = IF_I_mat(1,:); %取第一个chirp的中频采样信号
plot(IF_realtime)
title('一个chirp对应的中频信号时域采样')
xlabel('采样点数N');

%% 画中频信号频谱图
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %取中频矩阵一行1024点
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure(2);
plot(f_IF,IF_realfft);
title('中频信号fft')
xlabel('frequency(Hz)');


%% 信号时频分析 短时傅里叶变换STFT(暂定)
[tfr_IF,t_IF , f_IF_I]=  tfrstft(conj(IF_wonoise_mat(1,1:N_ADC)'));
figure(3);
imagesc(t_IF,f_IF_I(1:length(f_IF_I)/2)*Fs,abs(tfr_IF(1:N_ADC/2,:)));
title('无噪声中频信号短时傅里叶时频图');

[tfr_St,t_st,f_st] = tfrstft(conj(IF_mat(1,1:N_ADC)'));
figure(4);
imagesc(t_st,f_st(1:length(f_st)/2)*Fs,abs(tfr_St(1:N_ADC/2,:)));
title('受干扰中频信号短时傅里叶时频图');

figure(5);
tfr_IF = tfrstft(conj(IF_mat(1,1:N_ADC)'));
imagesc(abs(tfr_IF));


% figure;
% [B_IF,t_IF,f_IF_I]=tfrstft(conj(IF_mat(1,1:N_ADC)'));             %复数矩阵 若直接转置里面的符号也变了
% % imagesc(t_IF,f_IF_I(N_ADC/2:N_ADC-1)*Fs,abs(B_IF(N_ADC/2:N_ADC-1,:)));
% imagesc(t_IF/Fs,f_IF_I(1:N_ADC/2)*Fs,abs(B_IF));


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