
%实测数据中的一部分
%解决了正负距离的问题

% R=1;  %距离
c=3*10^8;   %光速
% td=2*R/c;   %回波信号的时延
tChirp=130.57*10^(-6);   %chirp的持续时间
k=24.043*10^12; %斜率，已改

B=k * tChirp; %带宽
fc=77*10^9; %载波频率，已改
N=256;  %采样点，已改
fs=10000*10^3; %采样频率，已改
tADC=(N - 1)/fs;  %采样持续时间，已改
tStart=5*10^(-6);  %采样开始时间，已改
t=tStart+(tADC/N)/2:tADC/N:tStart+tADC-(tADC/N)/2;


time = 0:1/fs:tADC;
x_r = result_r(1,1:256);
plot(time,x_r,'b');
hold on;
x_i = result_i(1,1:256);
plot(time,x_i,'r');
hold on;

din = result(1, 1:256); %result中的第一行1：256点数据赋值给din
% din = x_r + 1j* x_i; % complex value
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


% index = 1:1:N;
% freq_bin = (index - 1) * fs / N;
% range_bin = freq_bin * c / 2 / k;

% figure;
% plot(freq_bin,abs(temp_fft));
% title('after fft');

% range_win = hamming(N); 
% range_win = range_win';
% din_win = din .* range_win;
%  
% din_win_fft = fft(din_win);

Pout_fft = 20*log10(temp_fft);
b_adc = 16;
cf = 20*log10(2^(b_adc-1))+20*log10(wSum)-20*log10(sqrt(2));
Pout_dBFs = Pout_fft-cf;

figure
plot(range_bin, Pout_dBFs)
% yticks([-140 -120 -100 -80 -60 -40 -20 0])
title('fft(after window)')
xlabel('range(m)')


