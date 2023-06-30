
%ʵ�������е�һ����
%������������������

% R=1;  %����
c=3*10^8;   %����
% td=2*R/c;   %�ز��źŵ�ʱ��
tChirp=130.57*10^(-6);   %chirp�ĳ���ʱ��
k=24.043*10^12; %б�ʣ��Ѹ�

B=k * tChirp; %����
fc=77*10^9; %�ز�Ƶ�ʣ��Ѹ�
N=256;  %�����㣬�Ѹ�
fs=10000*10^3; %����Ƶ�ʣ��Ѹ�
tADC=(N - 1)/fs;  %��������ʱ�䣬�Ѹ�
tStart=5*10^(-6);  %������ʼʱ�䣬�Ѹ�
t=tStart+(tADC/N)/2:tADC/N:tStart+tADC-(tADC/N)/2;


time = 0:1/fs:tADC;
x_r = result_r(1,1:256);
plot(time,x_r,'b');
hold on;
x_i = result_i(1,1:256);
plot(time,x_i,'r');
hold on;

din = result(1, 1:256); %result�еĵ�һ��1��256�����ݸ�ֵ��din
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


