%--2022/10/2--%
%--�--%
%������Լ�FMCW���״���١������з��档�������״��źŷ�����endle_time ��֤�õ���ƵΪ��Ƶ�źŽ��д���Ϊ�������п����ŷ����̵�

%% ���⣺
%% 1.matlabд����ʱ��Ҫ��������﷢���ź�������źţ���ôʱ����Ӧ�ü�������  2022/10/6 �ѽ����
%% 2.ʱ���ӳ��Ƿ�Ҫ��ɱ������˳���д�ɱ���td(i)������ʦ��ĳ��������Ƕ�ֵtd.TI�ĵ���Ҳ���ƺ���IF�ź�Ƶ���������ٶ����������ɶ�ֵ 2022/10/6 �ѽ����
%% 3.��һ��main.m��������� ������ҲҪ�����������������������ƥ������
clear all;
close all;
clc;

%% �״��������
% maxR = 200;           % �״����̽��Ŀ��ľ���
% rangeRes = 1;         % �״�ľ������
% maxV = 70;            % �״������Ŀ����ٶ�
fc= 77e9;             % �״﹤��Ƶ�� ��Ƶ
c = 3e8;              % ����

%% Ŀ���������
r0 = 150; % Ŀ��������� (max = 200m)
v0 = 50; % Ŀ���ٶ����� (min =-70m/s, max=70m/s)


%% �״﷢��FMCW��������
T_chirp=7.33e-6;       %chirp����
B_chirp=150e6;         %�����źŵ�Ƶ����
slope=B_chirp/T_chirp; %���Ե�Ƶб��68MHz/us


%  slope = 29.306e12;
%  T_chirp = 29.56e-6; %chirp����
%  B_chirp = 750e6;    %�����źŵ�Ƶ����
%  fs=20e6 
N_chirp=128;                          %chirp���� 
N_ADC=1024;                           %ADC��������
fs=N_ADC/T_chirp;                     %��Ƶģ���źŲ����� 

t=linspace(0,N_chirp*T_chirp,N_ADC*N_chirp); %�����źźͽ����źŵĲ���ʱ�䣬��MATLAB�е�ģ���ź���ͨ�������ź����޲�������
Tx=zeros(1,length(t));                       %�����ź�
Rx=zeros(1,length(t));                       %�����ź�
Mix = zeros(1,length(t));                    %��Ƶ�����ġ���Ƶ����Ƶ�ź�

r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% ��Ŀ���ź�����
fint = fc;
Tint = N_chirp*T_chirp;
Bint = B_chirp;
slope_int = Bint/Tint;

for i=1:length(t) 
    
    r_t(i) = r0 + v0*t(i); % �������
    td(i) = 2*r_t(i)/c;    % �ӳ�ʱ��
    
    Tx(i) = cos(2*pi*(fc*t(i) + (slope*t(i)^2)/2));                 %�����ź� ʵ���ź�
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2/2))); %�����ź� ʵ���ź�
    int(i) = cos(2*pi*(fint*t(i) + slope_int*t(i)^2/2));            %FMCW�����ź�
    Rxint(i) = Rx(i) + int(i);                                      %�Ӹ���
    
    if i<=1024
         freq(i)=fc+slope*i;       %�����ź�ʱƵͼ ֻȡ��һ��chirp
         freq_echo(i)=fc+slope*i;  %�ز��ź�Ƶ���ӳ�
    end

    Mix(i) = Tx(i).*Rx(i);        %��Ƶ�����ġ���Ƶ����Ƶ�ź�
    Mixint(i) = Tx(i).* Rxint(i); %�������źŵ���Ƶ�����ź�
end
Tx;
Rx;
Rxint;

%% ʱƵ���� ���ڽ�����
sigint = reshape(Mixint,N_ADC,N_chirp);
wlen = 64; %��������С
hop = wlen/8;
nfft =1024;
win = blackman(wlen,'periodic');
% [S_if,F_if,T_if] = spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
figure;
spectrogram(sigint(:,128), win, wlen - hop, nfft, fs);
title('�ܸ�����Ƶ�źŶ�ʱ����ҶʱƵͼ');


%�����ź�ʱ��ͼ
figure;
plot(Tx(1:1024));
xlabel('����');
ylabel('����');
title('TX�����ź�ʱ��ͼ');
hold on;
%�����ź�ʱ��ͼ
plot(Rx(1:1024));
xlabel('����');
ylabel('����');
title('RX�����ź�ʱ��ͼ');

% %�����ź�ʱƵͼ
% figure;
% plot(t(1:1024),freq);
% xlabel('ʱ��');
% ylabel('Ƶ��');
% title('TX�����ź�ʱƵͼ');


% %�����ź�ʱ��ͼ
% figure;
% plot(Rx(1:1024));
% xlabel('����');
% ylabel('����');
% title('RX�����ź�ʱ��ͼ');

%�����ź��뷢���źŵ�ʱƵͼ
% figure;
% plot(t(1:1024),freq);
% hold on;
% plot(td(1:1024)+t(1:1024),freq);
% xlabel('ʱ��');
% ylabel('Ƶ��');
% title('�����ź��뷢���ź�ʱƵͼ');
% legend ('TX','RX');

%��Ƶ�ź�Ƶ�� ��Ƶ�źŹ۲�
%figure;
% plot(db(abs(fft(Mix(1:1024*256)))));%�鿴����ĺ�Ƶ�ź� ��chirp�ĵ�����Ϊ1024*256���ɿ�����һ�����źţ���ע�������ڴ档
% xlabel('Ƶ��');
% ylabel('����');
% title('��Ƶ�ź�Ƶ��');

% figure;
% plot(db(abs(fft(Mix(1:1024)))));%�鿴����ĺ�Ƶ�ź� ��chirp�ĵ�����Ϊ1024*256���ɿ�����һ�����źţ���ע�������ڴ档ֻ��Ե�һ��chirp��Ӧ����Ƶ1024��
% xlabel('Ƶ��');
% ylabel('����');
% title('��Ƶ�ź�Ƶ��');

%% ��ͨ�˲� ��ֹƵ��30MHz  ����Ƶ��120MHz
% Mix=lowpass(Mix(1:1024*256),30e6,120e6);
% plot(db(abs(fft(Mix(1:1024*256)))));
% xlabel('Ƶ��');
% ylabel('����');
% title('��Ƶ�źŵ�ͨ�˲���');

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
signal = reshape(Mix,N_ADC,N_chirp); %*�Ի�Ƶ�źž�������ع������1024��128�ľ��󣬷���������о���-�ٶ�άFFT

figure;
mesh(signal);
xlabel('������')
ylabel('��������');
title('��Ƶ�ź�ʱ��');

%% ����άFFT
sig_fft = fft(signal,N_ADC)./N_ADC;
sig_fft = abs(sig_fft);
sig_fft = sig_fft(1:(N_ADC/2),:);

figure;
plot(sig_fft(:,1));
xlabel('���루Ƶ�ʣ�');
ylabel('����')
title('����chirp��FTF���')


%% ����FFT����׾���
figure;
mesh(sig_fft);
xlabel('���루Ƶ�ʣ�');
ylabel('chirp������')
zlabel('����')
title('����άFTF���')

%% �ٶ�άFFT

Mix=reshape(Mix,[N_ADC,N_chirp]); %*MIX�źž�������ع�
sig_fft2 = fft2(Mix,N_ADC,N_chirp);

sig_fft2 = sig_fft2(1:N_ADC/2,1:N_chirp);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
%RDM = 10*log10(RDM) ;
doppler_axis = linspace(-100,100,N_chirp);
range_axis = linspace(-200,200,N_ADC/2)*((N_ADC/2)/400);

figure;
mesh(doppler_axis,range_axis,RDM);
xlabel('������ͨ��'); ylabel('����ͨ��'); zlabel('���ȣ�dB��');
title('�ٶ�άFFT �����������');