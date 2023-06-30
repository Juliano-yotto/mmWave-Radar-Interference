clear all;
close all;
clc;
%% �״�ϵͳ��������
fc= 77e9;             % �״﹤��Ƶ�� ��Ƶ
c = 3e8;              % ����

%% FMCW���β�������
maxR = 200;           % �״����̽��Ŀ��ľ���
rangeRes = 1;         % �״�ľ������ 
maxV = 20;            % �״������Ŀ����ٶ�
B = c/(2*rangeRes);       % �����źŴ��� (y-axis)  B = 150MHz
N_chirp=128;                          %chirp���� 
N_ADC=1024;                        %ADC��������
Tchirp = 5 * 2 * maxR/c;  % ɨƵʱ�� (x-axis), 5.5= sweep time should be at least 5 o 6 times the round trip time
endle_time=2.5e-6;          %����ʱ��

vres=(c/fc)/(2*N_chirp*(Tchirp+endle_time));%�ٶȷֱ���
%% �û��Զ���Ŀ�����
r0 = 20;  % Ŀ����� rangeRes*����ͨ���� ����ͨ���ŷ�Χ[1 256]
v0 = 12        % Ŀ���ٶ� vres*�ٶ�ͨ����[1 128]

slope = B / Tchirp;         %��Ƶб��
f_IFmax= (slope*2*maxR)/c ; %�����ƵƵ��
f_IF=(slope*2*r0)/c ;       %��ǰ��ƵƵ��
% Nr=1024*256;                %��Ƶ�źŵ�������
Fs=N_ADC/Tchirp;                 %ģ���źŲ���Ƶ��

%%
t=linspace(0,N_chirp*Tchirp,N_ADC*N_chirp); %�����źźͽ����źŵĲ���ʱ�䣬��MATLAB�е�ģ���ź���ͨ�������ź����޲������ɵġ�
Tx=zeros(1,length(t)); %�����ź�
Rx=zeros(1,length(t)); %�����ź�
Mix = zeros(1,length(t)); %��Ƶ�����ġ���Ƶ����Ƶ�ź�
r_t=zeros(1,length(t));
td=zeros(1,length(t));
%% ��Ŀ���ź�����
for i=1:length(t)
    r_t(i) = r0 + v0*t(i); % �������
    td(i) = 2*r_t(i)/c;    % �ӳ�ʱ��
    
    Tx(i) = cos(2*pi*(fc*t(i) + (slope*t(i)^2)/2)); % �����ź� ʵ���ź�
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2)/2)); %�����ź� ʵ���ź�
    
    if i<=1024
         freq(i)=fc+slope*i*Tchirp/N_ADC; %�����ź�ʱƵͼ ֻȡ��һ��chirp
         freq_echo(i)=fc+slope*i*Tchirp/N_ADC;%�ز��ź�Ƶ���ӳ�
    end

    Mix(i) = Tx(i).*Rx(i);%��Ƶ�����ġ���Ƶ����Ƶ�ź�
end

signal = reshape(Mix,N_ADC,N_chirp);

MixIQ=zeros(N_ADC,N_chirp);
for i=1:N_chirp

MixIQ(:,i) =hilbert(Mix(((i-1)*N_ADC+1):N_ADC*i)); % ϣ�����ر任 ��ΪIQ�����ź�

end

sig_fft1 = zeros(N_ADC,N_chirp);
for k=1:N_chirp
    sig_fft1(:,k)=fft(MixIQ(:,k));
end

sig_fft = abs(sig_fft1);
figure;
plot(sig_fft(:,1));
xlabel('���루Ƶ�ʣ�');
ylabel('����')
title('��һ��chirp��FFT���')

%% ����FFT����׾���
figure;
mesh(sig_fft);
ylabel('���루Ƶ�ʣ�');
xlabel('chirp������')
zlabel('����')
title('����άFTF���')

%% �ٶ�άFFT
sig_fft2 = zeros(N_ADC,N_chirp);
for k=1:N_ADC
    sig_fft2(k,:)=fft(sig_fft1(k,:));
end

% sig_fft2 = fft2(MixIQ);
RDM = abs(sig_fft2);
% RDM =10*log10(RDM);
doppler_axis = linspace(0,128,N_chirp)*vres;
range_axis = linspace(0,256,N_ADC)*rangeRes;

figure;
mesh(doppler_axis,range_axis,RDM);
xlabel('������ͨ��'); ylabel('����ͨ��'); zlabel('���ȣ�dB��');
title('�ٶ�άFFT �����������');

