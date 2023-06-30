% --2022/10/25-- %
% --�-- %
% FMCW����5.0.ver

%% �ɵ���Ŀ�����
%% �״���������� ����1.6G��chirp����62.48us��ʦ��ĳ���  ��
%% ��ʵ�źŽ�ģ ����������������������ϣ�����ر任ʵ�֣�    ��
%% ÿһ��chirp���ڶ��������ʱ�� idle_time               ��
%% 

clear all;
close all;
clc;

%% Ŀ���������
R_target = [45];
V_target = [3.5];
N_target = length(R_target);

%% FMCW�״�����������
c = 3e8;
fc = 77e9;                %chirp��ʼƵ��
B_chirp = 1.6e9           %chirp���� ��Դʦ��ĳ�������
T_chirp = 62.48*10^(-6);  %chirp����ʱ��
N_ADC = 1024;             %chirp�����ڲ�������
N_chirp = 128;            %chirp����
Fs = N_ADC/T_chirp;       %������ 
t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
slope = B_chirp/T_chirp;  %chirp��Ƶб��
lambda = c / fc ;         %����
IF_store = zeros(N_chirp,N_ADC);

%% IDEL TIME����(����ʱ���źŽ�ģ)
N_idel = 0;                       %�������
idel_time = N_idel / Fs;            %idel_time����ʱ��
t_zeros = 0:1/Fs:idel_time-1/Fs;    %ȷ�����������idel time�е�ʱ��λ��
idel_zeros = 0 * t_zeros;           %�������

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(2*slope);             %�״�̽��������
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���

%% �����źŲ���
AT = 10;                  %�����ź�����
% t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
% t = 0:1/Fs:T_chirp;
AR = 0.8;                 %�ز��ź�˥���ı���ֵ
Anoise = 10;              %�����ź�����
fnoise = 78e9;            %�����źŹ���

%% ���䡢������������
for i = 0:1:N_chirp-1  %chirp����ѭ��
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp) + i*idel_time*V_target(k); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
%     St = AT*exp((1i*2*pi)*(fc*(t+i*T_chirp)+slope/2*t.^2));   %�����ź�
    St = AT*cos(2*pi*(fc*(t-i*T_chirp-i*idel_time)+slope/2*t.^2));   %�����ź�
    St_full = [St,idel_zeros];                             %����ʱ���ڲ���
    Sr_sum = zeros(1,N_ADC);                                   %�ز��źŵ��ӿվ���
    Srfull_sum = zeros(1,N_ADC + N_idel);
    
    %% ���ɻز��ź�
    for m = 1:1:N_target
        Sr{m,1} = AR*AT*cos(2*pi*(fc*(t-t_d{m,1}-i*T_chirp-i*idel_time)+slope/2*(t-t_d{m,1}).^2));
        Srfull{m,1} = [Sr{m,1},idel_zeros];
        Sr_sum = Sr_sum + Sr{m,1};
        Srfull_sum = Srfull_sum +Srfull{m,1};
    end
    
%     Sr_itf = Anoise*exp(1i*2*pi*fnoise*(t+i*T_chirp));                   
%     Sr_itf_conj = conj(Snoise);                                           %���Ҹ����źŹ���
    
    %% ��Ƶ�õ���Ƶ
    IF_I = St .* Sr_sum;                      %��Ƶʵ�ź�
    IF_Q = imag(hilbert(IF_I));               %ϣ�����ر任�õ���Ƶ�źŵ�Q����      
    IF = IF_I + 1i*IF_Q;
%     IF = St .* (Sr_conj+Snoise_conj);         %��Ƶ�ź��и���ʱ
    SNR = 10;                                 %�����
    IF = awgn(IF,SNR,'measured');             %����Ƶ�źżӸ�˹�������������������ʱ��Ҫ���������Ĳ���
    IF_mat(i+1,:) = IF;                       %��������������Ƶ�źű���
    IF_I_mat(i+1,:) = IF_I;
    IF_Q_mat(i+1,:) = IF_Q;
    
    %% ������ȫ�ź�
    IF_full_I = St_full .* Srfull_sum;
    IF_full_Q = imag(hilbert(IF_full_I));
    IF_full = IF_full_I + 1i*IF_full_Q;
    IF_full_mat(i+1,:) = IF_full;
    IF_full_mat_I(i+1,:) = IF_full_I;
    
end
IF_full_tdomain = reshape(IF_full_mat_I',1,N_chirp*(N_ADC+N_idel));

%% ����������
range_win = hanning(N_ADC);      %����range��
doppler_win = hanning(N_chirp);  %����doppler��

%% range fft
for i = 1:1:N_chirp
    temp = IF_I_mat(i,:) .* range_win';
    temp_fft = fft(temp,N_ADC);
    IF_I_mat(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N_ADC
    temp = IF_I_mat(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,N_chirp));
    IF_I_mat(:,j) = temp_fft;
end

IF_I_mat = abs(IF_I_mat);
RDM = 10*log10(IF_I_mat);

%% ��ͼ
figure;
mesh(IF_full_mat_I);
xlabel('fast-time');
ylabel('slow_time');
zlabel('Amplitude');
title('idel time��Ƶ�ź�����');

figure;
plot(IF_full_tdomain);
xlabel('samples');
ylabel('Amplitude');
title('idel time��Ƶ�ź�ʱ��');

figure;
mesh(IF_I_mat);
xlabel('fast-time');
ylabel('slow_time');
zlabel('Amplitude');
title('��Ƶ�����ź�ʱ��');

figure;
plot(IF_I_mat(1,:));
hold on
plot(IF_Q_mat(1,:));
xlabel('samples');
ylabel('Amplitude');
title('I-Q·�ź�');
legend('I·','Q·');

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