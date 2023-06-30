% --2022/10/21-- %
% --�-- %
% FMCW����3.0.ver
%% �ɵ���Ŀ�����
%% ���ø��ź���ʽ
%% ������Ҹ��� �ٶȵ�Ԫ���ָ������� ���ٶȵ�Ԫ�ϵļ�����Ӱ��

clear all;
close all;
clc;

%% Ŀ���������
R_target = [45];
V_target = [3.5];
A_target = [0.8 0.1 0.5]; %Ŀ���Ӧ�Ľ�������
N_target = length(R_target);

%% FMCW�״�����������
c = 3e8;
fc = 77e9;                %chirp��ʼƵ��
fc_noise = 77.6e9;        %���Ҹ����ź�Ƶ��
B_chirp = 1.6e9;          %chirp���� ��Դʦ��ĳ�������
T_chirp = 62.48e-6;       %chirp����ʱ��
N_ADC = 1024;             %chirp�����ڲ�������
N_chirp = 128;            %chirp����
Fs = N_ADC/T_chirp;       %������ 
t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
slope = B_chirp/T_chirp;  %chirp��Ƶб��
lambda = c / fc ;         %����
IF_store = zeros(N_chirp,N_ADC);

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(2*slope);             %�״�̽��������
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���

%% �����źŲ���
AT = 1;                  %�����ź�����
t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
%AR = 0.8;                 %�ز��ź�˥���ı���ֵ
Anoise = 10;              %�����ź�����
fnoise = 78e9;            %�����źŹ���
%t = t - Tr/2; %��fc��Ϊ����Ƶ��

%% ���䡢������������
for i = 0:1:N_chirp-1  %chirp����ѭ��
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
    St = AT*exp((1i*2*pi)*(fc*(t-i*T_chirp)+slope/2*t.^2));    %�����ź�
    Sn = AT*exp(1i*2*pi*fc_noise*(t-i*T_chirp));               %�����ź�
    slope_N = 1.09e13;
    SN = AT*exp((1i*2*pi)*(fc_noise*(t-i*T_chirp)+slope_N/2*t.^2)); 
    Snconj = conj(Sn);                                         %�����źŹ���
    Sr_sum = zeros(1,N_ADC);                                   %�ز��źŵ��ӿվ���
    
    %% ���ɻز��ź�
    for m = 1:1:N_target
        Sr{m,1} = A_target(m)*AT*exp((1i*2*pi)*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2));
        Srconj{m,1} = conj(Sr{m,1});
        Sr_sum = Sr_sum + Srconj{m,1};
    end
    
    %% ��Ƶ�õ���Ƶ
    IF = St .* Sr_sum ;
%     IF = St .* (Sr_sum+Sn);    %��Ƶ�ź��и���ʱ
%     SNR = 10;                      %�����
%     IF = awgn(IF,SNR,'measured');  %����Ƶ�źżӸ�˹�������������������ʱ��Ҫ���������Ĳ���
    IF_mat_1(i+1,:) = IF;            %��������������Ƶ�źű���
end
save('Ego_vehicle_���ź�.mat', 'IF_mat_1');

%% ����������
range_win = hanning(N_ADC);      %����range��
doppler_win = hanning(N_chirp);  %����doppler��

%% range fft
for i = 1:1:N_chirp
    temp = IF_mat_1(i,:) .* range_win';
    temp_fft = fft(temp,N_ADC)./N_ADC;
    IF_mat_1(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N_ADC
    temp = IF_mat_1(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,N_chirp));
    IF_mat_1(:,j) = temp_fft;
end

IF_mat_1 = abs(IF_mat_1);
RDM = 10*log10(IF_mat_1);

%% ��ͼ
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
