% --2022/10/21-- %
% --�-- %
% FMCW����3.0.ver

%% �˳�����Ҫ����idel_time���ٶȲ�����Ӱ��
%% �ɵ���Ŀ�����
%% ���ø��ź���ʽ
%% �������ʱ��idel_time

clear all;
close all;
clc;

%% Ŀ���������
R_target = [1.3 7 60];
V_target = [0.5 3 10];
N_target = length(R_target);

%% FMCW�״�����������
c = 3e8;
fc = 77e9;                %chirp��ʼƵ��
B_chirp = 300e6%1.6e9;          %chirp���� ��Դʦ��ĳ�������
T_chirp = 30e-6%62.48*10^(-6);  %chirp����ʱ��
%N_ADC = 512%1024;             %chirp�����ڲ�������
N_chirp = 256%128;            %chirp����
Fs = 16.67e6;%N_ADC/T_chirp;       %������ 
N_ADC = floor(T_chirp * Fs);
t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
slope = B_chirp/T_chirp;  %chirp��Ƶб��
lambda = c / fc ;         %����
IF_store = zeros(N_chirp,N_ADC);

%% IDEL TIME����(����ʱ���źŽ�ģ)
idel_time = 4e-6;                   %idel_time����ʱ��
N_idel = floor(idel_time * Fs);     %�������
t_zeros = 0:1/Fs:idel_time-1/Fs;    %ȷ�����������idel time�е�ʱ��λ��
idel_zeros = 0 * t_zeros;           %�������

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(2*slope);             %�״�̽��������
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���

%% �����źŲ���
AT = 10;                  %�����ź�����
t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
AR = 0.8;                 %�ز��ź�˥���ı���ֵ
Anoise = 10;              %�����ź�����
fnoise = 78e9;            %�����źŹ���
%t = t - Tr/2; %��fc��Ϊ����Ƶ��

%% ���䡢������������
for i = 0:1:N_chirp-1  %chirp����ѭ��
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp+i*idel_time); %�������
        t_d{k,1} = 2*distant{k,1}/c;                                         %ʱ�Ӹ���
    end
    
    St = AT*exp((1i*2*pi)*(fc*(t-i*T_chirp-i*idel_time)+slope/2*t.^2));   %�����ź�
    St_full = [St,idel_zeros];                                            %����ʱ���ڲ���
    Sr_sum = zeros(1,N_ADC);                                              %�ز��źŵ��ӿվ���
    Srfull_sum = zeros(1,N_ADC + N_idel);
    
    %% ���ɻز��ź�
    for m = 1:1:N_target
        Sr{m,1} = AR*AT*exp((1i*2*pi)*(fc*(t-t_d{m,1}-i*T_chirp-i*idel_time)+slope/2*(t-t_d{m,1}).^2));
        Srfull{m,1} = [Sr{m,1},idel_zeros];
        Srfull_conj{m,1} = conj(Srfull{m,1});
        Srconj{m,1} = conj(Sr{m,1});
        Sr_sum = Sr_sum + Srconj{m,1};
        Srfull_sum = Srfull_sum +Srfull_conj{m,1};
    end
%     Sr_itf = Anoise*exp(1i*2*pi*fnoise*(t+i*T_chirp));                   
%     Sr_itf_conj = conj(Snoise);                                           %���Ҹ����źŹ���
    
    %% ��Ƶ�õ���Ƶ
    IF = St .* Sr_sum;
    IF_full = St_full .* Srfull_sum;
%     IF = St .* (Sr_conj+Snoise_conj);         %��Ƶ�ź��и���ʱ
    SNR = 10;                      %�����
    IF = awgn(IF,SNR,'measured');  %����Ƶ�źżӸ�˹�������������������ʱ��Ҫ���������Ĳ���
    IF_mat(i+1,:) = IF;          %��������������Ƶ�źű���
    IFfull_mat(i+1,:) = IF_full;
end


%% ����������
range_win = hanning(N_ADC);      %����range��
doppler_win = hanning(N_chirp);  %����doppler��

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

% %% ����������
% range_win = hanning(N_ADC);      %����range��1124��
% doppler_win = hanning(N_chirp);  %����doppler��
% 
% IFfull_mat1 = zeros(N_chirp,512); %���ڴ�Ŷ�άfftֵ

% %% range fft
% for i = 1:1:N_chirp
%     temp = IF_mat(i,:) .* range_win';
%     temp_fft = fft(temp,512); %���ǽ���1024��任
%     IFfull_mat1(i,:) = temp_fft;
% end
% %% doppler fft
% for j = 1:1:512
%     temp = IFfull_mat1(:,j) .* doppler_win;
%     temp_fft = fftshift(fft(temp,N_chirp));
%     IFfull_mat1(:,j) = temp_fft;
% end

%% ��ͼ
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