% --2022/10/25-- %
% --�-- %
% FMCW����5.0.ver

%% �ɵ���Ŀ�����
%% �״���������� ����1.6G��chirp����62.48us��ʦ��ĳ���  ��
%% ��ʵ�źŽ�ģ ����������������������ϣ�����ر任ʵ�֣�    ��

clear all;
close all;
clc;

%% Ŀ���������
R_target = [45];
V_target = [3];
N_target = length(R_target);

%% FMCW�״�����������
c = 3e8;
fc = 77e9;                %chirp��ʼƵ��
B_chirp = 1.6e9;          %chirp���� ��Դʦ��ĳ�������
T_chirp = 62.48*10^(-6);  %chirp����ʱ��
N_ADC = 1024;             %chirp�����ڲ�������
N_chirp = 128;            %chirp����
Fs = N_ADC/T_chirp;       %������ 
% t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
t = 0:1/Fs:(N_ADC-1)/Fs;  %**����
slope = B_chirp/T_chirp;  %chirp��Ƶб��
lambda = c / fc ;         %����
IF_store = zeros(N_chirp,N_ADC);

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(2*slope);             %�״�̽��������
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���

%% �����źŲ���
AT = 10;                  %�����ź�����
AR = 0.8;                 %�ز��ź�˥���ı���ֵ
Anoise = 10;              %�����ź�����
fnoise = 78e9;            %�����źŹ���

%% ���Ų�������
fint = [77.2e9 77.3e9 77.5e9 77.55e9 77.6e9  78e9 78.5e9]; %��Ƶ�������������Ҳ� Ƶ��
R_int = [1.5 20 10 7 8 45 50 ];                            %���ŵ�λ��
V_int = [1.5 8 2 10 12 4 5 ];                              %�����ٶ�
N_int = length(fint);                                      %���Ÿ���

%% ���䡢������������
for i = 0:1:N_chirp-1  %chirp����ѭ��
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2));   %�����ź�
    St = zeros(1,N_ADC);
    Sr_sum = zeros(1,N_ADC);                 %�ز��źŵ��ӿվ���
    St = AT*cos(2*pi*(fc*t+slope/2*t.^2));   %�����ź�
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2));
    
    %% ���Ŷ�������Ϣ
    for n = 1:1:N_int
        distant_int{n,1} = R_int(n) +V_int(n).*(t+i*T_chirp);    %���ž������
        td_int{n,1} = 2*distant_int{n,1}/c;                      %����ʱ�Ӹ���
    end
    
    %% ���ɸ����ź�
    S_int_sum = zeros(1,N_ADC);                           %�����źſվ��� 1024��
    
    for int = 1:1:N_int %���Ÿ���ѭ��
        
        S_int{int,1} = AT*exp(1i*2*pi*fint(int)*(t-i*T_chirp-td_int{int,1})); %���Ҹ����ź� 
%         S_int{int,1} = 2*AT*exp(1i*2*pi*(fc*(t-i*T_chirp-td_int{int,1})+slope/2*(t-td_int{int,1}).^2)); %���Ҹ����ź�
        S_int_sum = S_int_sum + S_int{int,1};
    
    end
   
    %% ���ɻز��ź�
    for m = 1:1:N_target
%         Sr{m,1} = AR*AT*cos(2*pi*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2));
        Sr{m,1} = AR*AT*cos(2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));
%         Sr{m,1} = AR*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));
        Sr_sum = Sr_sum + Sr{m,1};
    end
    
%     Sr_itf = Anoise*exp(1i*2*pi*fnoise*(t+i*T_chirp));                   
%     Sr_itf_conj = conj(Snoise);                                           %���Ҹ����źŹ���
    
    %% ��Ƶ�õ���Ƶ
%     IF_I = St .* (Sr_sum + S_int_sum);                      %��Ƶʵ�ź�
    IF_I = St .* Sr_sum ;                                   %�޸�����Ƶʵ�ź�
%     IF = IF_I + 1i*IF_I;

    %% ��ͨ�˲�
%     IF_I = filter(Hd,IF_I);
    
    %% ��ɸ��ź�
    IF_Q = imag(hilbert(IF_I));               %ϣ�����ر任�õ���Ƶ�źŵ�Q����      
    IF = IF_I + 1i*IF_Q;
%     IF = St .* (Sr_conj+Snoise_conj);         %��Ƶ�ź��и���ʱ
%     SNR = 10;                                 %�����
%     IF = awgn(IF,SNR,'measured');             %����Ƶ�źżӸ�˹�������������������ʱ��Ҫ���������Ĳ���
    IF_mat(i+1,:) = IF_I;                       %��������������Ƶ�źű���
    IF_I_mat(i+1,:) = IF_I;
    IF_Q_mat(i+1,:) = IF_Q;     
end

save('Ego_vehicle_���ź�.mat', 'IF_mat');

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


IF_mat = abs(IF_mat(:,1:N_ADC/2));

RDM = 10*log10(IF_mat);

%% ��ͼ
% figure;
% mesh(IF_mat);
% xlabel('fast-time');
% ylabel('slow_time');
% zlabel('Amplitude');
% title('��Ƶ�����ź�ʱ��');
% 
% figure;
% plot(IF_mat(1,:));
% hold on
% plot(IF_Q_mat(1,:));
% xlabel('samples');
% ylabel('Amplitude');
% title('I-Q·�ź�');
% legend('I·','Q·');

figure;
distance_temp = (0:N_ADC/2 - 1) * Fs * c / N_ADC / 2 / slope;
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