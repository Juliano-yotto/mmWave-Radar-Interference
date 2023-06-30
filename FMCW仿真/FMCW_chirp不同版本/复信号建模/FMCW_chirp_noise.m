% --2022/11/10-- %
% --�-- %
%% FMCW����
%% �ɵ���Ŀ�����
%% ���ø��ź���ʽ
%% ����˲�ͬλ�á���ͬ�ٶȡ���ͬƵ�ʵ����Ҹ��� �ٶȵ�Ԫ���ָ������� ���ٶȵ�Ԫ�ϵļ�����Ӱ��

clear all;
close all;
clc;

%% Ŀ���������
R_target = [45 50 80];
V_target = [10 6 14];
A_target = [1 0.1 0.5];      %��ͬĿ���Ӧ�Ľ���˥��
N_target = length(R_target); %Ŀ�����

%% FMCW�״�����������
c = 3e8;
fc = 77e9;                %chirp��ʼƵ��
B_chirp = 1.6e9;          %chirp���� ��Դʦ��ĳ�������
T_chirp = 62.48e-6;       %chirp����ʱ��
slope = B_chirp/T_chirp;  %chirp��Ƶб��
N_ADC = 1024;             %chirp�����ڲ�������
N_chirp = 128;            %chirp����
Fs = N_ADC/T_chirp;       %������ 
% t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ����������Tchirp�е�λ��
t = 0:1/Fs:(N_ADC-1)*1/Fs;  %ʱ������ ȷ����������Tchirp�е�λ��
lambda = c / fc ;         %����

IF_store = zeros(N_chirp,N_ADC);

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(2*slope);             %�״�̽��������
f_IF_max = slope * 2 * distance_max/c /2   %��Ƶ����
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���
f_d = 2 * fc * V_target / c;

%% �����źŲ���

AT = 1;  %�����ź�����

%% CW���Ų�������
fint1 = [77.2e9 77.3e9 77.5e9 77.55e9 77.6e9  78e9 78.5e9]; %CW����Ƶ�������������Ҳ� Ƶ��
R_int1 = [1.5 20 10 7 8 45 50 ];                            %���ŵ�λ��
V_int1 = [1.5 8 2 10 12 4 5 ];                              %�����ٶ�
N_int1 = length(fint1);                                     %CW���Ÿ���

%% FMCW�����źŲ�������
N_int_chirp = 7;                          % һ��FMCW���ſ�Խ�ķ���chirp������
T_int = N_int_chirp*T_chirp;                 % �����ź�����
N_int_division = ceil(N_chirp/N_int_chirp);  %����ȡ�� ����FMCW��������

fint2 = fc;                         % FMCW������ʼƵ��
B_int = B_chirp;                    % FMCW���Ŵ����뷢��chirp��ͬ

slope_int = B_int/T_int;        % FMCW���ŵĵ�Ƶб��
N_FMCWint = floor(T_int * Fs);  % һ�����ڸ����źŵ��� 1024*N_int_chirp
% t_int = 0:1/Fs:T_int-1/Fs;      % �����źŲ�����λ��
t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;% �����źŲ�����λ��

%% ����FMCW�����ź�
S_int_FMCW = zeros(1,N_FMCWint);
S_int_FMCW_sum = zeros(1,N_FMCWint *N_int_division);

for int_FMCW = 0:1:N_int_division-1
    
    S_int_FMCW = AT*exp(1i*2*pi*(fint2*(t_int-int_FMCW*T_int)+slope_int/2*t_int.^2)); %FMCW�����ź�
    S_int_FMCW_sum(int_FMCW*N_FMCWint+1:(int_FMCW+1)*N_FMCWint) = S_int_FMCW ;
    
end



%% ���䡢������������ ��Ӹ����ź�
for i = 0:1:N_chirp-1  %chirp����ѭ��
    
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
    %% ����λ�ö�������Ϣcell
    for n = 1:1:N_int1
        distant_int{n,1} = R_int1(n) +V_int1(n).*(t+i*T_chirp);    %���ž������
        td_int{n,1} = 2*distant_int{n,1}/c;                        %����ʱ�Ӹ���
    end
    
    St = AT*exp(1i*2*pi*(fc*(t-i*T_chirp)+slope/2*t.^2)); %�����ź�
    
    Sr_sum = zeros(1,N_ADC);                              %һ��chirp��Ӧ�ز��źŵ��ӿվ��� 1024��
    S_int_sum = zeros(1,N_ADC);                           %�����źſվ��� 1024��
    
    %% CW�����ź�
    for int_CW = 1:1:1 %���Ÿ���ѭ��
        
        S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint1(int_CW)*(t-i*T_chirp-td_int{int_CW,1})); %���Ҹ����ź�(��Ӧ��ͬλ�õĸ����ź�) 
        S_int_sum = S_int_sum + S_int_CW{int_CW,1};
    
    end
    
    
    %% ����û�и��ŵĻز��ź�
    for m = 1:1:1
        
        Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1}-i*T_chirp)+slope/2*(t-t_d{m,1}).^2)) + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %����FMCW����
        Sr_sum = Sr_sum + Sr{m,1};
        
    end
    
    %% ----�������ģ��-----
    %% CW����
    if i <= 50
        
        Sr_sum_withCW = Sr_sum + S_int_sum;
        
    else
        Sr_sum_withCW = Sr_sum; 
    end
    
    %% FMCW����
%     if i <= N_int_chirp-1 %FMCW����Ӱ���chirp����
%         
%         Sr_sum_withFMCW = Sr_sum + S_int_FMCW(i*N_ADC+1:(i+1)*N_ADC);  %��Ӱ�������ڼ���FMCW�����ź�
%         
%     else
%         Sr_sum_withFMCW = Sr_sum; %���ڸ�����������û�и���
%     end 
    Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %�������
    
    %% ��Ƶ�õ���Ƶ
    IF = St .*conj(Sr_sum);            %�޸���
%     IF = St .*conj(Sr_sum_withCW);     %CW����
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW����
    

    %% ���õ�ͨ�˲�
    %% ��ư�����˹��ͨ�˲��� ģ��or���֣�
%     omg_p = 8e6;
%     omg_s = 8.15e6;
%     w_p = 2*pi*omg_p;
%     w_s = 2*pi*omg_s;
%     Rp = 1;
%     As = 30
%     [N,wc] = buttord(w_p,w_s,Rp,As,'s'); %'s'����ģ���˲���
%     [B,A] = butter(N,wc);
    
%     %% �˲�����Ƶ����
%     f_filter = 0:10:1.65e7;
%     w = 2*pi * f_filter;
%     Hk = freqs(B,A,w);
%     figure;
%     plot(f_filter/1e6,abs(Hk));
    
    %% �˲�
%     IF = filter(B,A,IF);
   

    IF_mat(i+1,:) = IF;              %����ʱƵ����
    IF_mat_withnoise(i+1,:) = IF;    %��������������Ƶ�źű���
    St_mat(i+1,:) = St;
    Sr_mat(i+1,:) = Sr_sum;
    
end


%% ����Ƶ�ź�Ƶ��ͼ
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %ȡ��Ƶ����һ��1024��
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure;
plot(f_IF,IF_realfft);
title('��Ƶ�ź�fft')
xlabel('frequency(Hz)');


% %% �ź�ʱƵ���� ��ʱ����Ҷ�任(����Щ����)
% tfr_St = tfrstft(St_mat(1,1:1024)');
% figure;
% imagesc(abs(tfr_St));
% title('�����źŶ�ʱ����ҶʱƵͼ');
% 
% tfr_IF =  tfrstft(IF_mat(1,1:1024)');
% figure;
% imagesc(abs(tfr_IF));
% title('��Ƶ�źŶ�ʱ����ҶʱƵͼ');



%% �Ӻ���������
range_win = hanning(N_ADC);      %����range��
doppler_win = hanning(N_chirp);  %����doppler��

%% range fft
for i = 1:1:N_chirp
    
    temp = IF_mat_withnoise(i,:); %.* range_win'; %��ÿһ�н���fft
    temp_fft = fft(temp,N_ADC)./N_ADC;

    IF_mat_withnoise(i,:) = temp_fft;           %����IF_mat_withnoise
    
end

%% doppler fft
for j = 1:1:N_ADC
    
    temp = IF_mat_withnoise(:,j); %.* doppler_win; %�ٶ�ÿһ�н���fft
    temp_fft = fftshift(fft(temp,N_chirp))./N_chirp;      
    IF_mat_withnoise(:,j) = temp_fft;            %����IF_mat_withnoise
    
end

IF_mat_withnoise = abs(IF_mat_withnoise);
RDM = 10*log10(IF_mat_withnoise); %����dB

Gs_I = 10*log10(N_ADC)%���������

%% ��RDͼ
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

% % % %% ****************
% % % sig_fft2_i = fft2(IF_mat,N_chirp,N_ADC);
% % % sig_fft2_i = sig_fft2_i(1:N_chirp,1:N_ADC/2);
% % % sig_fft2_i = fftshift (sig_fft2_i);
% % % sig_fft2_i = fftshift(sig_fft2_i);
% % % RDM_i = abs(sig_fft2_i);
% % % RDM_i = 10*log10(RDM_i);
% % % 
% % % doppler_axis_i = linspace(-100,100,N_chirp);
% % % range_axis_i = linspace(-200,200,N_ADC/2)*((N_ADC/2)/400);
% % % figure;
% % % 
% % % mesh(range_axis_i,doppler_axis_i,RDM_i);