% --2022/11/10-- %
% --�-- %
%% FMCW����
%% �ɵ���Ŀ�����
%% �ź���ʽ��ѡ *�ǵøĳɸ��ź���ʽ��������ķ���
%% �������ġ�A Peer-to-Peer Interference Analysis for Automotive Chirp Sequence Radars��
%% �½������FMCW�������ң� ��N_intchirp������һ������FMCW chirp�������뷢��chirp���ڱ�ֵ�� ����1ʱΪƽ�и��ţ�������Ŀ�ꡣ 

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
B_chirp = 150e6;          %chirp���� ��Դʦ��ĳ�������
T_chirp = 26e-6;          %chirp����ʱ��
slope = B_chirp/T_chirp;  %chirp��Ƶб��
Fs = 16e6;                %������ 
N_ADC = T_chirp *Fs;      %һ�����ڵĲ������� 26*16=416��
N_chirp = 128;             %chirp������
% t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ����������Tchirp�е�λ��
t = 0:1/Fs:(N_ADC-1)*1/Fs;  %ʱ������ ȷ����������Tchirp�е�λ��
lambda = c / fc ;         %����

IF_store = zeros(N_chirp,N_ADC);

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(4*slope);             %�״�̽�������� *ע����������/4������ԭ������д��/2 
f_IF_max = slope * 2 * distance_max/c /2   %��Ƶ����
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���
f_d = 2 * fc * V_target / c;

%% �����źŲ���

AT = 1;  %�����ź�����

%% CW���Ų�������
fint_CW = [77.2e9 77.3e9 77.5e9 77.55e9 77.6e9  78e9 78.5e9]; %CW����Ƶ�������������Ҳ� Ƶ��
R_int_CW = [1.5 20 10 7 8 45 50 ];                            %���ŵ�λ��
V_int_CW = [1.5 8 2 10 12 4 5 ];                              %�����ٶ�
N_int_CW = length(fint_CW);                                     %CW���Ÿ���

%% FMCW�����źŲ�������
R_int_FMCW = [80];
V_int_FMCW = [25];
N_int_FMCW = length(R_int_FMCW);

N_intchirp = 128;                             % һ��FMCW���ſ�Խ�ķ���chirp������
T_int_FMCW = N_intchirp*T_chirp;            % �����ź�����
N_int_division = ceil(N_chirp/N_intchirp);  %����ȡ�� ����FMCW��������

fint_FMCW = fc;                         % FMCW������ʼƵ��
B_int = B_chirp;                        % FMCW���Ŵ����뷢��chirp��ͬ

slope_int = B_int/T_int_FMCW;        % FMCW���ŵĵ�Ƶб��
N_FMCWint = floor(T_int_FMCW * Fs);  % һ�����ڸ����źŵ��� 416* N_int_chirp
% t_int = 0:1/Fs:T_int-1/Fs;      % �����źŲ�����λ��
t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;   % �����źŲ�����λ��

%% ����FMCW�����ź�
% S_int_FMCW = zeros(1,N_FMCWint);
% S_int_FMCW_targetsum = zeros(1,N_FMCWint);
S_int_FMCW_sum = zeros(1,N_FMCWint *N_int_division);

for int_FMCW = 0:1:N_int_division-1 %fmcw����������ѭ��
    
    S_int_FMCW_targetsum = zeros(1,N_FMCWint); % !**ע��**�������λ�� 
    
    %% FMCW����λ�ö�������Ϣcell
    for n_fmcw = 1:1:N_int_FMCW
        distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw).*(t_int+int_FMCW *T_int_FMCW);  %CW���ž������
        td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                       %����ʱ�Ӹ���   *ע������ʱ�Ӳ��������� ������ֱ�ӵ�����ն˵�
%         S_int_FMCW{n_fmcw,1} = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW�����ź�
        S_int_FMCW{n_fmcw,1} = AT*cos(2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2));      %FMCW�����ź� *ʵ�ź�
        S_int_FMCW_targetsum = S_int_FMCW_targetsum + S_int_FMCW{n_fmcw,1};
    end
    
%     S_int_FMCW = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW�����ź�
        S_int_FMCW_sum(int_FMCW*N_FMCWint+1:(int_FMCW+1)*N_FMCWint) = S_int_FMCW_targetsum;
end



%% ���䡢������������ ��Ӹ����ź�
for i = 0:1:N_chirp-1  %chirp����ѭ��
    
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
    %% CW����λ�ö�������Ϣcell
    for n_cw = 1:1:N_int_CW
        distant_int_CW{n_cw,1} = R_int_CW(n_cw) +V_int_CW(n_cw).*(t+i*T_chirp);    %CW���ž������
        td_int_CW{n_cw,1} = distant_int_CW{n_cw,1}/c;                              %����ʱ�Ӹ���   *ע������ʱ�Ӳ��������� ������ֱ�ӵ�����ն� û�о���Ŀ�귴��
    end
    
%     %% FMCW����λ�ö�������Ϣcell
%     for n_fmcw = 1:1:N_int_FMCW
%         distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw) .*(t+i*T_chirp); %FMCW���ž������
%         td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                %FMCW����ʱ�Ӹ���  *ע������ʱ�Ӳ���������
%     end

  %% �����ź�
%     St = AT*exp(1i*2*pi*(fc*(t-i*T_chirp)+slope/2*t.^2)); %�����ź�
    St = AT*cos(2*pi*(fc*t+slope/2*t.^2)); %�����ź�
    
    Sr_sum = zeros(1,N_ADC);                              %һ��chirp��Ӧ�ز��źŵ��ӿվ��� 1024��
    S_int_sum = zeros(1,N_ADC);                           %�����źſվ��� 1024��
    
    %% CW�����ź�
    for int_CW = 1:1:1 %���Ÿ���ѭ��
        
        S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %���Ҹ����ź�(��Ӧ��ͬλ�õĸ����ź�) 
        S_int_sum = S_int_sum + S_int_CW{int_CW,1};
    
    end
    
    
    %% ����û�и��ŵĻز��ź�
    for m = 1:1:1
        
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
        Sr{m,1} = A_target(m)*AT*cos(2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
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

    Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %�������
    
    %% ��Ƶ�õ���Ƶ
%     IF = St .*conj(Sr_sum);            %�޸���
%     IF = St .*conj(Sr_sum_withCW);     %CW����
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW����

    IF_I = St .* Sr_sum_withFMCW;   %FMCW����
    
    IF = IF_I + 1i*Sr_sum_withFMCW .*sin(2*pi*(fc*t+slope/2*t.^2));             %��ɽ����ź� 90�������� ��ȥ��Ƶ�ʲ���
    
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
    IF_I_mat(i+1,:) = IF_I;          %��ƵI·�ź�
    
    St_mat(i+1,:) = St;
    Sr_mat(i+1,:) = Sr_sum;
    
end

N_fastfft = 512;
N_slowfft =128;


%% ����Ƶ�ź�Ƶ��ͼ
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %ȡ��Ƶ����һ��1024��
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure;
plot(f_IF,IF_realfft);
title('��Ƶ�ź�fft')
xlabel('frequency(Hz)');


%% �ź�ʱƵ���� ��ʱ����Ҷ�任(����Щ����)
wlen = 416; %��������С
hop = wlen/4;
nfft =416;
win = blackman(wlen,'periodic');
% [S_if,F_if,T_if] = spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
figure;
spectrogram(IF_I_mat(1,:), win, wlen - hop, nfft, Fs);
title('�ܸ�����Ƶ�źŶ�ʱ����ҶʱƵͼ');



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
