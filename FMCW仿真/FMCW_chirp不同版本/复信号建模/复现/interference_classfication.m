%% 2022/11/18
%% ���������ڿ���������㷨�о�
%% ����1����Ƶ��ʱ�������ź� ����ſ�ʼʱ�������ʱ�� ��ÿ��chirp�������
%% ��Գ�ʱ����Ƶ�ʲ���һ�����Ҹ����źţ���Խ�Ĵ���ܴ� �������״����150M 
%% �������Ҹ���[77.04G, 77.08G, 77.12G]��Խ����0.08G=80M
%% 2022.11.21 ������������1.�²�ͬ�Ľ�����š������г��������FMCW����(�۲�б��) 2.��б��Ϊ�����Ը����źŽ��з��� 3.���³��� �ܸо����������⣿

%% �շ�����ʵ�źŽ�ģ

clear all;
close all;
clc;

%% Ŀ���������
R_target = [40 50 80];
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
N_chirp = 75;             %chirp������
% t = 0:1/Fs:T_chirp-1/Fs;  %ʱ������ ȷ����������Tchirp�е�λ��
t = 0:1/Fs:(N_ADC-1)*1/Fs;%ʱ������ ȷ����������Tchirp�е�λ��
lambda = c / fc ;         %����

BAAF = 8e6;               %��Ƶ�˲�������8M
IF_store = zeros(N_chirp,N_ADC);

%% ����FMCW�״����ܲ���
distance_max = c*Fs/(4*slope);             %�״�̽�������� *ע����������/4������ԭ������д��/2 
f_IF_max = slope * 2 * distance_max/c /2;  %��Ƶ����
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���
f_d = 2 * fc * V_target / c;

%% �����źŲ���

AT = 1;  %�����ź�����

%% CW���Ų�������
fint_CW = [ 77.04e9 77.08e9 77.12e9]; %CW����Ƶ�������������Ҳ� Ƶ�� 
R_int_CW = [ 0  0 0];                 %���ŵ�λ��
V_int_CW = [ 0  0 0];                 %�����ٶ�
N_int_CW = length(fint_CW);           %CW���Ÿ���


t_intCW_start = (fint_CW - fc - BAAF)/slope;   %CW���ſ�ʼӰ��ʱ��                     
n_intCW_start = floor(t_intCW_start*Fs);       %��ʼ��Ӧ�ĵ���
t_intCW_end = (fint_CW - fc + BAAF)/slope;     %���Ž���Ӱ��ʱ��
n_intCW_end = floor(t_intCW_end*Fs);           %������Ӧ�ĵ���
t_intCW_consist = 2*BAAF/slope;                %����Ӱ����ʱ��
n_intCW_consist = floor(t_intCW_consist*Fs);   %����Ӱ���ܵ���


%% FMCW�����źŲ�������
R_int_FMCW = [20  40 60];
V_int_FMCW = [5  2 5];
N_int_FMCW = length(R_int_FMCW);

N_FMCWintchirp = 1;                             % �˳����о�����ռ��������(ֻ����б��) �ʱ�ֵ��Ϊ1 һ��FMCW���ſ�Խ�ķ���chirp������
T_int_FMCW = N_FMCWintchirp*T_chirp;            % �����ź�����
N_int_division = ceil(N_chirp/N_FMCWintchirp);  % ����ȡ�� ����FMCW��������

fint_FMCW = [77.04e9   77.025e9 77.12e9];     % FMCW������ʼƵ��
B_int = [0   100e6 0];                        % FMCW���Ŵ����뷢��chirp��ͬ
slope_int = B_int./T_int_FMCW;                % FMCW���ŵĵ�Ƶб��
% N_FMCWint = floor(T_int_FMCW * Fs);  % һ�����ڸ����źŵ��� 416* N_int_chirp
% t_int = 0:1/Fs:T_int-1/Fs;         % �����źŲ�����λ��
% t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;   % �����źŲ�����λ��


t_intFMCW_start = (fint_FMCW -fc -BAAF)./(slope-slope_int);   % FMCW���ſ�ʼʱ��
n_intFMCW_start = floor(t_intFMCW_start*Fs);                  % FMCW���ſ�ʼ��Ӧ�ĵ���
t_intFMCW_end = (fint_FMCW - fc + BAAF)./(slope-slope_int);    % ����ʱ��
n_intFMCW_end = floor(t_intFMCW_end*Fs);                      % ������Ӧ�ĵ���
t_intFMCW_consist = 2*BAAF./(slope - slope_int);              % ���ų���ʱ��
n_intFMCW_consist = floor(t_intFMCW_consist*Fs);              % ���ų����ĵ���


%% ����FMCW�����ź�
% S_int_FMCW = zeros(1,N_FMCWint);
% S_int_FMCW_targetsum = zeros(1,N_FMCWint);
% S_int_FMCW_sum = zeros(1,N_FMCWint *N_int_division);

% for int_FMCW = 0:1:N_int_division-1 %fmcw����������ѭ��
%     
%     S_int_FMCW_targetsum = zeros(1,N_FMCWint); % !**ע��**�������λ�� 
%     
%     %% FMCW����λ�ö�������Ϣcell
%     for n_fmcw = 1:1:N_int_FMCW
%         distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw).*(t_int+int_FMCW *T_int_FMCW);  %CW���ž������
%         td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                       %����ʱ�Ӹ���   *ע������ʱ�Ӳ��������� ������ֱ�ӵ�����ն˵�
% %         S_int_FMCW{n_fmcw,1} = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW�����ź�
%         S_int_FMCW{n_fmcw,1} = AT*cos(2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2));      %FMCW�����ź� *ʵ�ź�
%         S_int_FMCW_targetsum = S_int_FMCW_targetsum + S_int_FMCW{n_fmcw,1};
%     end
%     
% %     S_int_FMCW = AT*exp(1i*2*pi*(fint_FMCW*(t_int-td_int_FMCW{n_fmcw,1})+slope_int/2*(t_int-td_int_FMCW{n_fmcw,1}).^2)); %FMCW�����ź�
%         S_int_FMCW_sum(int_FMCW*N_FMCWint+1:(int_FMCW+1)*N_FMCWint) = S_int_FMCW_targetsum;
% end


%% ���䡢������������ ��Ӹ����ź�
for i = 0:1:N_chirp-1  %chirp����ѭ��
    
    %% Ŀ��λ�ö�������Ϣcell
    for k = 1:1:N_target
        distant{k,1} = R_target(k) + V_target(k).*(t+i*T_chirp); %�������
        t_d{k,1} = 2*distant{k,1}/c;                             %ʱ�Ӹ���
    end
    
    %% CW����λ�ö�������Ϣcell
    for n_cw = 1:1:N_int_CW
        distant_int_CW{n_cw,1} = R_int_CW(n_cw) + V_int_CW(n_cw).*(t+i*T_chirp);    %CW���ž������
        td_int_CW{n_cw,1} = distant_int_CW{n_cw,1}/c;                               %����ʱ�Ӹ���   *ע������ʱ�Ӳ��������� ������ֱ�ӵ�����ն� û�о���Ŀ�귴��
    end
    
    %% FMCW����λ�ö�������Ϣcell
    for n_fmcw = 1:1:N_int_FMCW
        distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) + V_int_FMCW(n_fmcw) .*(t+i*T_chirp); %FMCW���ž������
        td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                 %FMCW����ʱ�Ӹ���  *ע������ʱ�Ӳ���������
    end

  %% �����ź�
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2)); %�����ź�
    St = AT*cos(2*pi*(fc*t+slope/2*t.^2)); %�����ź�
    
    Sr_sum = zeros(1,N_ADC);                              %һ��chirp��Ӧ�ز��źŵ��ӿվ��� 1024��
    S_int_sum = zeros(1,N_ADC);                           %�����źſվ��� 1024��
    
%     %% CW�����ź�
%     for int_CW = 1:1:3 %���Ÿ���ѭ��
%         
% %         S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %���Ҹ����ź�(��Ӧ��ͬλ�õĸ����ź�)
%         S_int_CW{int_CW,1} = 10*cos(2*pi*fint_CW(int_CW)*(t-td_int_CW{int_CW,1})); %CW���� *ʵ�ź�
%         S_int_sum = S_int_sum + S_int_CW{int_CW,1};
%         
%     end
    
    
    %% ����û�и��ŵĻز��ź�
    for m = 1:1:1
        
%         Sr{m,1} = A_target(m)*AT*exp(1i*2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
        Sr{m,1} = A_target(m)*AT*cos(2*pi*(fc*(t-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
        Sr_sum = Sr_sum + Sr{m,1};
        
    end
    
    
     %% CW�����ź�
%     Sr_sum_withCW = zeros(1,N_ADC);
%     Sr_sum_addCW = zeros(1,N_ADC);
%     for int_CW = 1:1:N_int_CW %���Ÿ���ѭ��
%         
% %         S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %���Ҹ����ź�(��Ӧ��ͬλ�õĸ����ź�)
%         S_int_CW{int_CW,1} = 20*cos(2*pi*fint_CW(int_CW)*(t-td_int_CW{int_CW,1})); %CW���� *ʵ�ź�
% %         S_int_sum = S_int_sum + S_int_CW{int_CW,1};
%         Sr_sum_addCW = [Sr_sum(1:n_intCW_start(int_CW)-1),Sr_sum(n_intCW_start(int_CW):n_intCW_end(int_CW))+S_int_CW{int_CW,1}(n_intCW_start(int_CW):n_intCW_end(int_CW)),Sr_sum(n_intCW_end(int_CW)+1:N_ADC)];
%         Sr_sum_withCW = Sr_sum_withCW + Sr_sum_addCW;
%     end
    
    %% FMCW�����ź�
    Sr_sum_withFMCW = zeros(1,N_ADC);
    Sr_sum_addFMCW = zeros(1,N_ADC);
    for int_FMCW = 1:1:N_int_FMCW %���Ÿ���ѭ��
        
        S_int_FMCW{int_FMCW,1} = 20*cos(2*pi*(fint_FMCW(int_FMCW)*(t-td_int_FMCW{int_FMCW,1}) + slope_int(int_FMCW)/2*(t-td_int_FMCW{int_FMCW,1}).^2)); %FMCW���� *ʵ�ź�
        Sr_sum_addFMCW =[Sr_sum(1:n_intFMCW_start(int_FMCW)-1),Sr_sum(n_intFMCW_start(int_FMCW):n_intFMCW_end(int_FMCW))+S_int_FMCW{int_FMCW,1}(n_intFMCW_start(int_FMCW):n_intFMCW_end(int_FMCW)),Sr_sum(n_intFMCW_end(int_FMCW)+1:N_ADC)]; 
        Sr_sum_withFMCW = Sr_sum_withFMCW + Sr_sum_addFMCW;
        
    end
    
    %% ----�������ģ��-----
    %% CW����
%     if i <= 80
%         
% %         Sr_sum_withCW = Sr_sum + S_int_sum;
%         Sr_sum_withCW = [Sr_sum(1:n_int_start-1),Sr_sum(n_int_start:n_int_end)+S_int_sum(n_int_start:n_int_end),Sr_sum(n_int_end+1:N_ADC)];
%         
%     else
%         Sr_sum_withCW = Sr_sum; 
%     end
    
%     %% FMCW����
% 
%     Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %�������
    
    %% ��Ƶ�õ���Ƶ
%     IF = St .*conj(Sr_sum);            %�޸���
%     IF = St .*conj(Sr_sum_withCW);     %CW����
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW����
    
    %% ��Ƶ�õ���Ƶ
%     IF = St .* Sr_sum;
%     IF = IF + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));           %�޸���
    
    IF_I = St .* Sr_sum_withFMCW; 
    IF = IF_I + 1i*Sr_sum_withFMCW .*sin(2*pi*(fc*t+slope/2*t.^2));  %FMCW����
    IF_I_mat(i+1,:) = IF_I;                                          %* ���ջ���Ƶ�ź�I���� ����ʱƵ����
    
%     IF_I = St .*Sr_sum_withCW;                                       %CW���Ž��ջ���Ƶ�õ�I����
%     IF = IF_I + 1i*Sr_sum_withCW .*sin(2*pi*(fc*t+slope/2*t.^2));    %������Ҹ���
%     IF_I_mat(i+1,:) = IF_I;                                          %* ���ջ���Ƶ�ź�I���� ����ʱƵ����
    
%     IF = St .*conj(Sr_sum);                                        %���źŽ�ģ���  *�޸���                                
    
    IF_wonoise = St .*Sr_sum + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));  %* ��������Ƶ�źŽ���ʱƵ������
    IF_wonoise_mat(i+1,:) = IF_wonoise;                                    
    IF_wonoise_I_mat(i+1,:) = St .*Sr_sum;                                 %* ��������Ƶ�ź�I���� ����ʱƵ������

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

% N_fastfft = 512;
% N_slowfft =128;
% 
% IF_mat_trans = IF_mat';
% % ����fft
% sig_fft = fft(IF_mat_trans,N_ADC)./N_ADC;
% sig_fft = abs(sig_fft);
% sig_fft = sig_fft(1:(N_ADC),:);
% %������fft
% sig_fft2 = fft2(IF_mat_trans,N_ADC,N_chirp);
% sig_fft2 = sig_fft2(1:N_ADC,1:N_chirp);
% sig_fft2 = fftshift(sig_fft2);
% RDM = abs(sig_fft2);
% RDM = 10*log10(RDM) ;
% doppler_axis = linspace(-100,100,N_chirp);
% range_axis = linspace(-200,200,N_ADC)*((N_ADC)/400);
% %��ͼ
% mesh(doppler_axis,range_axis,RDM);
% xlabel('�ٶ�'); ylabel('����'); zlabel('���ȣ�dB��');
% title('�ٶ�άFFT �����������');


%% ����Ƶ�ź�Ƶ��ͼ
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %ȡ��Ƶ����һ��1024��
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure(1);
plot(f_IF,IF_realfft);
title('��Ƶ�ź�fft')
xlabel('frequency(Hz)');


%% �ź�ʱƵ���� ��ʱ����Ҷ�任(����Щ����)
tfr_St = tfrstft(conj(IF_mat(1,1:N_ADC)'));
figure(2);
% imagesc(abs(tfr_St(N_ADC/2-1:N_ADC,:)));
imagesc(abs(tfr_St));
title('�ܸ�����Ƶ�źŶ�ʱ����ҶʱƵͼ');

tfr_IF =  tfrstft(IF_wonoise_I_mat(1,1:N_ADC)');
figure(3);
imagesc(abs(tfr_IF(N_ADC/2-1:N_ADC,:)));
title('��Ƶ�źŶ�ʱ����ҶʱƵͼ');

figure(4);
[B_IF,t_IF,f_IF_I]=tfrstft(IF_I_mat(1,1:N_ADC)');
imagesc(t_IF,f_IF_I(N_ADC/2:N_ADC-1)*Fs,abs(B_IF(N_ADC/2:N_ADC-1,:)));


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