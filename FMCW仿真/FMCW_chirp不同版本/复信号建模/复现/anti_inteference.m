%% 2022/11/18
%% ���������ڿ���������㷨�о�
%% ����1����Ƶ��ʱ�������ź� ����ſ�ʼʱ�������ʱ�� ��ÿ��chirp�������
%% ��Գ�ʱ����Ƶ�ʲ���һ�����Ҹ����źţ���Խ�Ĵ���ܴ� �������״����150M 
%% �������Ҹ���[77.04G, 77.08G, 77.12G]��Խ����0.08G=80M

%% �շ�����ʵ�źŽ�ģ


clear all;
close all;
clc;

%% Ŀ���������
R_target = [40 50 80];
V_target = [1.5 6 14];
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
f_IF_max = slope * 2 * distance_max/c /2   %��Ƶ����
distance_res = c/(2*B_chirp);              %����ֱ���
velocity_max = lambda/(4*T_chirp);         %�״��������ٶ�
velocity_res = lambda/(2*N_chirp*T_chirp); %�ٶȷֱ���
f_d = 2 * fc * V_target / c;

%% �����źŲ���

AT = 1;  %�����ź�����

%% CW���Ų�������
fint_CW = [ 77.040e9 77.08e9]; %CW����Ƶ�������������Ҳ� Ƶ�� 
R_int_CW = [ 0 0 0];                            %���ŵ�λ��
V_int_CW = [ 0 0 0];                              %�����ٶ�
N_int_CW = length(fint_CW);                                   %CW���Ÿ���


t_int_start = (fint_CW - fc - BAAF)/slope;   %���ſ�ʼӰ��ʱ��                     
n_int_start = floor(t_int_start*Fs);       %��ʼ��Ӧ�ĵ���
t_int_end = (fint_CW - fc + BAAF)/slope;     %���Ž���Ӱ��ʱ��
n_int_end = floor(t_int_end*Fs);           %������Ӧ�ĵ���
t_int_consist = 2*BAAF/slope;              %����Ӱ����ʱ��
n_int_consist = floor(t_int_consist*Fs);       %����Ӱ���ܵ���


% %% FMCW�����źŲ�������
% R_int_FMCW = [80];
% V_int_FMCW = [25];
% N_int_FMCW = length(R_int_FMCW);
% 
% N_intchirp = 5;                             % һ��FMCW���ſ�Խ�ķ���chirp������
% T_int_FMCW = N_intchirp*T_chirp;            % �����ź�����
% N_int_division = ceil(N_chirp/N_intchirp);  %����ȡ�� ����FMCW��������
% 
% fint_FMCW = fc;                         % FMCW������ʼƵ��
% B_int = B_chirp;                        % FMCW���Ŵ����뷢��chirp��ͬ
% 
% slope_int = B_int/T_int_FMCW;        % FMCW���ŵĵ�Ƶб��
% N_FMCWint = floor(T_int_FMCW * Fs);  % һ�����ڸ����źŵ��� 416* N_int_chirp
% % t_int = 0:1/Fs:T_int-1/Fs;      % �����źŲ�����λ��
% t_int = 0:1/Fs:(N_FMCWint-1)*1/Fs;   % �����źŲ�����λ��
% 
% %% ����FMCW�����ź�
% % S_int_FMCW = zeros(1,N_FMCWint);
% % S_int_FMCW_targetsum = zeros(1,N_FMCWint);
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
        distant_int_CW{n_cw,1} = R_int_CW(n_cw) +V_int_CW(n_cw).*(t+i*T_chirp);    %CW���ž������
        td_int_CW{n_cw,1} = distant_int_CW{n_cw,1}/c;                              %����ʱ�Ӹ���   *ע������ʱ�Ӳ��������� ������ֱ�ӵ�����ն� û�о���Ŀ�귴��
    end
    
%     %% FMCW����λ�ö�������Ϣcell
%     for n_fmcw = 1:1:N_int_FMCW
%         distant_int_FMCW{n_fmcw,1} = R_int_FMCW(n_fmcw) +V_int_FMCW(n_fmcw) .*(t+i*T_chirp); %FMCW���ž������
%         td_int_FMCW{n_fmcw,1} = distant_int_FMCW{n_fmcw,1}/c;                                %FMCW����ʱ�Ӹ���  *ע������ʱ�Ӳ���������
%     end

  %% �����ź�
%     St = AT*exp(1i*2*pi*(fc*t+slope/2*t.^2)); %�����ź�
    St = AT*cos(2*pi*(fc*(t-i*T_chirp)+slope/2*t.^2)); %�����ź�
    
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
        Sr{m,1} = A_target(m)*AT*cos(2*pi*(fc*(t-i*T_chirp-t_d{m,1})+slope/2*(t-t_d{m,1}).^2));  %Ŀ��ز��ź�
        Sr_sum = Sr_sum + Sr{m,1};
        
    end
    
    Sr_sum_withCW = zeros(1,N_ADC);
    Sr_sum_addCW = zeros(1,N_ADC);
    
     %% CW�����ź�
    for int_CW = 1:1:N_int_CW %���Ÿ���ѭ��
        
%         S_int_CW{int_CW,1} = AT*exp(1i*2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %���Ҹ����ź�(��Ӧ��ͬλ�õĸ����ź�)
        S_int_CW{int_CW,1} = 20*cos(2*pi*fint_CW(int_CW)*(t-i*T_chirp-td_int_CW{int_CW,1})); %CW���� *ʵ�ź�
%         S_int_sum = S_int_sum + S_int_CW{int_CW,1};
        Sr_sum_addCW = [Sr_sum(1:n_int_start(int_CW)-1),Sr_sum(n_int_start(int_CW):n_int_end(int_CW))+S_int_CW{int_CW,1}(n_int_start(int_CW):n_int_end(int_CW)),Sr_sum(n_int_end(int_CW)+1:N_ADC)];
        Sr_sum_withCW = Sr_sum_withCW + Sr_sum_addCW;
        
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
    
    %% FMCW����

%     Sr_sum_withFMCW = Sr_sum + S_int_FMCW_sum(i*N_ADC+1:(i+1)*N_ADC); %�������
    
    %% ��Ƶ�õ���Ƶ
%     IF = St .*conj(Sr_sum);            %�޸���
%     IF = St .*conj(Sr_sum_withCW);     %CW����
%     IF = St .*conj(Sr_sum_withFMCW);   %FMCW����
    
    %% ��Ƶ�õ���Ƶ
%     IF = St .* Sr_sum;
%     IF = IF + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));           %�޸���
    
%     IF = St .* Sr_sum_withFMCW; 
%     IF = IF + 1i*Sr_sum_withFMCW .*sin(2*pi*(fc*t+slope/2*t.^2));  %FMCW����
    
    IF_I = St .*Sr_sum_withCW;                                       %���ջ���Ƶ�õ�I����
    IF = IF_I + 1i*Sr_sum_withCW .*sin(2*pi*(fc*t+slope/2*t.^2));    %������Ҹ���
%     IF = hilbert(IF_I);
    IF_I_mat(i+1,:) = IF_I;                                          %* ���ջ���Ƶ�ź�I���� ����ʱƵ����
    
%     IF = St .*conj(Sr_sum);                                        %���źŽ�ģ���  *�޸���                                
    
    IF_wonoise = St .*Sr_sum + 1i*Sr_sum .*sin(2*pi*(fc*t+slope/2*t.^2));  %* ��������Ƶ�źŽ���ʱƵ������
    IF_wonoise_mat(i+1,:) = IF_wonoise;                                    
    IF_wonoise_I_mat(i+1,:) = St .*Sr_sum;                                 %* ��������Ƶ�ź�I���� ����ʱƵ������

   

    IF_mat(i+1,:) = IF;              %����ʱƵ����
    IF_mat_withnoise(i+1,:) = IF;    %��������������Ƶ�źű���
    St_mat(i+1,:) = St;
    Sr_mat(i+1,:) = Sr_sum;
    
end


%% ����Ƶ�ź�ʱ��ͼ
figure(1);
IF_realtime = IF_I_mat(1,:); %ȡ��һ��chirp����Ƶ�����ź�
plot(IF_realtime)
title('һ��chirp��Ӧ����Ƶ�ź�ʱ�����')
xlabel('��������N');

%% ����Ƶ�ź�Ƶ��ͼ
IF_realfft =abs(fft(IF_mat(1,:),N_ADC)./N_ADC); %ȡ��Ƶ����һ��1024��
f_IF = 0:N_ADC-1;
f_IF = f_IF *Fs/N_ADC;
figure(2);
plot(f_IF,IF_realfft);
title('��Ƶ�ź�fft')
xlabel('frequency(Hz)');


%% �ź�ʱƵ���� ��ʱ����Ҷ�任STFT(�ݶ�)
[tfr_IF,t_IF , f_IF_I]=  tfrstft(conj(IF_wonoise_mat(1,1:N_ADC)'));
figure(3);
imagesc(t_IF,f_IF_I(1:length(f_IF_I)/2)*Fs,abs(tfr_IF(1:N_ADC/2,:)));
title('��������Ƶ�źŶ�ʱ����ҶʱƵͼ');

[tfr_St,t_st,f_st] = tfrstft(conj(IF_mat(1,1:N_ADC)'));
figure(4);
imagesc(t_st,f_st(1:length(f_st)/2)*Fs,abs(tfr_St(1:N_ADC/2,:)));
title('�ܸ�����Ƶ�źŶ�ʱ����ҶʱƵͼ');

figure(5);
tfr_IF = tfrstft(conj(IF_mat(1,1:N_ADC)'));
imagesc(abs(tfr_IF));


% figure;
% [B_IF,t_IF,f_IF_I]=tfrstft(conj(IF_mat(1,1:N_ADC)'));             %�������� ��ֱ��ת������ķ���Ҳ����
% % imagesc(t_IF,f_IF_I(N_ADC/2:N_ADC-1)*Fs,abs(B_IF(N_ADC/2:N_ADC-1,:)));
% imagesc(t_IF/Fs,f_IF_I(1:N_ADC/2)*Fs,abs(B_IF));


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