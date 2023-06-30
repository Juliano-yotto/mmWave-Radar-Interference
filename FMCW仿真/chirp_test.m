%% ������
% c = 3e8; %����
% fc = 76.5e9;  %�����ź���Ƶ ����Ƶ��
% bw = 500e6;  %�����źŴ���
% Tr = 10e-6;  %ɨƵʱ�� Ҳ��������
% N = 256;  %������
% Fs = 25.6e6; %������
% M = 256;  %chirp����Ŀ
% k = bw/Tr; %chirpб��
% index = 1:1:N;  %����������
% IF_mat = zeros(M,N);  %�洢������������Ƶ�ź�

%% ��������
c = 3e8;
fc = 77e9; %chirp��ʼƵ��
bw = 1.6e9; %chirp����
Tr = 62.48*10^(-6) %chirp����ʱ��
N = 1024; %chirp�����ڲ�������
M = 128; %chirp����
Fs = N/Tr; %������ 
k = bw/Tr; %chirp��Ƶб��
index = 1:1:N;
IF_mat = zeros(M,N);

%% �����źŲ���
AT = 1;  %�����ź�����
t = 0:1/Fs:Tr-1/Fs;  %ʱ������ ȷ��256������һ��Tr�е�ÿ��ʱ��
%t = t - Tr/2; %��fc��Ϊ����Ƶ��

%% �ز��źŲ���
distance = 80;  %Ŀ������״�50m�ľ���
t_d = 2 * distance / c;  %Ŀ������״���ӳ�
velocity = 10;  %Ŀ����״������ٶ�Ϊ30m/s
f_d = 2 * (fc - bw/2) * velocity / c;  %������Ƶ��
AR = 0.8;  %�ز��ź�˥���ı���ֵ

%% ��������
for i = 1:1:M  %chirp��ѭ��
%     St = AT*exp((1i*2*pi)*(fc*(t-i*Tr)+k/2*t.^2)); %�����ź�
    St = AT*cos(2*pi*(fc*(t-i*Tr)+k/2*t.^2));
%     Sr = AR*AT*exp((1i*2*pi)*((fc-f_d)*(t-t_d-i*Tr)+k/2*(t-t_d).^2));  %�ز��ź�
    Sr = AR*AT*cos(2*pi*((fc-f_d)*(t-t_d-i*Tr)+k/2*(t-t_d).^2));  %�ز��ź�
    %% ��ز��źŵĹ���
%     Sr = conj(Sr);  %��ز��źŵĹ���
    
    %% ����Ƶ�ź�
    IF = St .* Sr + 1i*Sr.*sin(2*pi*(fc*(t-i*Tr)+k/2*t.^2));  %����Ƶ�ź�
    SNR = 1e15;  %�����
    IF_with_Noise = awgn(IF,SNR,'measured');  %����Ƶ�źżӸ�˹�������������������ʱ��Ҫ���������Ĳ���
    IF_mat(i,:) = IF_with_Noise;  %��������������Ƶ�źű���
    IF_mat_1(i,:) = IF_with_Noise;
    IF_mat1(i,:) = IF;%����Ƶ�ʷ��� ���ù�
end
save('Ego_vehicle.mat', 'IF_mat');  %�������ݵı���

%% ��������
IF_mat = cell2mat(struct2cell(load('Ego_vehicle.mat','IF_mat')));
%% ������
% c = 3e8; %����
% fc = 76.5e9;  %�����ź���Ƶ ����Ƶ��
% bw = 500e6;  %�����źŴ���
% Tr = 10e-6;  %ɨƵʱ�� Ҳ��������
% N = 256;  %������
% Fs = 25.6e6; %������
% M = 256;  %chirp����Ŀ
% k = bw/Tr; %chirpб��

lambda = c / (fc - bw/2);
point = 1:1:N;  %����������

%% ��Ƶ�ź�ʱƵ����
tfr_IF =  tfrstft(IF_mat_1(1,1:1024)');
figure;
imagesc(abs(tfr_IF));
title('��ʱ����ҶʱƵͼ');

% %% Ƶ�ʷ���
% figure;
% IF_fft = fft(IF_mat1(1,:),1024);
% plot(IF_fft);

%% ���ɴ�
range_win = hamming(N);  %����range��
doppler_win = hamming(M);  %����doppler��
%% range fft
for i = 1:1:M
    temp = IF_mat(i,:) .* range_win';
    temp_fft = fft(temp,N);
    IF_mat(i,:) = temp_fft;
end
%% doppler fft
for j = 1:1:N
    temp = IF_mat(:,j) .* doppler_win;
    temp_fft = fftshift(fft(temp,M));
    IF_mat(:,j) = temp_fft;
end
%% ��ͼ
figure;
distance_temp = (0:N - 1) * Fs * c / N / 2 / k;
speed_temp = (-M / 2:M / 2 - 1) * lambda / Tr / M / 2;
[X,Y] = meshgrid(distance_temp,speed_temp);
mesh(X,Y,(abs(IF_mat)));
xlabel('����(m)');
ylabel('�ٶ�(m/s)');
zlabel('�źŷ�ֵ');
title('2άFFT������ά��ͼ');

figure;
speed_temp = -speed_temp;
imagesc(distance_temp,speed_temp,abs(IF_mat));
title('����-��������ͼ');
xlabel('����(m)');
ylabel('�ٶ�(m/s)');