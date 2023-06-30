% ������
% һ��Ŀ��
% ��Ϊ���Ҹ���Ӱ�����

%% ���̽�����Լ51.89m��c*fs/2s��������ֱ���Լ0.09m(c/2B)
%% �������ٶ�   ��
%% ��ֹ��Ŀ��
%% 139�� ΪʲôҪ�Ӵ���
% fs=200
% t=cell(2,1);
% t(1,1)={1:1/fs:2};
% t(2,1)={1:1/fs:3};

clear all;
close all;
clc;

c = 3*10^8; %����
R1 = 1.3; %�״���Ŀ������1�ľ��룬��λΪm


td_1 = 2*R1/c;   %Ŀ������1��Ӧ�Ļز��ź�ʱ��

tChirp=62.48*10^(-6);   %chirp�ĳ���ʱ��
k=26.023*10^12; %б��
B=k*tChirp; %����  Լ1.6GHz
fc=77*10^9; %�ز�Ƶ�� 

N = 512;    %��Ƶ�źŲ�����
fs = 9000*10^3;    %����Ƶ��9MHz
tADC=(N-1)/fs;  %��������ʱ�� 56.7us
t_start=4*10^(-6);  %������ʼʱ��
t=t_start:1/fs:(t_start+tADC);

% ����ز��źŵĹ���
Pt_dBm = 12;    %�״�ķ��书�ʣ���λΪdBm
Pt_W = 0.01585; %�״�ķ��书�ʣ���λΪW
% Pt_dBm = 18;    %�״�ķ��书�ʣ���λΪdBm
% Pt_W = 0.06309573444801932; %�״�ķ��书�ʣ���λΪW
Pt_dB = 10*log10(Pt_W);
Gt_dB = 10.5;  %�״�ķ����������棬��λΪdBi��dB��
Gt = 11.22; %�״�ķ�����������
Gr_dB = 10.5; %�״�������ߵ�����,��λΪdB
Gr = 11.22;  %�״�������ߵ�����

lambda = c/fc;  %��������λΪm
RCS_1 = 0.2;   %Ŀ������1��Ӧ���״�ɢ������������λΪm^2


Pr_W_1 = Pt_W*Gt*Gr*lambda^2*RCS_1/((4*pi)^3*R1^4);


R_r = 50;
Ar_1 = sqrt(Pr_W_1 * R_r * 2);


Ar_1_Amplified = 31.62*Ar_1;


VRef = 1.8;
ADCCode_1 = Ar_1_Amplified*(2^16)/VRef; %Ŀ������1��Ӧ��ADC����������


%��Ϊ���Ҹ����źŲ���
% PI_W = 5*10^(-7); %���յ�����Ϊ���Ҹ����źŹ���

% ���������Ϣ������յ��ĸ����źŹ���
P_itx = 10; %��Ϊ���Ҹ��ŷ��书�ʣ���λΪdBm
txAntGain = 20; %�����������棬��λΪdB
rxAntGain = 10.5; %�״�����������棬��λΪdB
Ri = 1.3; %�����������״�֮��ľ���
PI_Calculation_dBm = P_itx + txAntGain + rxAntGain - 10*log10((4*pi*Ri/lambda)^2); %��λΪdBm
PI_W = 10^(PI_Calculation_dBm / 10)*10^(-3);


% PI_W = 3.16*10^(-7); %���յ�����Ϊ���Ҹ����źŹ���
A_I = sqrt(PI_W * R_r * 2);
A_I_Amplified = 31.62*A_I;
ADCCode_I = A_I_Amplified*(2^15)/VRef; %��Ϊ���Ҹ��Ŷ�Ӧ����Ƶ�źŲ���������
fi = 78*10^9; %��Ϊ���Ҹ����źŵ�Ƶ��
% fi = 78*10^9; %��Ϊ���Ҹ����źŵ�Ƶ��
% rand����һ�������� (0,1) �ھ��ȷֲ��������
phi = 2*180*rand;


BIF = 15*10^6; %��Ƶ����Ϊ15MHz ��ͨ�˲�������15MHz
tGlitch = BIF/k; %�ܵ����Ҹ���Ӱ��ĳ���ʱ��
t_start_i = (fi - fc)/k; %��ʼ�ܵ����Ҹ���Ӱ���ʱ��

%δ�ܵ�����Ӱ���I��Q��·��Ƶ�źŵ�ʱ����ʽ
s_I = ADCCode_1*cos((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);
s_Q = ADCCode_1*sin((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);

% s_I = ADCCode_1*cos((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);
% s_Q = ADCCode_1*sin((2*pi*k*td_1)*t+2*pi*fc*td_1-pi*k*td_1^2);

figure
plot(t  - t_start,s_I,'Color','blue','LineWidth',4);
xlabel('ʱ��','FontSize',36,'FontWeight','bold');
ylabel('ADCֵ','FontSize',36,'FontWeight','bold');
ylim([-500 500])
hold on;
plot(t  - t_start,s_Q,'Color','red','LineWidth',4);
hold on;
legend('I·��Ƶ�ź�','Q·��Ƶ�ź�','FontSize',28,'FontWeight','bold')



%�ܵ�����Ӱ���I��Q��·��Ƶ�źŵ�ʱ����ʽ
t1 = t_start:1/fs:t_start_i;
t2 = t_start_i:1/fs:(t_start_i+tGlitch);
t3 = (t_start_i+tGlitch):1/fs:(t_start+tADC); %ʱ���Ϊ����ǰ�������źŶΡ����ź�

s_I_1 = ADCCode_1*cos((2*pi*k*td_1)*t1+2*pi*fc*td_1-pi*k*td_1^2);
s_I_2 = ADCCode_1*cos((2*pi*k*td_1)*t2+2*pi*fc*td_1-pi*k*td_1^2) + ADCCode_I*sin(2*pi*(fi-fc)*t2-pi*k*t2.^2 + phi); %���������λ��
s_I_3 = ADCCode_1*cos((2*pi*k*td_1)*t3+2*pi*fc*td_1-pi*k*td_1^2);

s_Q_1 = ADCCode_1*sin((2*pi*k*td_1)*t1+2*pi*fc*td_1-pi*k*td_1^2);
s_Q_2 = ADCCode_1*sin((2*pi*k*td_1)*t2+2*pi*fc*td_1-pi*k*td_1^2) + ADCCode_I*cos(2*pi*(fi-fc)*t2-pi*k*t2.^2 + phi); %�ܵ����ŵ��źŶ�
s_Q_3 = ADCCode_1*sin((2*pi*k*td_1)*t3+2*pi*fc*td_1-pi*k*td_1^2);

t_sum = [t1 t2 t3];
s_I_sum = [s_I_1 s_I_2 s_I_3];
s_Q_sum = [s_Q_1 s_Q_2 s_Q_3];


figure
plot(t_sum - t_start,s_I_sum,'Color','blue','LineWidth',4);
xlabel('ʱ��','FontSize',36,'FontWeight','bold');
ylabel('ADCֵ','FontSize',36,'FontWeight','bold');
hold on
plot(t_sum - t_start,s_Q_sum,'Color','red','LineWidth',4);
hold on
legend('I·��Ƶ�ź�','Q·��Ƶ�ź�','FontSize',28,'FontWeight','bold')

din = s_I + 1j* s_Q; % δ�ܵ����ŵ���ƵIQ�ź�
din_i = s_I_sum + 1j* s_Q_sum; % �ܵ���Ϊ���Ҹ��ŵ���Ƶ�����ź�

range_win = hanning(N); %�ӳ�ΪN=512�ĺ����� 
range_win = range_win';
wSum=sum(range_win(1,:));

temp = din.*range_win;  %�޸����źżӴ�fft
temp_fft =  abs(fftshift(fft(temp, N))); 

temp_i = din_i.*range_win;  %%�����źżӴ�fft
temp_i_fft =  abs(fftshift(fft(temp_i, N)));

index = -N/2 + 0.5: 1: N/2 - 0.5;
freq_bin = (index - 1) * fs / N;
range_bin = freq_bin * c / 2 / k;

Pout_fft = 20*log10(temp_fft);
Pout_fft_i = 20*log10(temp_i_fft);

b_adc = 16;
cf = 20*log10(2^(b_adc-1))+20*log10(wSum)-20*log10(sqrt(2));

Pout_dBFs = Pout_fft - cf;
Pout_dBFs_i = Pout_fft_i - cf;

% Pout_threshold = 20*log10(XT);
% Pout_threshold_dBFs = Pout_threshold - cf;
% 
% 
% i1 = -i_length/2 + 0.5: 1: i_length/2 - 0.5;
% freq_bin_CFAR = (i1 - 1) * fs / N;
% range_bin_CFAR = freq_bin_CFAR * c / 2 / k;
% Pout_CFAR = 20*log10(XT) - cf;

% M=12;
% pro_N=3;
% Q = N-M-pro_N;
% 
% PAD=10^(-4);
% [ind, XT] = cfar_whyu(temp_fft, M, pro_N, PAD);
% Pout_XT = 20*log10(XT);
% Pout_XT_dBFs = Pout_XT-cf;
% Q_ind = -Q/2+0.5:1:Q/2-0.5;
% freq_Q=(Q_ind-1)*fs /Q;
% range_Q = freq_Q*c/2/k;

% %��������ʱ��FFT���
% figure
% plot(range_bin, Pout_dBFs,'b','LineWidth',4);
% % title('fft(after window)')
% xlabel('���루�ף�','FontSize',36,'FontWeight','bold')
% ylabel('FFT�����dBFS��','FontSize',36,'FontWeight','bold')
% axis([-30 30, -250 0])
% hold on

%�ܵ�����ʱ��FFT���
figure
plot(range_bin, Pout_dBFs_i,'b','LineWidth',4);
% title('fft(after window)')
xlabel('���루�ף�','FontSize',36,'FontWeight','bold')
ylabel('FFT�����dBFS��','FontSize',36,'FontWeight','bold')
axis([-30 30, -150 0])
hold on


% plot(range_Q,Pout_XT_dBFs,'r');
% hold on

% result = [range_bin.' Pout_dBFs' Pout_dBFs_i.'];
% result_s = xlswrite('D:\matlabprojects\DifferentTransmittedPower\18dBm-DifferentTransmittedPower-fft.xls', result);    % ��resultд�뵽wind.xls�ļ���
