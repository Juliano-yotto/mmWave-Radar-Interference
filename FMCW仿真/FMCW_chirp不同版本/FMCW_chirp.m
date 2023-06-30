%--2022/10/2--%
%--�--%
%�������FMCW���״���١������з��档Ϊ�������п����ŷ����̵�

%% ���⣺
%% 1.ʵ�����ʱ�����ջ��յ��ز��źŲſ�ʼ������������ʼʱ�����ȷ����ֻ���м�һ�Σ�
%% 2.matlabд����ʱ��Ҫ��������﷢���ź�������źţ���ôʱ����Ӧ�ü�������      2022/10/6 �ѽ����
%% 3.ʱ���ӳ�tao�Ƿ�Ҫ��ɱ�����������������Ϣ���˳���д�ɱ���td(i)������ʦ��ĳ��������Ƕ�ֵtd.TI�ĵ���Ҳ���ƺ���IF�ź�Ƶ���������ٶ����������ɶ�ֵ
%% 4.��θ�ÿ��chirp�źŵ��������ұ��������ź������䣬������ɱ䣿
%% 5.�ÿ���ʦ������� �������
%% 2022/10/23 ���˽����źŴ�����wr_first;wr
%% sr���˸����Ҫȥ����

clear all;
close all;
clc;

%% Ŀ������
R_target_ini=100;
V_target=10;

%% �״��������   
fc=77e9; %�ز�Ƶ��
c=3e8; %����

%%-----------�����������----------------
% % slope = 25.49*10^12;
% % T_chirp=125.07e-6; %chirp����
% % B_chirp=slope*T_chirp; %�����źŵ�Ƶ����

Endle_time=6.3e-6; %chirp����ʱ��
T_chirp=7.33e-6; %chirp����
B_chirp=150e6; %�����źŵ�Ƶ����
slope=B_chirp/T_chirp; %���Ե�Ƶб��68MHz/us

N_chirp=128; %chirp����
N_ADC=1024; %chirpһ��������ADC��������
fs=N_ADC/T_chirp; %��Ƶģ���źŲ����� 

dres=c/(2*B_chirp); %����ֱ���
vres=(c/fc)/(2*N_chirp*(T_chirp+Endle_time)); %�ٶȷֱ���

num_zeros=round(Endle_time*fs); %�������
N_total=N_ADC+num_zeros; %һ�������� chirp�����ź�+����ʱ���ܵ���

t_0=linspace(0,Endle_time,num_zeros); %�����ź�ʱ����
t_use=linspace(0,T_chirp,N_ADC); %�����źŲ���ʱ����
t=linspace(0,T_chirp+Endle_time,N_total); %һ���������ź����� ��ʱ����
t_target=linspace(0,(T_chirp+Endle_time)*N_chirp,N_chirp*N_total); %Ŀ���ƶ�����ʱ����
%t=0:1/fs:N_chirp*(T_chirp+Endle_time); %����ʱ����

st=zeros(1,length(t_target));
rt=zeros(1,length(t_target));
st_zeros=0*t_0; %���貹��㻯Ϊ0�����飩�����㲹��ȥ
mix_IF=zeros(1,length(t_target));

%% ���䡢�����źŽ�ģ
R_target = R_target_ini+V_target*t_target; %Ŀ��λ�ø�������
td = 2*R_target/c; %�ز��ź��ӳ�����
fd = 2 * (fc - B_chirp/2) * V_target / c; %������Ƶ��
for i=1:N_chirp %��i��chirp�ź�
%     for j=1:N_ADC
%     R_target(j)=R_target_ini+V_target*t(j);
%     td(j)=2*R_target(j)/c;
%     sig_tx_use(i) = cos(2*pi*(fc*t(j)+slope*t(j)^2/2));
%     sig_rx_use(i) = cos(2*pi*(fc*(t(j)-td(j))+slope*(t(j)-td(j))^2/2));
%     sig_tx{1,i} = [Sig_tx_use{1,i},st_zeros];
%     end
if i==1 %��һ�����䡢����chirp�ź�
    Wt_first = rectpuls(t_target-T_chirp/2,T_chirp); %��һ�����䴰���� 
    st_use_first = cos(2*pi*(fc*t_target+slope*t_target.^2/2)).*Wt_first; %�����źŵ�һ������
%    st_temp_first = [st_use_first(1:N_ADC),st_zeros]; %��һ������chirp�źŲ���
    st = st + st_use_first; %�źŵ���
%    st=st_temp_first;
    
    Wr_first = rectpuls(t_target-td-T_chirp/2,T_chirp); %��һ�������źŵĴ�����
    %Wr_first = rectpuls(t_target-T_chirp/2,T_chirp);
    rt_use_first = cos(2*pi*(fc*(t_target-td)+slope*(t_target-td).^2/2)).*Wr_first; %�ز��ź� ����������targetʱ��
%    rt_temp_first = [rt_use_first(1:N_ADC),st_zeros]; %��һ���ز��źŲ���
    rt = rt + rt_use_first;
%    rt=rt_temp_first;
    %IF_first = st_use_first.*rt_use_first;
    IF_first = st_use_first.*conj(rt_use_first);
    IF_first_use = IF_first(1:1024);
else 
    Wt = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2),T_chirp); %����chirp�Ĵ������ӳ�
    %st_temp = cos(2*pi*(fc*(t_target-(i-1)*(T_chirp+Endle_time))+slope*(t_target-(i-1)*(T_chirp+Endle_time)).^2/2)).*Wt; %�����ź��ӳ� ����ʱ��������target
    st_temp = cos(2*pi*(fc*(t_target-(i-1)*(T_chirp+Endle_time))+slope*t_target.^2/2)).*Wt; %�����ź��ӳ� ����ʱ��������target
    
    st = st+ st_temp; %�����ź�����ֱ�ӵ���
%    st_temp = [st_temp_use((i-1)*(N_ADC+num_zeros)+1:i*(N_ADC+num_zeros)),st_zeros]; %���ò�������������ٲ��� ���һ������1024+52
%    st = [st,st_temp]; %�����ź�stʱ�������ӳ�����
    
    Wr = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2+td),T_chirp);
%    Wr = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2),T_chirp);
    %rt_temp = cos(2*pi*(fc*(t_target-((i-1)*(T_chirp+Endle_time)+td))+slope*(t_target-((i-1)*(T_chirp+Endle_time)+td)).^2/2)).*Wr; %�ز��ź��ӳ�
    rt_temp = cos(2*pi*(fc*(t_target-((i-1)*(T_chirp+Endle_time)+td))+slope*(t_target-td).^2/2)).*Wr;
    rt = rt + rt_temp; %�ز��ź����ڵ���

%     sig_tx{1,i} = [sig_tx_use,st_zeros];
%     sig_rx{1,i} = [sig_rx_use,st_zeros];
end
end

% ���䡢�����ź�ʱ��
figure (1)
plot(t_target,st);
title('�����ź�ʱ��');
figure (2)
plot(t_target,rt);
title('�ز��ź�ʱ��');

%% �����źŻ�Ƶ�õ���Ƶ�ź�
mix_IF=st.*rt;
%mix_IF=st.*conj(rt);
%mix_IF = rt.*(st').'; 
signal_IF_use = reshape(mix_IF,N_ADC+num_zeros,N_chirp); %�ź��ع� 137728=(1024+52)*128 ��Ƶ�źž��� Ϊ[��������(�������)��chirp����] (1024+52)*128
signal_IF = signal_IF_use(1:1024,:) ;

figure (3)
mesh(signal_IF);
xlabel('������')
ylabel('��������');
title('��Ƶ�ź�ʱ��');

%% ����άFFT
sig_fft = fft(signal_IF,N_ADC)./N_ADC;
sig_fft = abs(sig_fft);
sig_fft = sig_fft(1:(N_ADC/2),:);

figure (4)
plot(sig_fft(:,1));
xlabel('���루Ƶ�ʣ�');
ylabel('����')
title('��һ��chirp��FTF���')

figure (5)
mesh(sig_fft);
xlabel('���루Ƶ�ʣ�');
ylabel('chirp������')
zlabel('����')
title('����άFTF���')

%% �ٶ�άFFT
sig_fft2 = fft2(mix_IF,N_ADC,N_chirp); %�Ծ��롢�ٶ������������FFT
sig_fft2 = sig_fft2(1:N_ADC/2,1:N_chirp);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM); %ת��Ϊ�ֱ�
doppler_axis = linspace(-100,100,N_chirp);
range_axis = linspace(-200,200,N_ADC/2)*((N_ADC/2)/400);

figure (6)
mesh(doppler_axis,range_axis,RDM);
xlabel('������ͨ��'); ylabel('����ͨ��'); zlabel('���ȣ�dB��');
title('�ٶ�άFFT �����������');

%figure (7)
%plot(db(abs(fft(IF_first_use(1:1024)))));
