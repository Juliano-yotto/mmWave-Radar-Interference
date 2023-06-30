%--2022/10/2--%
%--杨帆--%
%本程序对FMCW波雷达测速、测距进行仿真。为后续进行抗干扰仿真铺垫

%% 问题：
%% 1.实际情况时，接收机收到回波信号才开始采样。采样开始时间如何确定？只采中间一段？
%% 2.matlab写程序时，要想完整表达发射信号与接收信号，那么时间轴应该加以区分      2022/10/6 已解决√
%% 3.时间延迟tao是否要算成变量？包含多普勒信息？此程序写成变量td(i)，但是师姐的程序里面是定值td.TI文档中也近似忽略IF信号频率与物体速度依赖，看成定值
%% 4.如何给每个chirp信号点添加零点且保持有用信号数不变，零点数可变？
%% 5.用可意师姐的数据 结果错误
%% 2022/10/23 改了接收信号窗函数wr_first;wr
%% sr加了复共轭（要去掉）

clear all;
close all;
clc;

%% 目标设置
R_target_ini=100;
V_target=10;

%% 雷达参数设置   
fc=77e9; %载波频率
c=3e8; %光速

%%-----------有问题的数据----------------
% % slope = 25.49*10^12;
% % T_chirp=125.07e-6; %chirp周期
% % B_chirp=slope*T_chirp; %发射信号调频带宽

Endle_time=6.3e-6; %chirp空闲时间
T_chirp=7.33e-6; %chirp周期
B_chirp=150e6; %发射信号调频带宽
slope=B_chirp/T_chirp; %线性调频斜率68MHz/us

N_chirp=128; %chirp个数
N_ADC=1024; %chirp一个周期内ADC采样点数
fs=N_ADC/T_chirp; %中频模拟信号采样率 

dres=c/(2*B_chirp); %距离分辨率
vres=(c/fc)/(2*N_chirp*(T_chirp+Endle_time)); %速度分辨率

num_zeros=round(Endle_time*fs); %补零点数
N_total=N_ADC+num_zeros; %一个周期中 chirp有用信号+空闲时间总点数

t_0=linspace(0,Endle_time,num_zeros); %补零信号时间轴
t_use=linspace(0,T_chirp,N_ADC); %有用信号采样时间轴
t=linspace(0,T_chirp+Endle_time,N_total); %一整个发射信号周期 的时间轴
t_target=linspace(0,(T_chirp+Endle_time)*N_chirp,N_chirp*N_total); %目标移动连续时间轴
%t=0:1/fs:N_chirp*(T_chirp+Endle_time); %连续时间轴

st=zeros(1,length(t_target));
rt=zeros(1,length(t_target));
st_zeros=0*t_0; %将需补零点化为0（数组），方便补进去
mix_IF=zeros(1,length(t_target));

%% 发射、接收信号建模
R_target = R_target_ini+V_target*t_target; %目标位置更新数组
td = 2*R_target/c; %回波信号延迟数组
fd = 2 * (fc - B_chirp/2) * V_target / c; %多普勒频移
for i=1:N_chirp %第i个chirp信号
%     for j=1:N_ADC
%     R_target(j)=R_target_ini+V_target*t(j);
%     td(j)=2*R_target(j)/c;
%     sig_tx_use(i) = cos(2*pi*(fc*t(j)+slope*t(j)^2/2));
%     sig_rx_use(i) = cos(2*pi*(fc*(t(j)-td(j))+slope*(t(j)-td(j))^2/2));
%     sig_tx{1,i} = [Sig_tx_use{1,i},st_zeros];
%     end
if i==1 %第一个发射、接收chirp信号
    Wt_first = rectpuls(t_target-T_chirp/2,T_chirp); %第一个发射窗函数 
    st_use_first = cos(2*pi*(fc*t_target+slope*t_target.^2/2)).*Wt_first; %发射信号第一个周期
%    st_temp_first = [st_use_first(1:N_ADC),st_zeros]; %第一个周期chirp信号补零
    st = st + st_use_first; %信号叠加
%    st=st_temp_first;
    
    Wr_first = rectpuls(t_target-td-T_chirp/2,T_chirp); %第一个接收信号的窗函数
    %Wr_first = rectpuls(t_target-T_chirp/2,T_chirp);
    rt_use_first = cos(2*pi*(fc*(t_target-td)+slope*(t_target-td).^2/2)).*Wr_first; %回波信号 长度是整个target时长
%    rt_temp_first = [rt_use_first(1:N_ADC),st_zeros]; %第一个回波信号补零
    rt = rt + rt_use_first;
%    rt=rt_temp_first;
    %IF_first = st_use_first.*rt_use_first;
    IF_first = st_use_first.*conj(rt_use_first);
    IF_first_use = IF_first(1:1024);
else 
    Wt = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2),T_chirp); %发射chirp的窗函数延迟
    %st_temp = cos(2*pi*(fc*(t_target-(i-1)*(T_chirp+Endle_time))+slope*(t_target-(i-1)*(T_chirp+Endle_time)).^2/2)).*Wt; %发射信号延迟 点数时长是整个target
    st_temp = cos(2*pi*(fc*(t_target-(i-1)*(T_chirp+Endle_time))+slope*t_target.^2/2)).*Wt; %发射信号延迟 点数时长是整个target
    
    st = st+ st_temp; %发射信号周期直接叠加
%    st_temp = [st_temp_use((i-1)*(N_ADC+num_zeros)+1:i*(N_ADC+num_zeros)),st_zeros]; %有用采样点数拆出来再补零 变成一个周期1024+52
%    st = [st,st_temp]; %发射信号st时域序列延长叠加
    
    Wr = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2+td),T_chirp);
%    Wr = rectpuls(t_target-((i-1)*(T_chirp+Endle_time)+T_chirp/2),T_chirp);
    %rt_temp = cos(2*pi*(fc*(t_target-((i-1)*(T_chirp+Endle_time)+td))+slope*(t_target-((i-1)*(T_chirp+Endle_time)+td)).^2/2)).*Wr; %回波信号延迟
    rt_temp = cos(2*pi*(fc*(t_target-((i-1)*(T_chirp+Endle_time)+td))+slope*(t_target-td).^2/2)).*Wr;
    rt = rt + rt_temp; %回波信号周期叠加

%     sig_tx{1,i} = [sig_tx_use,st_zeros];
%     sig_rx{1,i} = [sig_rx_use,st_zeros];
end
end

% 发射、接收信号时域
figure (1)
plot(t_target,st);
title('发射信号时域');
figure (2)
plot(t_target,rt);
title('回波信号时域');

%% 接收信号混频得到拍频信号
mix_IF=st.*rt;
%mix_IF=st.*conj(rt);
%mix_IF = rt.*(st').'; 
signal_IF_use = reshape(mix_IF,N_ADC+num_zeros,N_chirp); %信号重构 137728=(1024+52)*128 中频信号矩阵 为[采样点数(含补零点)×chirp个数] (1024+52)*128
signal_IF = signal_IF_use(1:1024,:) ;

figure (3)
mesh(signal_IF);
xlabel('脉冲数')
ylabel('距离门数');
title('中频信号时域');

%% 距离维FFT
sig_fft = fft(signal_IF,N_ADC)./N_ADC;
sig_fft = abs(sig_fft);
sig_fft = sig_fft(1:(N_ADC/2),:);

figure (4)
plot(sig_fft(:,1));
xlabel('距离（频率）');
ylabel('幅度')
title('第一个chirp的FTF结果')

figure (5)
mesh(sig_fft);
xlabel('距离（频率）');
ylabel('chirp脉冲数')
zlabel('幅度')
title('距离维FTF结果')

%% 速度维FFT
sig_fft2 = fft2(mix_IF,N_ADC,N_chirp); %对距离、速度两个方向进行FFT
sig_fft2 = sig_fft2(1:N_ADC/2,1:N_chirp);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM); %转化为分贝
doppler_axis = linspace(-100,100,N_chirp);
range_axis = linspace(-200,200,N_ADC/2)*((N_ADC/2)/400);

figure (6)
mesh(doppler_axis,range_axis,RDM);
xlabel('多普勒通道'); ylabel('距离通道'); zlabel('幅度（dB）');
title('速度维FFT 距离多普勒谱');

%figure (7)
%plot(db(abs(fft(IF_first_use(1:1024)))));
