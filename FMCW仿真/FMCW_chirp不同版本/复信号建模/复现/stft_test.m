%% 用这个语句出图
% 时频分析工具箱测试使用
clc
clear
close all
% 构造一个简单信号
T=1;%1s
fs =1024;
N=fs*T;
ts = T/N;
t = 0:ts:N*ts;
% yt1=exp(-j*70*sin(3*pi*t));
yt2=exp(j*2*pi*(-70*t+150*t.^2/2));
%yt = exp(-j*70*sin(3*pi*t))+exp(j*90*pi*t)+exp(-j*90*pi*t);
%yt = yt - mean(yt);
yh=yt2;
[tfr2,t2,f2]=tfrstft(yh(1:N).',1:N,N,hamming(127));
df=fs*(0:N-1)/N - fs/2;
tfr=fftshift(tfr2,1);
figure
imagesc(t(1:end-1), df(2:end), abs(tfr)); axis xy
ylim([-150 150]);
axis([0,1024*ts,-200,200]);

%% 
% [tfr,t,f] = tfrstft(yh.', 1:length(yh),500); % 1024为频率分辨率，就是想把频率轴分多细
% figure,surf(t*ts,f*fs,abs(tfr),'edgecolor','none');
% axis tight, view(0,90);
% ylim([-150 150]);%确定频率范围