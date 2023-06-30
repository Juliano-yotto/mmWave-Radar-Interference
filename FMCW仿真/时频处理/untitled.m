clc;
clear;
close all;
 
%% 仿真信号
figure;
z=amgauss(160,90,40); 
subplot(3,1,1);plot(z);title('常用产生信号的函数--高斯幅值调制信号');
 
z=fmconst(128,0.05,50);
subplot(3,1,2);plot(real(z));title('常用产生信号的函数--固定频率的频率调制信号');
 
[z, f]=fmlin(128,0.05,0.3,50); 
subplot(3,1,3);plot(real(z));title('常用产生信号的函数--线性调频信号');
 
%% Wigner-Ville时频分布图
sig=amgauss(160,90,40);
figure;
tfrwv(sig);
 
sig=fmconst(128,0.05,50);
figure;
tfrwv(sig);
 
sig = fmlin(128,0.1,0.4);
figure;
tfrwv(sig);
 
sig1 = fmlin(128,0.1,0.4);
sig2 = fmlin(128,0.2,0.5);
sig = sig1+sig2;
figure;
tfrwv(sig);
 
%% 伪Wigner-Ville时频分布图
sig = fmlin(128,0.1,0.4);
figure;
tfrpwv(sig);
 
sig1 = fmlin(128,0.1,0.4);
sig2 = fmlin(128,0.2,0.5);
sig = sig1+sig2;
figure;
tfrpwv(sig);