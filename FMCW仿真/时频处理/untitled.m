clc;
clear;
close all;
 
%% �����ź�
figure;
z=amgauss(160,90,40); 
subplot(3,1,1);plot(z);title('���ò����źŵĺ���--��˹��ֵ�����ź�');
 
z=fmconst(128,0.05,50);
subplot(3,1,2);plot(real(z));title('���ò����źŵĺ���--�̶�Ƶ�ʵ�Ƶ�ʵ����ź�');
 
[z, f]=fmlin(128,0.05,0.3,50); 
subplot(3,1,3);plot(real(z));title('���ò����źŵĺ���--���Ե�Ƶ�ź�');
 
%% Wigner-VilleʱƵ�ֲ�ͼ
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
 
%% αWigner-VilleʱƵ�ֲ�ͼ
sig = fmlin(128,0.1,0.4);
figure;
tfrpwv(sig);
 
sig1 = fmlin(128,0.1,0.4);
sig2 = fmlin(128,0.2,0.5);
sig = sig1+sig2;
figure;
tfrpwv(sig);