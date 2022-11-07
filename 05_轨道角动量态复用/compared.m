clc;
clear;
close all;
warning off;

load R.mat

figure;
loglog(Cn2,P1,'b-o');
hold on
loglog(Cn2(3:end),P2,'r-s');
hold on
grid on
legend('�޲�ǰУ��ģ��','����ǰУ��ģ��');

xlabel('Cn2');
ylabel('Ber');
grid on

axis([5e-16,1e-13,1e-6,1]);