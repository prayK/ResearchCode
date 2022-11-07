%���źͷ�����
clc
clear all
close all
%% ��Դ����
N = 200;
Diameter = 10;
x = linspace(-5,5,N); y = linspace(-5,5,N);
[x1,y1] = meshgrid(x,y);
E1 = ones(N);
E1(sqrt(x1.^2+y1.^2)>Diameter/2) = 0;
I1 = E1.*conj(E1);
I1 = I1/(max(max(I1)));
figure;mesh(x1,y1,I1)
set(gca,'fontname','times new roman','fontsize',16);
title('��Դ��ǿ�ֲ�','fontname','��������','fontsize',16);
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
%% ���ź̷�����(��λ��mm)
f = 5000;            %%͸������
lambda = 1064e-6;   %%����
k = 2*pi/lambda;    
z = 1000;
A = 0;  B = f;  C = -1/f;   D = 1-1/f;  %%ABCD����
xoo = linspace(-1,1,N);  yoo = linspace(-1,1,N);
[x2,y2] = meshgrid(xoo,yoo);
[THETA,R] = cart2pol(x2,y2);
E2 = -1i*k*Diameter^2/8/f*exp(1i*k*z).*exp(1i*k/2/f.*R.^2).*(2*(besselj(1,k*Diameter*R/2*f))./(k*Diameter*R/2*f));
I2 = E2.*conj(E2);  
I2 = I2/max(max(I2));
figure;mesh(x2,y2,I2)
set(gca,'fontname','times new roman','fontsize',16);
title('���ź̷������Ĺ�ǿ�ֲ�','fontname','��������','fontsize',16);
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('��һ��ǿ��','fontname','��������','fontsize',16);