clc
clear all
close all
%%
N = 200;
lambda = 1064e-6;               %����1064nm
row = linspace(-1,1,N); col = linspace(-1,1,N);
[x,y] = meshgrid(row,col);
[theta,r] = cart2pol(x,y);
w = 3;                          %��˹�����������
k = 2*pi/lambda;                %����
k_r = 20;                       %����ʸ - ����
k_z = sqrt(k^2-k_r^2);          %����ʸ
z = 0;
for n = 0 : 3                    %��������������n = 0,1,2,3�ȵ�
    E = besselj(n,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*n*theta);
    I = E.*conj(E);
    I = I/max(max(I));          %��һ��
    figure;mesh(x,y,I)
    set(gca,'fontname','times new roman','fontsize',16);
    title([num2str(n),'�ױ�����-��˹������ǿ�ֲ�'],'fontname','��������','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
%% ���ɴ��䲿��(�������������)
n = 1;
E1 = besselj(n,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*n*theta);
Z = 100;   %�������
x00 = linspace(-0.5,0.5,N); y00 = linspace(-0.5,0.5,N);
[x2,y2] = meshgrid(x00,y00);
for a=1:N
    for b=1:N
        E2(a,b) = -1i/lambda/Z*exp(1i*k*Z)*sum(sum(E1.*exp(1i*k/2/Z.*((x00(a)-x).^2+(y00(b)-y).^2))));
    end
    a
end
I2 = E2.*conj(E2);      I2 = I2/max(max(I2));
figure;mesh(x2,y2,I2)
set(gca,'fontname','times new roman','fontsize',16);
title([num2str(n),'�ױ�����-��˹�������ɴ�����ǿ�ֲ�'],'fontname','��������','fontsize',16);
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('��һ��ǿ��','fontname','��������','fontsize',16);