clc
clear all
close all
%%
N = 200;
lambda = 1064e-6;               %波长1064nm
row = linspace(-1,1,N); col = linspace(-1,1,N);
[x,y] = meshgrid(row,col);
[theta,r] = cart2pol(x,y);
w = 3;                          %高斯光束束腰宽度
k = 2*pi/lambda;                %波数
k_r = 20;                       %径向波矢 - 常量
k_z = sqrt(k^2-k_r^2);          %轴向波矢
z = 0;
for n = 0 : 3                    %贝塞尔函数阶数n = 0,1,2,3等等
    E = besselj(n,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*n*theta);
    I = E.*conj(E);
    I = I/max(max(I));          %归一化
    figure;mesh(x,y,I)
    set(gca,'fontname','times new roman','fontsize',16);
    title([num2str(n),'阶贝塞尔-高斯光束光强分布'],'fontname','华文中宋','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('归一化强度','fontname','华文中宋','fontsize',16);
end
%% 自由传输部分(菲涅尔衍射积分)
n = 1;
E1 = besselj(n,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*n*theta);
Z = 100;   %传输距离
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
title([num2str(n),'阶贝塞尔-高斯光束自由传输后光强分布'],'fontname','华文中宋','fontsize',16);
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);