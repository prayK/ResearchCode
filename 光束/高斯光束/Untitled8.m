%光场传输的各种算法
clc
close all
clear all
%% 光源部分
No = 200;           %取样点数
N  = 128;
lambda = 1064e-6;   %波长1064nm
k = 2*pi/lambda;    %波矢
w = 3;              %高斯光束的束宽
d = 1000;           %传输距离1000mm
[x0,y0] = meshgrid(linspace(-1.5*w,1.5*w,No));
E0 = exp(-(x0.^2+y0.^2)/w^2);           %高斯光束
I0 = E0.*conj(E0);   I0 = I0/max(max(I0));
figure;mesh(x0,y0,I0)
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);
%%  光传输矩阵相乘算法
dx0 = x0(1,2)-x0(1,1);  dy0 = dx0;
[x1,y1] = meshgrid(linspace(-1.5*w,1.5*w,N));
Mx = exp(-1i*k/d*x1(1,:)'*x0(1,:));     
My = exp(-1i*k/d*y0(:,1)*y1(:,1)');          
M  = E0.*exp(1i*k/2/d*(x0.^2+y0.^2));     
E1 = -1i/lambda/d*exp(1i*k*d)*exp(1i*k/2/d*(x1.^2+y1.^2)).*(Mx*M*My);
I1 = E1.*conj(E1);   I1 = I1/max(max(I1));
figure;mesh(x1,y1,I1);
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);
%% 积分算法
for a = 1:N
    for b = 1:N
        E2(a,b) = -1i/lambda/d*exp(1i*k*d)*sum(sum(E0.*exp(1i*k/2/d*((x1(a,b)-x0).^2+(y1(a,b)-y0).^2))));
    end
    a
end
I2 = E2.*conj(E2);      I2 = I2/max(max(I2));
figure;mesh(x1,y1,I2)
set(gca,'fontname','times new roman','fontsize',16);
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);