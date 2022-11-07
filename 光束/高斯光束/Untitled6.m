%高斯光束
clc
clear all
close all
%%  
N = 100;
lambda = 1064e-6;   %波长为1064nm
k = 2*pi/lambda;    %波矢
w = 3;              %高斯光束的束宽
[x1,y1] = meshgrid(linspace(-10,10,N));
E1 = exp(-(x1.^2+y1.^2)/w^2);   
I1 = E1.*conj(E1);
I1 = I1/max(max(I1));
figure;mesh(x1,y1,I1)
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);
%% 传输程序
f = 500;
A = 0;  B = f;  D = 0;
x = linspace(-0.1,0.1,N); y = linspace(-0.1,0.1,N);
[x2,y2] = meshgrid(x,y);
for a = 1:N
    for b = 1:N
        E2(a,b) = sum(sum(E1.*exp(i*k/2/B*(A*(x1.^2+y1.^2)+D*(x(a).^2+y(b).^2)-2*(x1.*x(a)+y1.*y(b))))));
    end
    a
end
I2 = E2.*conj(E2);
I2 = I2/max(max(I2));
figure;mesh(x2,y2,I2)
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
xlabel('x/mm','fontname','times new roman','fontsize',16);
ylabel('y/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);