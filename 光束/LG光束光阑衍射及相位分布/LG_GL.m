%拉盖尔-高斯光束光阑衍射
function main()
clc;
clear;
close all;
%% /* 参数设置 */
Nxy = 200;          %x,y坐标采样点
lambda = 632e-9;    %波长为632nm
k = 2*pi/lambda;    %波数
w = 0.5e-3;         %光斑尺寸 0.5mm
p = 1;              %径向标量 p = 0, 1, 2......
%   /* 坐标设置 */
[x,y] = meshgrid(linspace(-3*w,3*w,Nxy));       %近场直角坐标
[theta,r] = cart2pol(x,y);                      %近场极坐标
[xf,yf] = meshgrid(linspace(-3*w,3*w,Nxy));     %远场坐标
z = 0;
Z_R = pi*w^2/lambda;           %瑞利长度
w_z = w*sqrt(1+(z/Z_R)^2);     %光束在z位置的半径
%% 衍射部分
Z = 0.1;            %传输距离0.5m
Length = 0.2e-3;    %矩形光阑0.2mm*0.2mm
Width = Length;
for m = 1:6
    tic
    %   /* 不同拓扑荷数的L-G光束 */
    E0 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*Laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    I0 = E0.*conj(E0);  I0 = I0/max(max(I0));
    figure(1);subplot(2,3,m);
    mesh(x*1e3,y*1e3,I0);
    view(2);
    set(gca,'fontname','times new roman','fontsize',16);
    title(['l = ',num2str(m)],'fontname','华文中宋','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('归一化强度','fontname','华文中宋','fontsize',16);
    if m == 6
        suptitle(['p = ',num2str(p),'时，不同拓扑荷数对应的L-G光束']);
    end
    %   /* 矩形光阑衍射 */
    %     E = E0.*RectGrating(E0,Length,Width,6*w,6*w);
    %     for a = 1:Nxy
    %         for b = 1:Nxy
    %             E1(a,b) = (-1i/lambda/Z)*exp(1i*k*Z)*...
    %                 sum(sum(E.*exp(1i*k/2/Z*((xf(a,b)-x).^2+(yf(a,b)-y).^2))));
    %         end
    %     end
    %     I1 = E1.*conj(E1);  I1 = I1/max(max(I1));
    %     figure(2);subplot(2,3,m)
    %     mesh(xf*1e3,yf*1e3,I1);
    %     colormap gray;
    %     view(2);
    %     set(gca,'fontname','times new roman','fontsize',16);
    %     title(['m = ',num2str(m)],'fontname','华文中宋','fontsize',16);
    %     xlabel('x/mm','fontname','times new roman','fontsize',16);
    %     ylabel('y/mm','fontname','times new roman','fontsize',16);
    %     zlabel('归一化强度','fontname','华文中宋','fontsize',16);
    %     if m == 6
    %         suptitle('不同拓扑荷数的L-G光束的矩形光阑衍射图样');
    %     end
    %   /* 圆形光阑衍射 */
    E = E0.*CircularGrating(E0,6*w,6*w,w);
    for a = 1:Nxy
        for b = 1:Nxy
            E2(a,b) = (-1i/lambda/Z)*exp(1i*k*Z)*...
                sum(sum(E.*exp(1i*k/2/Z*((xf(a,b)-x).^2+(yf(a,b)-y).^2))));
        end
    end
    I2 = E2.*conj(E2);  I2 = I2/max(max(I2));
    figure(3);subplot(2,3,m)
    mesh(xf*1e3,yf*1e3,I2);
    colormap gray;
    view(2);
    set(gca,'fontname','times new roman','fontsize',16);
    title(['m = ',num2str(m)],'fontname','华文中宋','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('归一化强度','fontname','华文中宋','fontsize',16);
    if m == 6
        suptitle('不同拓扑荷数的L-G光束的圆形光阑衍射图样');
    end
    toc
end
end
%% 拉盖尔多项式
function [result] = Laguerre(p,l,x)
if p == 0
    result = 1;
elseif p == 1
    result = 1+abs(l)-x;
else
    result = (1/p)*((2*p+l-1-x).*Laguerre(p-1,abs(l),x)-(p+l-1)*Laguerre(p-2,abs(l),x));
end
end
%% 光阑透过率函数
%   /* 矩形光阑 */
function [ result ] = RectGrating(E0,Length,Width,Lx,Ly)
%Length Width 为光阑长宽 Lx Ly为采样点x,y方向长度
[Nx,Ny] = size(E0);
result = zeros(Nx,Ny);
X_begin = ceil((1-Length/Lx)*Nx/2);
X_end = ceil((1+Length/Lx)*Nx/2);
Y_begin = ceil((1-Width/Ly)*Ny/2);
Y_end = ceil((1+Width/Ly)*Ny/2);
for a = X_begin:X_end
    for b = Y_begin:Y_end
        result(a,b) = 1;
    end
end
end
%   /* 圆形光阑 */
%Lx Ly为采样点x,y方向长度 R为光阑半径
function [ result ] = CircularGrating(E0,Lx,Ly,R)
[Nx,Ny] = size(E0);
[x,y] = meshgrid(linspace(-Lx/2,Lx/2,Nx),linspace(-Ly/2,Ly/2,Ny));
[theta,r] = cart2pol(x,y);
result = zeros(Nx,Ny);
result(r<R) = 1;
end