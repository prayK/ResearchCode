%% 拉盖尔-高斯光束相位
function main()
%%
clc
clear
close all
%% 参数设置
Nxy = 512;          %x,y采样点
lambda = 632e-9;    %波长为632nm
k = 2*pi/lambda;    %波数
w = 3e-3;           %光斑尺寸
%   /* 坐标设置 */
[x,y] = meshgrid(linspace(-5*w,5*w,Nxy));
[theta,r] = cart2pol(x,y);
%   /* L-G光束参数 */
l = input('请输入拓扑荷数l：');         %拓扑荷数
p = input('请输入p：');                 %p = 0, 1, 2...;
z = input('请输入传输距离l：');         %传输距离-
Z_R = pi*w^2/lambda;        %瑞利长度
w_z = w*sqrt(1+(z/Z_R)^2);  %光束在z位置的半径
E = sqrt(2*factorial(p)/pi/(p+factorial(abs(l))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(l)...
    .*exp(-r.^2/w_z^2).*laguerre(p,abs(l),2*r.^2/w_z^2).*exp(-1i*l*theta).*exp(-1i*k*z)...
    .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(l)+1)*atan(z/Z_R));
I = E.*conj(E);  I = I/max(max(I));
figure; mesh(x*1e3,y*1e3,I); view(2);
set(gca,'fontname','times new roman','fontsize',16);
title(['{\itl}=',num2str(l),',{\itp}=',num2str(p),',{\itz}=',num2str(z),'m处的光强'],'fontname','华文中宋','fontsize',16);
xlabel('\itx/mm','fontname','times new roman','fontsize',16);
ylabel('\ity/mm','fontname','times new roman','fontsize',16);
zlabel('归一化强度','fontname','华文中宋','fontsize',16);
phase = angle(E);   phase = 2*pi*phase/max(max(phase));
figure; mesh(x*1e3,y*1e3,phase); view(2);
colormap gray;
xlabel('\itx/mm','fontname','times new roman','fontsize',16);
ylabel('\ity/mm','fontname','times new roman','fontsize',16);
set(gca,'fontname','times new roman','fontsize',16);
title(['{\itl}=',num2str(l),',{\itp}=',num2str(p),',{\itz}=',num2str(z),'m处的相位'],'fontname','华文中宋','fontsize',16);
zlabel('phase','fontname','华文中宋','fontsize',16);
end
%% 拉盖尔多项式
function result = laguerre(p,l,x)
result = 0;
if p == 0
    result = 1;
elseif p == 1
    result = 1+abs(l)-x;
else
    result = (1/p)*((2*p+l-1-x).*laguerre(p-1,abs(l),x)-(p+l-1)*laguerre(p-2,abs(l),x));
end
end