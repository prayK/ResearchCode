clc;    clear;  close all;
%% 参数设置
Nxy = 200;                      %采样点
lambda = 1064e-9;               %波长1064nm
k = 2*pi/lambda;                %波数
w = 3e-3;                       %高斯光束束腰宽度
[x,y] = meshgrid(linspace(-1.5*w,1.5*w,Nxy));
[theta,r] = cart2pol(x,y);
k_r = 7000;                     %径向波矢 - 常量
k_z = sqrt(k^2-k_r^2);          %轴向波矢
z = 0;
c = 3e8;
T = lambda/c;   %周期
m = 0;          %拓扑荷数
omega = 2*pi*c/lambda;
for jt = 1:9
    theta = imrotate(theta, 90*(jt-1));      %每个周期相位旋转90° 本来是想旋转任意角度的 可惜相位旋转之后矩阵维度不匹配了
    E = besselj(m,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*m*theta).*exp(1i*omega*T*(jt-1));
    I = E.*conj(E);     I = I/max(max(I));
    m = m + 1;          %因为涡旋相位是旋转对称的 相位旋转也看不出来 所以改变拓扑荷数看着更直观一些
    figure;mesh(x*1e3,y*1e3,I);
    view(2);
    set(gca,'fontname','times new roman','fontsize',16);
    title([num2str(jt-1),'T时的光强分布'],'fontname','华文中宋','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('归一化强度','fontname','华文中宋','fontsize',16);
    %set( gca, 'XTick', [], 'YTick', [] );
    F1 = getframe;
    %imwrite(F1.cdata,[ 'C:\Users\Administrator\Desktop\贝塞尔\t = ',num2str(jt-1),'T.png']);  %保存图片
end
