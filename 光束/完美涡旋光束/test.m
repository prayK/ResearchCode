clc
clear all
close all
%% 振幅光学相位元件产生完美涡旋光场
N = 300;
row = linspace(-1.5e-6,1.5e-6,N);   
col = linspace(-1.5e-6,1.5e-6,N);
[x,y] = meshgrid(row,col);
[phi,rho] = cart2pol(x,y);
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
f = 0.5;            %透镜焦距0.5m
R = 0.5;            %圆孔光阑半径0.3m
X = k*rho/f;
alpha = 20;          %径向波矢 - 常数

figure;
for m = 0 : 8
    subplot(3,3,m+1)
    E1 = (-1i)^(m+1)*(k*R/f)*exp(1i*m*phi).*...
            ((alpha*besselj(m+1,alpha*R)*besselj(m,X*R)-X.*besselj(m,alpha*R).*besselj(m+1,X*R))./(alpha^2-X.^2));  
    I1 = E1.*conj(E1);  I1 = I1/max(max(I1));
    %二维
    h1 = pcolor(x,y,I1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    %colormap(gray);        %输出灰度图像
    axis square;
    
    %三维
%     mesh(x,y,I1)           %三维
%     set(gca,'fontname','times new roman','fontsize',16);
%     title(['m = ',num2str(m)],'fontname','华文中宋','fontsize',16);
%     xlabel('x/m','fontname','times new roman','fontsize',16);
%     ylabel('y/m','fontname','times new roman','fontsize',16);
%     zlabel('归一化强度','fontname','华文中宋','fontsize',16);
end
suptitle('振幅光学相位元件产生完美涡旋光场')   %为图一添加总标题
%% 利用锥透镜产生完美涡旋光场
row = linspace(-0.8,0.8,N);   
col = linspace(-0.8,0.8,N);
[x,y] = meshgrid(row,col);
[phi,rho] = cart2pol(x,y);
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
f = 0.5;            %透镜焦距0.5m
R = 0.5;            %圆孔光阑半径0.3m
alpha = 20;          %径向波矢 - 常数
w0 = 0.3;           %透镜焦平面上的束腰
w_g = 0.5;          %高斯项的束腰
figure;
for m = 0 : 8
    subplot(3,3,m+1)
    E2 =  (1i)^(m+1)*(w_g/w0)*exp(1i*m*phi).*exp(-(rho.^2+R^2)/w0^2).*besseli(m,2*R*rho/w0^2);
    I2 = E2.*conj(E2);  I2 = I2/max(max(I2));
    %二维
    h2 = pcolor(x,y,I2);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    %colormap(gray);        %输出灰度图像
    axis square;
    
    %三维
%     mesh(x,y,I2)           %三维
%     set(gca,'fontname','times new roman','fontsize',16);
%     title(['m = ',num2str(m)],'fontname','华文中宋','fontsize',16);
%     xlabel('x/m','fontname','times new roman','fontsize',16);
%     ylabel('y/m','fontname','times new roman','fontsize',16);
%     zlabel('归一化强度','fontname','华文中宋','fontsize',16);
end
suptitle('利用锥透镜产生完美涡旋光场')   %为图二添加总标题