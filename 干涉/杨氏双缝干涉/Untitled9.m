%杨氏双缝干涉
clc
close all
clear all
%%
lambda = 500e-9;    %%波长500nm
d = 2e-3;           %%双缝间距2mm
D=1;                %%双缝距离观察屏之间的距离1m
ym=5*lambda*D/d;
xs=ym;
n=101;
ys=linspace(-ym,ym,n);              %%观察屏面
for i=1:n
   r1=sqrt((ys(i)-d/2).^2+D^2);     %%光程1
   r2=sqrt((ys(i)+d/2).^2+D^2);     %%光程2
   phi=2*pi*(r2-r1)./lambda;        %%相位差
   B(i,:)=sum(4*cos(phi/2).^2);     %%计算光强
end
N=255;
Br=(B/4.0)*N;
subplot(1,2,1)      %%一行两列中的第一幅图
image(xs,ys,Br);
colormap(gray(N));  %%使图像以灰度图像显示
subplot(1,2,2)      %%一行两列中的第二幅图
plot(B,ys)