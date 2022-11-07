clc
clear all
close all
lambda=1.064e-6;
k=2*pi/lambda;
x=(-1:0.01:1)*1e-4;
y=x;            
[X,Y]=meshgrid(x,y);
phi=atan2(Y,X);

for l=-4:4
    E1=exp(-i*k*X);
    E2=exp(i*l*phi);
    c=E1+E2;
    E=c.*conj(c);
    E=(E-min(min(E)))/(max(max(E))-min(min(E)));
    subplot(3,3,l+5)
    pcolor(x,y,E)
    
    shading interp%通过在每个线条或面中对颜色图索引或真彩色值进行插值来改变该线条或面中的颜色。
    colormap hot%颜色
    axis off%去掉坐标轴
    title(num2str(l))%将数字转换为字符数组。
    %waitforbuttonpress
end
