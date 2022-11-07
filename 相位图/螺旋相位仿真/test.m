clc;clear all;close all;

R1=20;
l=3;
r=linspace(0,R1,24);
for j=1:l
    fai=linspace(2*(j-1)*pi/l,2*j*pi/l,100);
    [R,T]=meshgrid(r,fai);
    x=cos(fai')*r;
    y=sin(fai')*r;
    z=l*T/(2*pi)-j+1;
    m=ones(size(z))*100;
    figure(2),surf(x,y,z,m);
    hold on;
end
hold off;
colormap([1 0 0;1 0 0;1 0 0]);
% axis off;

%^ÂÝÐýÏàÎ»Îª0
% R1=20;
% h=0;
% r=linspace(0,R1,24);
% fai=linspace(0,2*pi,100);
% [R,T]=meshgrid(r,fai);
% x=cos(fai')*r;
% y=sin(fai')*r;
% z=h*T/(2*pi);
% m=ones(size(z))*100;
% surf(x,y,z,m);
% colormap([1 0 0;1 0 0;1 0 0]);
% axis off