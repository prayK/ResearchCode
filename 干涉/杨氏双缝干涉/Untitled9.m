%����˫�����
clc
close all
clear all
%%
lambda = 500e-9;    %%����500nm
d = 2e-3;           %%˫����2mm
D=1;                %%˫�����۲���֮��ľ���1m
ym=5*lambda*D/d;
xs=ym;
n=101;
ys=linspace(-ym,ym,n);              %%�۲�����
for i=1:n
   r1=sqrt((ys(i)-d/2).^2+D^2);     %%���1
   r2=sqrt((ys(i)+d/2).^2+D^2);     %%���2
   phi=2*pi*(r2-r1)./lambda;        %%��λ��
   B(i,:)=sum(4*cos(phi/2).^2);     %%�����ǿ
end
N=255;
Br=(B/4.0)*N;
subplot(1,2,1)      %%һ�������еĵ�һ��ͼ
image(xs,ys,Br);
colormap(gray(N));  %%ʹͼ���ԻҶ�ͼ����ʾ
subplot(1,2,2)      %%һ�������еĵڶ���ͼ
plot(B,ys)