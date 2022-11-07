clc;
close all;
clear all;
warning off;


%%
%光线通过大气湍流的变换仿真分析
N         = 300;
w0        = 0.03;
s         = 5;
z         = 1000;
lamda     = 1.550e-6;
k         = 2*pi/lamda;
z0        = k*w0^2/2;
%图像坐标范围
b         = 0.2;
dx        = b/N;
Cn2       = 1e-13;
%传输距离
Numz      = 10;%把传播距离分成Numz段
dz        = z/Numz;%每段的距离
m2        = [-N/2:N/2-1];
[r,theta] = meshgrid(linspace(0,b,N),linspace(0,2*pi,N));
[x,y]     = pol2cart(theta,r);
%窗口宽度
L         = 0.7;
df        = 1/L;
%空间频率
[fx,fy]   = meshgrid(m2*df);
fr        = sqrt(fx.^2+fy.^2);
kx        = 2*pi*fx;
ky        = 2*pi*fy;
kr        = 2*pi*fr;
A         = sqrt(2/(pi*gamma(abs(s)+1)));
%发射光束的表达式
u0        = A*exp(-r.^2/w0^2).*(sqrt(2)*r/w0).^s.*exp(1i*s*theta);
u         = u0;
us        = u;
phi       = func_influence(lamda,z,Cn2);
for j=1:Numz
      u1 = fft2(exp(1i*phi).*u);
      u1 = fftshift(u1);
      u  = ifft2(ifftshift(exp(1i*k*dz)*exp(-1i*dz*kr.^2/(2*k)).*u1));
end
 
%强度和相位图
[X,Y] = meshgrid(m2*dx);
uz    = griddata(x,y,u,X,Y);

figure
subplot(121);
imshow(abs(uz).^2,[]);
axis square;
xlabel('x/mm');ylabel('y/mm')
title('通过大气湍流');
subplot(122);
imshow(angle(uz),[]);
axis square;
xlabel('x/mm');ylabel('y/mm')

%%
%校正
phi0      = -phi;
for j=1:Numz
      u1 = fft2(exp(1i*phi0).*u);
      u1 = fftshift(u1);
      u  = ifft2(ifftshift(exp(-1i*k*dz)*exp(1i*dz*kr.^2/(2*k)).*u1));
end
%强度和相位图
[X,Y] = meshgrid(m2*dx);
uz    = griddata(x,y,u,X,Y);
figure
subplot(121);
imshow(abs(uz).^2,[]);
axis square;
xlabel('x/mm');ylabel('y/mm')
title('校正后效果');
subplot(122);
imshow(angle(uz),[]);
axis square;
xlabel('x/mm');ylabel('y/mm')



% 
% %%
% %光线通过大气湍流的变换仿真分析
% N         = 300;
% w0        = 0.03;
% s         = 5;
% z         = 1000;
% lamda     = 1.550e-6;
% k         = 2*pi/lamda;
% z0        = k*w0^2/2;
% %图像坐标范围
% b         = 0.2;
% dx        = b/N;
% Cn2       = 1e-15;
% %传输距离
% Numz      = 40;%把传播距离分成Numz段
% dz        = z/Numz;%每段的距离
% m2        = [-N/2:N/2-1];
% [r,theta] = meshgrid(linspace(0,b,N),linspace(0,2*pi,N));
% [x,y]     = pol2cart(theta,r);
% %窗口宽度
% L         = 0.7;
% df        = 1/L;
% %空间频率
% [fx,fy]   = meshgrid(m2*df);
% fr        = sqrt(fx.^2+fy.^2);
% kx        = 2*pi*fx;
% ky        = 2*pi*fy;
% kr        = 2*pi*fr;
% A         = sqrt(2/(pi*gamma(abs(s)+1)));
% %发射光束的表达式
% u0        = A*exp(-r.^2/w0^2).*(sqrt(2)*r/w0).^s.*exp(1i*s*theta);
% u         = u0;
% us        = u;
% phi       = func_influence(lamda,z,Cn2);
% for j=1:Numz
%       u1 = fft2(exp(1i*phi).*u);
%       u1 = fftshift(u1);
%       u  = ifft2(ifftshift(exp(1i*k*dz)*exp(-1i*dz*kr.^2/(2*k)).*u1));
% end
%  
% %强度和相位图
% [X,Y] = meshgrid(m2*dx);
% uz    = griddata(x,y,u,X,Y);
% 
% figure
% subplot(221);
% imshow(abs(uz).^2,[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% title('通过大气湍流');
% subplot(222);
% imshow(angle(uz),[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% 
% %%
% %校正
% phi0      = -phi;
% for j=1:Numz
%       u1 = fft2(exp(1i*phi0).*u);
%       u1 = fftshift(u1);
%       u  = ifft2(ifftshift(exp(-1i*k*dz)*exp(1i*dz*kr.^2/(2*k)).*u1));
% end
% %强度和相位图
% [X,Y] = meshgrid(m2*dx);
% uz    = griddata(x,y,u,X,Y);
% subplot(223);
% imshow(abs(uz).^2,[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% title('校正后效果');
% subplot(224);
% imshow(angle(uz),[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% 
% 
% 
% 
% %%
% %光线通过大气湍流的变换仿真分析
% N         = 300;
% w0        = 0.03;
% s         = 5;
% z         = 1000;
% lamda     = 1.550e-6;
% k         = 2*pi/lamda;
% z0        = k*w0^2/2;
% %图像坐标范围
% b         = 0.2;
% dx        = b/N;
% Cn2       = 1e-15;
% %传输距离
% Numz      = 100;%把传播距离分成Numz段
% dz        = z/Numz;%每段的距离
% m2        = [-N/2:N/2-1];
% [r,theta] = meshgrid(linspace(0,b,N),linspace(0,2*pi,N));
% [x,y]     = pol2cart(theta,r);
% %窗口宽度
% L         = 0.7;
% df        = 1/L;
% %空间频率
% [fx,fy]   = meshgrid(m2*df);
% fr        = sqrt(fx.^2+fy.^2);
% kx        = 2*pi*fx;
% ky        = 2*pi*fy;
% kr        = 2*pi*fr;
% A         = sqrt(2/(pi*gamma(abs(s)+1)));
% %发射光束的表达式
% u0        = A*exp(-r.^2/w0^2).*(sqrt(2)*r/w0).^s.*exp(1i*s*theta);
% u         = u0;
% us        = u;
% phi       = func_influence(lamda,z,Cn2);
% for j=1:Numz
%       u1 = fft2(exp(1i*phi).*u);
%       u1 = fftshift(u1);
%       u  = ifft2(ifftshift(exp(1i*k*dz)*exp(-1i*dz*kr.^2/(2*k)).*u1));
% end
%  
% %强度和相位图
% [X,Y] = meshgrid(m2*dx);
% uz    = griddata(x,y,u,X,Y);
% 
% figure
% subplot(221);
% imshow(abs(uz).^2,[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% title('通过大气湍流');
% subplot(222);
% imshow(angle(uz),[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% 
% %%
% %校正
% phi0      = -phi;
% for j=1:Numz
%       u1 = fft2(exp(1i*phi0).*u);
%       u1 = fftshift(u1);
%       u  = ifft2(ifftshift(exp(-1i*k*dz)*exp(1i*dz*kr.^2/(2*k)).*u1));
% end
% %强度和相位图
% [X,Y] = meshgrid(m2*dx);
% uz    = griddata(x,y,u,X,Y);
% subplot(223);
% imshow(abs(uz).^2,[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% title('校正后效果');
% subplot(224);
% imshow(angle(uz),[]);
% axis square;
% xlabel('x/mm');ylabel('y/mm')
% 
% 
% 
% 
% 
