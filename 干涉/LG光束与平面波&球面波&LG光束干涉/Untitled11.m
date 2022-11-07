%拉盖尔高斯光束的干涉
clc
clear all
close all
%%  拉盖尔-高斯光束与平面光波干涉
N = 300;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
w0 = 3;
x = linspace(-1e-4,1e-4,N); y = x;
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);

p = 0; 
Z_R = pi*w0^2/lambda;      %瑞利长度
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%光束在z位置的半径
figure;
for m = -4:4
    E1 = exp(-1i*k*X);      %平面波
    E2 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    c1 = E1+E2;
    E_1 = c1.*conj(c1);
    subplot(3,3,m+5)
    h1 = pcolor(X,Y,E_1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    colormap(gray);
    axis square;
end
suptitle('拉盖尔-高斯光束与平面光波干涉')   %为图一添加总标题
%% 涡旋光束与球面光波干涉
N = 200;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-2e-3,2e-3,N);  y = x;
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
Z = 1;
Z_R = pi*w0^2/lambda;      %瑞利长度
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%光束在z位置的半径
figure;
for m = -4:4
    E3 = exp(-1i*k*Z*(1+0.5*X.^2/Z^2+0.5*Y.^2/Z^2));      %球面波
    E4 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    c2 = E3+E4;
    E_2 = c2.*conj(c2);
    subplot(3,3,m+5)
    h2 = pcolor(X,Y,E_2);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    colormap(gray);
    axis square;
end
suptitle('拉盖尔高斯光束与球面光波干涉')   %为图二添加总标题
%%  拓扑荷数相反的拉盖尔高斯光束互相干涉
N = 200;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-8,8,N); y = x;
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
Z_R = pi*w0^2/lambda;      %瑞利长度
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%光束在z位置的半径
figure;
for m = 1:9
    E5 = sqrt(2*factorial(p)/pi/(p+factorial(abs(-m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(-m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(-m),2*r.^2/w_z^2).*exp(-1i*(-m)*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(-m)+1)*atan(z/Z_R));
    E6 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    c2 = E5+E6;
    E_3 = c2.*conj(c2);
    subplot(3,3,m)
    h3 = pcolor(X,Y,E_3);
    colorbar;
    set(h3,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    colormap(gray);
    axis square;
end
suptitle('拓扑荷数相反的拉盖尔高斯光束互相干涉')   %为图三添加总标题
%% 拉盖尔多项式(文献3中的公式)
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