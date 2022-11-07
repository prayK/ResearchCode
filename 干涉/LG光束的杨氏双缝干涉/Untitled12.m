%拉盖尔高斯光束杨氏双缝干涉
clc
clear all
close all
%%  L-G光束双缝干涉
N = 300;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-2e-5,2e-5,N);
y = linspace(-2e-5,2e-5,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
w0 = 3e-1;          %束腰
d = 2e-4;           %双缝间距200μm
D = 9e-4;           %双缝与观察屏之间的距离900μm
p = 1;
Z_R = pi*w0^2/lambda;      %瑞利长度
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%光束在z位置的半径
figure;
for m = -4:4
    E = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    I = E.*conj(E);
    I_1 = 4*I.*cos(pi*X*d/lambda/D+delta_phi(m,Y)/2);
    subplot(3,3,m+5)
    h1 = pcolor(X,Y,I_1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    colormap(gray);
    axis square;
end
suptitle('拉盖尔-高斯光束双缝干涉')   %为图一添加总标题
%% delta_phi
function result = delta_phi(m,y)
    %delphs=@(L,y)L.*(2.*pi-2.*atan(a./y)).*(y>0&y<8)+L.*(2.*atan(-a./y)).*(y>-8&y<0);
    %m为拓扑荷数，改变拓扑荷数会使图形中的明暗条纹分布发生极大的改变
    result = m*2*(0.5*pi+atan(1e7*y));
end
%% 拉盖尔多项式(文献5中的公式)
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