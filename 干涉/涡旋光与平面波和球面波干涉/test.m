clc
clear all
close all
%%  涡旋光束与平面光波干涉
N = 300;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-1e-2,1e-2,N);
y = linspace(-1e-2,1e-2,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);

figure;
for m = -4:4
    E1 = exp(-1i*k*X);      %平面波
    E2 = exp(1i*m*theta);   %涡旋光
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
suptitle('涡旋光束与平面光波干涉')   %为图一添加总标题
%% 涡旋光束与球面光波干涉
N = 200;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-2e-3,2e-3,N);
y = linspace(-2e-3,2e-3,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
Z = 1;
figure;
for m = -4:4
    E3 = exp(-1i*k*Z*(1+0.5*X.^2/Z^2+0.5*Y.^2/Z^2));      %球面波
    E4 = exp(1i*m*theta);   %涡旋光
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
suptitle('涡旋光束与球面光波干涉')   %为图二添加总标题