%倍频过程中复合涡旋拓扑荷数的倍频效应
clc
clear all
close all
%%
N = 300;            %取样点数
lambda = 632e-9;    %波长632nm
k = 2*pi/lambda;    %波数
x = linspace(-4,4,N);
y = linspace(-4,4,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
w0 = 3;
p = 0;              %径向量子数
%%  m_A = 0, m_B = 1
m_A = 0;            %两光束的拓扑荷数
m_B = 1;
figure;             %光强部分
F = 1;              %subplot画图
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c1 = E_A+E_B;
    E_1 = c1.*conj(c1);
    subplot(3,3,F)
    h1 = pcolor(X,Y,E_1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['相位差 delta = ',num2str(F-1),'*π/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' 理论模拟得到的倍频光强图'])   %添加总标题

figure;             %相位部分
F = 1;              %subplot画图
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c2 = E_A+E_B;
    E_2_angle = angle(c2).*conj(angle(c2));
    subplot(3,3,F)
    h2 = pcolor(X,Y,E_2_angle);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['相位差 delta = ',num2str(F-1),'*π/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' 理论模拟得到的倍频相位图'])   %添加总标题
%%  m_A = 0, m_B = 2
m_A = -1;            %两光束的拓扑荷数
m_B = 2;
figure;             %光强部分
F = 1;              %subplot画图
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c3 = E_A+E_B;
    E_3 = c3.*conj(c3);
    subplot(3,3,F)
    h3 = pcolor(X,Y,E_3);
    colorbar;
    set(h3,'edgecolor','none','facecolor','interp');
    title(['相位差 delta = ',num2str(F-1),'*π/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' 理论模拟得到的倍频光强图'])   %添加总标题

figure;             %相位部分
F = 1;              %subplot画图
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c4 = E_A+E_B;
    E_4_angle = angle(c4).*conj(angle(c4));
    subplot(3,3,F)
    h4 = pcolor(X,Y,E_4_angle);
    colorbar;
    set(h4,'edgecolor','none','facecolor','interp');
    title(['相位差 delta = ',num2str(F-1),'*π/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' 理论模拟得到的倍频相位图'])   %添加总标题