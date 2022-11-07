%��Ƶ�����и����������˺����ı�ƵЧӦ
clc
clear all
close all
%%
N = 300;            %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
x = linspace(-4,4,N);
y = linspace(-4,4,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
w0 = 3;
p = 0;              %����������
%%  m_A = 0, m_B = 1
m_A = 0;            %�����������˺���
m_B = 1;
figure;             %��ǿ����
F = 1;              %subplot��ͼ
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c1 = E_A+E_B;
    E_1 = c1.*conj(c1);
    subplot(3,3,F)
    h1 = pcolor(X,Y,E_1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['��λ�� delta = ',num2str(F-1),'*��/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' ����ģ��õ��ı�Ƶ��ǿͼ'])   %����ܱ���

figure;             %��λ����
F = 1;              %subplot��ͼ
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c2 = E_A+E_B;
    E_2_angle = angle(c2).*conj(angle(c2));
    subplot(3,3,F)
    h2 = pcolor(X,Y,E_2_angle);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['��λ�� delta = ',num2str(F-1),'*��/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' ����ģ��õ��ı�Ƶ��λͼ'])   %����ܱ���
%%  m_A = 0, m_B = 2
m_A = -1;            %�����������˺���
m_B = 2;
figure;             %��ǿ����
F = 1;              %subplot��ͼ
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c3 = E_A+E_B;
    E_3 = c3.*conj(c3);
    subplot(3,3,F)
    h3 = pcolor(X,Y,E_3);
    colorbar;
    set(h3,'edgecolor','none','facecolor','interp');
    title(['��λ�� delta = ',num2str(F-1),'*��/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' ����ģ��õ��ı�Ƶ��ǿͼ'])   %����ܱ���

figure;             %��λ����
F = 1;              %subplot��ͼ
for delta = 0:pi/4:2*pi
    E_A = sqrt(2/pi/factorial(abs(m_A)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_A).*exp(-r.^2/w0).*exp(-1i*m_A*theta);
    E_B = sqrt(2/pi/factorial(abs(m_B)))*(1/w0)*(sqrt(2)*r/w0).^abs(m_B).*exp(-r.^2/w0).*exp(-1i*m_B*theta)*exp(1i*delta);
    c4 = E_A+E_B;
    E_4_angle = angle(c4).*conj(angle(c4));
    subplot(3,3,F)
    h4 = pcolor(X,Y,E_4_angle);
    colorbar;
    set(h4,'edgecolor','none','facecolor','interp');
    title(['��λ�� delta = ',num2str(F-1),'*��/4']);
    %colormap(gray);
    axis square;
    F = F+1;
end
suptitle(['m_A=',num2str(m_A),' m_b=',num2str(m_B), ' ����ģ��õ��ı�Ƶ��λͼ'])   %����ܱ���