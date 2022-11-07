clc;    clear;  close all;
%% ��������
Nxy = 200;                      %������
lambda = 1064e-9;               %����1064nm
k = 2*pi/lambda;                %����
w = 3e-3;                       %��˹�����������
[x,y] = meshgrid(linspace(-1.5*w,1.5*w,Nxy));
[theta,r] = cart2pol(x,y);
k_r = 7000;                     %����ʸ - ����
k_z = sqrt(k^2-k_r^2);          %����ʸ
z = 0;
c = 3e8;
T = lambda/c;   %����
m = 0;          %���˺���
omega = 2*pi*c/lambda;
for jt = 1:9
    theta = imrotate(theta, 90*(jt-1));      %ÿ��������λ��ת90�� ����������ת����Ƕȵ� ��ϧ��λ��ת֮�����ά�Ȳ�ƥ����
    E = besselj(m,k_r*r).*exp(-r.^2/w^2).*exp(1i*k_z*z).*exp(1i*m*theta).*exp(1i*omega*T*(jt-1));
    I = E.*conj(E);     I = I/max(max(I));
    m = m + 1;          %��Ϊ������λ����ת�ԳƵ� ��λ��תҲ�������� ���Ըı����˺������Ÿ�ֱ��һЩ
    figure;mesh(x*1e3,y*1e3,I);
    view(2);
    set(gca,'fontname','times new roman','fontsize',16);
    title([num2str(jt-1),'Tʱ�Ĺ�ǿ�ֲ�'],'fontname','��������','fontsize',16);
    xlabel('x/mm','fontname','times new roman','fontsize',16);
    ylabel('y/mm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    %set( gca, 'XTick', [], 'YTick', [] );
    F1 = getframe;
    %imwrite(F1.cdata,[ 'C:\Users\Administrator\Desktop\������\t = ',num2str(jt-1),'T.png']);  %����ͼƬ
end
