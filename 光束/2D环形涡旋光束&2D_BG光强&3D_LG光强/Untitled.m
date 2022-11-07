clc
clear all
close all
%% ������������
N = 200;
lambda = 632e-9;    %����Ϊ632nm
k = 2*pi/lambda;    %����
w0 = 3;             %�����뾶
x = linspace(-10,10,N);
y = linspace(-10,10,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);

beta = 50*pi/180;
figure;
for m = -4 : 4
    subplot(3,3,m+5)
    E1 = (r/w0).^abs(m).*exp(-r.^2/w0^2)*exp(1i*beta).*exp(-1i*m*theta);
    I1 = E1.*conj(E1);  I1 = I1/max(max(I1));
    
    %��ά
    h1 = pcolor(X,Y,I1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    axis square;
    
    %     %��ά
    %     mesh(X,Y,I1)
    %     set(gca,'fontname','times new roman','fontsize',16);
    %     title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
    %     %xlabel('x/mm','fontname','times new roman','fontsize',16);
    %     %ylabel('y/mm','fontname','times new roman','fontsize',16);
    %     %zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
suptitle('����������������ͬ���˺���(m)')   %Ϊͼһ����ܱ���
%% ������-��˹����
N = 200;
lambda = 632e-9;    %����Ϊ632nm
k = 2*pi/lambda;    %����
w0 = 3;             %�����뾶
x = linspace(-5,5,N);
y = linspace(-5,5,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
figure;
alpha = 5;
for m = -4 : 4
    subplot(3,3,m+5)
    E2 = besselj(m,alpha.*r).*exp(-r.^2/w0^2).*exp(-1i*m*theta); %ʹ����matlab���õı���������
    I2 = E2.*conj(E2);  I2 = I2/max(max(I2));
    
    %��ά
    h2 = pcolor(X,Y,I2);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    %colormap(gray);        %����Ҷ�ͼ��
    axis square;
    
    %     %��ά
    %     mesh(X,Y,I2)           %��ά
    %     set(gca,'fontname','times new roman','fontsize',16);
    %     title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
    %     %xlabel('x/mm','fontname','times new roman','fontsize',16);
    %     %ylabel('y/mm','fontname','times new roman','fontsize',16);
    %     %zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
suptitle('������-��˹��������ͬ���˺���(m)')   %Ϊͼ������ܱ���
%% ���Ƕ�-��˹����
N = 200;
lambda = 632e-9;    %����Ϊ632nm
k = 2*pi/lambda;    %����
w0 = 3e-3;          %��߳ߴ�
x = linspace(-3*w0,3*w0,N);     y = x;
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
Z_R = pi*w0^2/lambda;      %��������
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%������zλ�õİ뾶
figure;
p = 2;      %p = 0, 1, 2...;
for m = -4 : 4
    subplot(3,3,m+5)
    E3 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    I3 = E3.*conj(E3);  I3 = I3/max(max(I3));
    
%     %��ά
%     h3 = pcolor(X,Y,I3);
%     colorbar;
%     set(h3,'edgecolor','none','facecolor','interp');
%     title(['m = ',num2str(m)]);
%     %colormap(gray);        %����Ҷ�ͼ��
%     axis square;
    
        %��ά
        mesh(X,Y,I3)           %��ά
        set(gca,'fontname','times new roman','fontsize',16);
        title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
        xlabel('x/mm','fontname','times new roman','fontsize',16);
        ylabel('y/mm','fontname','times new roman','fontsize',16);
        zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
suptitle(['���Ƕ�-��˹��������ͬ���˺���(m)    p = ',num2str(p)])   %Ϊͼ������ܱ���
%% ���Ƕ�����ʽ(����5�еĹ�ʽ)
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
