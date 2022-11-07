%���Ƕ���˹������������
function main()
clc;
clear;
close all;
%% ��������
%   /* �������� */
Nxy = 200;          %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
p = 0;
z = 0;
%% ��ͬ���˺����ĵ�������ͼ��
w0 = 1e-3;          %��˹������� 1mm
d = 1e-3;           %��� 1mm
zf = 0.5;           %�������0.5m
Z_R = pi*w0^2/lambda;      %��������
w_z = w0*sqrt(1+(z/Z_R)^2);%������zλ�õİ뾶
%   /* �������� */
[x,y] = meshgrid(linspace(-1.5*w0,1.5*w0,Nxy));   %��������
[theta,r] = cart2pol(x,y);
[xf,yf] = meshgrid(linspace(-2e-3,2e-3,Nxy)); %��������
Ef = zeros(Nxy);
for m = -1:1              %���˺���-1,0,1
    tic
    %   /* ��Դ���� z=0ʱ��L-G���� */
    E0 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    I0 = E0.*conj(E0);      I0 = I0/max(max(I0));
    figure(1);subplot(1,3,m+2);
    mesh(x*1e2,y*1e2,I0);
    view(2);
    axis square;
    set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
    title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
    xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
    ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    if m == 1
        suptitle('��ͬ���˺�����Ӧ�Ĺ�Դ�ⳡ');
    end
    %   /* �������� */
    for a = 1:Nxy
        for b = 1:Nxy
            Ef(a,b) = (-1i/lambda/zf)*exp(1i*k*zf)*sum(sum(E0.*SingleSeam(E0,d,3*w0).*exp(1i*k/2/zf*((xf(a,b)-x).^2+(yf(a,b)-y).^2))));
        end
    end
    If = Ef.*conj(Ef);      If = If/max(max(If));
    figure(2);subplot(1,3,m+2);
    mesh(xf*1e2,yf*1e2,If);
    colormap gray;
    view(2);
    axis square;
    set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
    title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
    xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
    ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    if m == 1
        suptitle(['{\it\sigma}=',num2str(w0*1e3),'mm,{\ita}=',num2str(d*1e3),'mm,{\itz}= ',num2str(zf),'mʱ,����ͼ�������˺����ı仯']);
    end
    toc
end
%% ��ͬ���ʱ������ͼ��
w0 = 2e-3;          %��˹������� 2mm
m = 1;              %���˺��� 1
zf = 0.5;           %�������0.5
Z_R = pi*w0^2/lambda;      %��������
w_z = w0*sqrt(1+(z/Z_R)^2);%������zλ�õİ뾶
%   /* �������� */
[x,y] = meshgrid(linspace(-1.5*w0,1.5*w0,Nxy));   %��������
[theta,r] = cart2pol(x,y);
[xf,yf] = meshgrid(linspace(-3e-3,3e-3,Nxy)); %��������
E0 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
    .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
    .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
I0 = E0.*conj(E0);      I0 = I0/max(max(I0));
figure(3);mesh(x*1e2,y*1e2,I0);
view(2);
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
for d = 0.2e-3:0.2e-3:0.6e-3        %d = 0.2,0.4,0.6mm
    tic
    %   /* �������� */
    for a = 1:Nxy
        for b = 1:Nxy
            Ef(a,b) = (-1i/lambda/zf)*exp(1i*k*zf)*sum(sum(E0.*SingleSeam(E0,d,3*w0).*exp(1i*k/2/zf*((xf(a,b)-x).^2+(yf(a,b)-y).^2))));
        end
    end
    If = Ef.*conj(Ef);      If = If/max(max(If));
    figure(4);subplot(1,3,round(d*1e3/0.2));
    mesh(xf*1e2,yf*1e2,If);
    colormap gray;
    view(2);
    axis square;
    set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
    title(['d = ',num2str(d*1e3),'mm'],'fontname','��������','fontsize',16);
    xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
    ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    if round(d*1e3/0.2) == 3
        suptitle(['{\it\sigma}=',num2str(w0*1e3),'mm,{\itl}=',num2str(m),',{\itz}= ',num2str(zf),'mʱ,����ͼ������ı仯']);
    end
    toc
end
%% ��ͬ��߰뾶ʱ������ͼ��
d = 1e-3;           %���1mm
m = 1;              %���˺��� 1
zf = 0.5;           %�������0.5m
for w0 = 1e-3:1e-3:3e-3        %d = 1,2,3mm
    tic
    Z_R = pi*w0^2/lambda;      %��������
    w_z = w0*sqrt(1+(z/Z_R)^2);%������zλ�õİ뾶
    %   /* �������� */
    [x,y] = meshgrid(linspace(-1.5*w0,1.5*w0,Nxy));   %��������
    [theta,r] = cart2pol(x,y);
    [xf,yf] = meshgrid(linspace(-4e-3,4e-3,Nxy)); %��������
    E0 = sqrt(2*factorial(p)/pi/(p+factorial(abs(m))))*(1/w_z)*(sqrt(2)*r/w_z).^abs(m)...
        .*exp(-r.^2/w_z^2).*laguerre(p,abs(m),2*r.^2/w_z^2).*exp(-1i*m*theta).*exp(-1i*k*z)...
        .*exp(-1i*k*r.^2*z/2/(z^2+Z_R^2))*exp(-1i*(2*p+abs(m)+1)*atan(z/Z_R));
    I0 = E0.*conj(E0);      I0 = I0/max(max(I0));
    figure(5);subplot(1,3,w0*1e3);
    mesh(x*1e2,y*1e2,I0);
    view(2);
    axis square;
    set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
    title(['\it\sigma = ',num2str(w0*1e3),'mm'],'fontname','��������','fontsize',16);
    xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
    ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    if m == 3e-3
        suptitle('��ͬ��߰뾶��Ӧ�Ĺ�Դ�ⳡ');
    end
    %   /* �������� */
    for a = 1:Nxy
        for b = 1:Nxy
            Ef(a,b) = (-1i/lambda/zf)*exp(1i*k*zf)*sum(sum(E0.*SingleSeam(E0,d,3*w0).*exp(1i*k/2/zf*((xf(a,b)-x).^2+(yf(a,b)-y).^2))));
        end
    end
    If = Ef.*conj(Ef);      If = If/max(max(If));
    figure(6);subplot(1,3,w0*1e3);
    mesh(xf*1e2,yf*1e2,If);
    colormap gray;
    view(2);
    axis square;
    set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
    title(['d = ',num2str(d*1e3),'mm'],'fontname','��������','fontsize',16);
    xlabel('{\itx}/cm','fontname','times new roman','fontsize',16);
    ylabel('{\ity}/cm','fontname','times new roman','fontsize',16);
    zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
    if w0 == 3e-3
        suptitle(['{\ita}=',num2str(d*1e3),'mm,{\itl}=',num2str(m),',{\itz}= ',num2str(zf),'mʱ,����ͼ�����߰뾶�ı仯']);
    end
    toc
end

end
%% ���Ƕ�����ʽ
function result = laguerre(p,l,x)
if p == 0
    result = 1;
elseif p == 1
    result = 1+abs(l)-x;
else
    result = (1/p)*((2*p+l-1-x).*laguerre(p-1,abs(l),x)-(p+l-1)*laguerre(p-2,abs(l),x));
end
end
%% ����͸�亯��
function [result] = SingleSeam(E0,d,L)
% d,L�ֱ��Ƿ���Լ�x�����ȡֵ����
[Nx,Ny] = size(E0);
result = zeros(Nx,Ny);
X_begin = ceil((1-d/L)*Nx/2);
X_end = ceil((1+d/L)*Nx/2);
for a = X_begin:X_end
    result(:,a) = 1;
end

end