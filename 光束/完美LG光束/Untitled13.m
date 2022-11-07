%�������Ƕ���˹����
clc
clear all
close all
%% �����ѧ��λԪ���������������ⳡ
N = 300;
row = linspace(-1.5e-6,1.5e-6,N);   
col = linspace(-1.5e-6,1.5e-6,N);
[x,y] = meshgrid(row,col);
[phi,rho] = cart2pol(x,y);
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
f = 0.5;            %͸������0.5m
R = 0.5;            %Բ�׹����뾶0.3m
X = k*rho/f;
alpha = 20;          %����ʸ - ����

figure;
for m = 0 : 8
    subplot(3,3,m+1)
    E1 = (-1i)^(m+1)*(k*R/f)*exp(1i*m*phi).*...
            ((alpha*besselj(m+1,alpha*R)*besselj(m,X*R)-X.*besselj(m,alpha*R).*besselj(m+1,X*R))./(alpha^2-X.^2));  
    I1 = E1.*conj(E1);  I1 = I1/max(max(I1));
    %��ά
    h1 = pcolor(x,y,I1);
    colorbar;
    set(h1,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    %colormap(gray);        %����Ҷ�ͼ��
    axis square;
    
    %��ά
%     mesh(x,y,I1)           %��ά
%     set(gca,'fontname','times new roman','fontsize',16);
%     title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
%     xlabel('x/m','fontname','times new roman','fontsize',16);
%     ylabel('y/m','fontname','times new roman','fontsize',16);
%     zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
suptitle('�����ѧ��λԪ���������������ⳡ')   %Ϊͼһ����ܱ���
%% ����׶͸���������������ⳡ
row = linspace(-0.8,0.8,N);   
col = linspace(-0.8,0.8,N);
[x,y] = meshgrid(row,col);
[phi,rho] = cart2pol(x,y);
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
f = 0.5;            %͸������0.5m
R = 0.5;            %Բ�׹����뾶0.3m
alpha = 20;          %����ʸ - ����
w0 = 0.3;           %͸����ƽ���ϵ�����
w_g = 0.5;          %��˹�������
figure;
for m = 0 : 8
    subplot(3,3,m+1)
    E2 =  (1i)^(m+1)*(w_g/w0)*exp(1i*m*phi).*exp(-(rho.^2+R^2)/w0^2).*besseli(m,2*R*rho/w0^2);
    I2 = E2.*conj(E2);  I2 = I2/max(max(I2));
    %��ά
    h2 = pcolor(x,y,I2);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['m = ',num2str(m)]);
    %colormap(gray);        %����Ҷ�ͼ��
    axis square;
    
    %��ά
     mesh(x,y,I2)           %��ά
     set(gca,'fontname','times new roman','fontsize',16);
     title(['m = ',num2str(m)],'fontname','��������','fontsize',16);
     xlabel('x/m','fontname','times new roman','fontsize',16);
     ylabel('y/m','fontname','times new roman','fontsize',16);
     zlabel('��һ��ǿ��','fontname','��������','fontsize',16);
end
suptitle('����׶͸���������������ⳡ')   %Ϊͼ������ܱ���