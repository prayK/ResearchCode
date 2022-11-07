%����������ƽ��Ⲩ������Ⲩ����
clc
clear 
close all
%%  ����������ƽ��Ⲩ����
N = 300;            %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
x = linspace(-1e-2,1e-2,N);
y = linspace(-1e-2,1e-2,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);

figure;
for m = -4:4
    E1 = exp(-1i*k*X);      %ƽ�沨
    E2 = exp(1i*m*theta);   %������
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
suptitle('����������ƽ��Ⲩ����')   %Ϊͼһ����ܱ���

%% ��������������Ⲩ����
N = 200;            %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
x = linspace(-2e-3,2e-3,N);
y = linspace(-2e-3,2e-3,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
Z = 1;

figure;
for m = -4:4
    E3 = exp(-1i*k*Z*(1 + (0.5 * X .^ 2 / Z^2) + (0.5 * Y .^ 2 / Z^2)));%�����      %���沨
    E4 = exp(1i*m*theta);   
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
suptitle('��������������Ⲩ����')   %Ϊͼ������ܱ���

%% ����������������������
N = 300;            %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
x = linspace(-1e-3,1e-3,N);
y = linspace(-1e-3,1e-3,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
figure;

for m = -4:4
    E3 = exp(-1i*m*theta);  
    E4 = exp(1i*m*theta);   
    c2 = E3+E4;
    E_2 = c2.*conj(c2);
    subplot(3,3,m+5)
    h2 = pcolor(X,Y,E_2);
    colorbar;
    set(h2,'edgecolor','none','facecolor','interp');
    title(['  m = ',num2str(m)]);
    colormap(gray);
    axis square;
end
suptitle('����������������������')   %Ϊͼ������ܱ���