%���Ƕ���˹��������˫�����
clc
clear all
close all
%%  L-G����˫�����
N = 300;            %ȡ������
lambda = 632e-9;    %����632nm
k = 2*pi/lambda;    %����
x = linspace(-2e-5,2e-5,N);
y = linspace(-2e-5,2e-5,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
w0 = 3e-1;          %����
d = 2e-4;           %˫����200��m
D = 9e-4;           %˫����۲���֮��ľ���900��m
p = 1;
Z_R = pi*w0^2/lambda;      %��������
z = 0;
w_z = w0*sqrt(1+(z/Z_R)^2);%������zλ�õİ뾶
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
suptitle('���Ƕ�-��˹����˫�����')   %Ϊͼһ����ܱ���
%% delta_phi
function result = delta_phi(m,y)
    %delphs=@(L,y)L.*(2.*pi-2.*atan(a./y)).*(y>0&y<8)+L.*(2.*atan(-a./y)).*(y>-8&y<0);
    %mΪ���˺������ı����˺�����ʹͼ���е��������Ʒֲ���������ĸı�
    result = m*2*(0.5*pi+atan(1e7*y));
end
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