clc;
close all;
clear all;
warning off;

 
%����
Scales = 201;
x      = linspace(-0.008,0.008,Scales );
y      = linspace(-0.008,0.008,Scales );
[X,Y]  = meshgrid(x,y);
rho    = sqrt(X.*X+Y.*Y);
Fais   = mod(atan2(Y,X)+pi,2*pi);

%��Ų���Դ����
lamda = 1550*10^(-9);
z     = 1;
w0    = 0.001;
k     = 2*pi/lamda;
p     = 0;
l1    = 1/2; 
%�����������
zR    = k*w0^2/2;
w     = w0*(1+(z/zR)^2)^(1/2);
t     = rho./w;
d     = 0;
alaf  = 1/w0^2-1i*k/2/z;

for m=1:length(x)
    for n=1:length(y)
        c       = sqrt(1/pi/abs(l1))*((sqrt(2))^abs(l1));
        c1(m,n) = c*(exp(1i*k*(z+eps))/(1i*lamda*(z+eps)*(w0^l1)))*(exp(1i*k*(x(m).^2+y(n).^2)/2/(z+eps)))*exp(-d^2/w0^2);
        c2(m,n) = c1(m,n)*exp(-1i*d*x(m)*k/((w0^2)*(z+eps)*alaf))*exp(d^2/(w0^4)/alaf)*(exp(-(k^2)*(x(m).^2+y(n).^2)/4/alaf/((z+eps)^2)));
        %�ֲ�����
        E(m,n)  =(c2(m,n)*pi*((-1i)^l1)*(k^l1)/alaf/((2*alaf*(z+eps))^l1))*((x(m)-d)-1i*y(n))^l1;
        I(m,n)  = abs(E(m,n)).^2;
    end
end

%������Ų���ǿ�ֲ�
ph=angle(E);
 
figure
subplot(221);
surf(x,y,ph);
colormap hot;
xlabel('x/m');
ylabel('y/m');
title('��λ')
shading interp;
view(2);
axis square

subplot(222);
surf(x,y,I);
colormap hot;
view(2);
xlabel('x/m');
ylabel('y/m');
title('ǿ��')
shading interp;
axis square

subplot(223);
dete_I=I((Scales-1)/2+1,:);
plot(x,dete_I);
xlabel('x/m');
ylabel('ǿ��/a.u.');
title('ǿ�ȷֲ�')
axis square 


N = Scales-1;  
s = zeros(1,N);
r = 50;   %�����뾶
p =(2*pi-0.0001)/N;

for th=0.0001:p:2*pi;
    e    = round(th/p+1);
    s(e) = ph(round(r*cos(th))+2*r,round(r*sin(th))+2*r)/pi*180;
end

th  = 0.0001:p:2*pi;
th1 = th/pi*180;

subplot(224);
plot(th1,s);
axis([0 360 -100 250]);
axis square
hold off
xlabel('Detecting plane Angle');
ylabel('Phase');
title('��λ�ֲ�')
