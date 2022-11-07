clc;
close all;
clear all;
warning off;

%%
%�����������������Է���
Dist = [0:2000:20000];

figure;
indx = 0;
rad1 = [];
for L=Dist
    indx      = indx + 1;
    Cn        = 8*10^(-9); 
    th        = 4.03*Cn.^(6/5)*0.6328.^(-1/5)*L.^(3/5); 
    rad1(indx)= th*L; 
end 
plot(Dist,rad1,'b-o')  
% title('����������߳߶Ⱥʹ������Ĺ�ϵ') 
xlabel('�������') 
ylabel('����������߳߶�') 
 
figure;
indx = 0;
rad2 = [];
Cset = [0:20*10^(-9):500*10^(-9)];
for Cn = Cset
    indx       = indx + 1;
    L          = 100; 
    th         = 4.03*Cn.^(6/5)*0.6328.^(-1/5)*L.^(3/5); 
    rad2(indx) = th*L; 
end 
plot(Cset,rad2,'b-o') 
% title('����������߳߶Ⱥ�����ǿ�ȵĹ�ϵ') 
xlabel('����ǿ��') 
ylabel('����������߳߶�') 
 

figure;
indx = 0;
alf  = [];
Dist = [0:2000:20000];
for L= Dist 
    indx      = indx + 1;
    Cn        = 8*10^(-9); 
    alf(indx) = sqrt(1.75*Cn*Cn*L*3.2^(-1/3)*10^(-18)); 
end 
plot(Dist,alf,'b-o') 
% title('�����������ƫ�ƽǶȺʹ������Ĺ�ϵ') 
xlabel('�������') 
ylabel('�����������ƫ�ƽǶ�') 

figure;
indx = 0;
Cset = [0:20*10^(-9):500*10^(-9)];
alf  = [];
for Cn=Cset 
    indx      = indx + 1;
    L         = 1000; 
    alf(indx) = sqrt(1.75*Cn*Cn*L*3.2^(-1/3)*10^(-18)); 
end 
plot(Cset,alf,'b-o') 
xlabel('����ǿ��') 
ylabel('�����������ƫ�ƽǶ�') 

figure;
indx = 0;
B    = 0.49;
Dist = [0:500:20000];
I    = [];
for  L=Dist 
     indx   = indx + 1;
     Cn     = 8*10^(-9); 
     I(indx)= B*(2*pi/0.6328).^(7/6)*L.^(11)*Cn.^2*10^(-18); 
end 
plot(Dist,I,'b-o') 
xlabel('�������') 
ylabel('����������ǿ���')

figure;
indx = 0;
B    = 0.49; 
Cset = [0:20*10^(-9):500*10^(-9)];
I    = [];
for Cn=Cset
    indx   = indx + 1;
    L      = 1000; 
    I(indx)= B*(2*pi/0.6328).^(7/6)*L.^(11)*Cn.^2*10^(-18); 
end 
plot(Cset,I,'b-o') 
xlabel('����ǿ��') 
ylabel('����������ǿ���') 

%%
%����ͨ�����������ı任�������
N         = 300;
w0        = 0.03;
s         = 3;
z         = 1000;
lamda     = 1.550e-6;
k         = 2*pi/lamda;
z0        = k*w0^2/2;
%ͼ�����귶Χ
b         = 0.2;
dx        = b/N;
Cn2       = 1e-15;
%�������
Numz      = 20;%�Ѵ�������ֳ�Numz��
dz        = z/Numz;%ÿ�εľ���
m2        = [-N/2:N/2-1];
[r,theta] = meshgrid(linspace(0,b,N),linspace(0,2*pi,N));
[x,y]     = pol2cart(theta,r);
%���ڿ��
L         = 0.7;
df        = 1/L;
%�ռ�Ƶ��
[fx,fy]   = meshgrid(m2*df);
fr        = sqrt(fx.^2+fy.^2);
kx        = 2*pi*fx;
ky        = 2*pi*fy;
kr        = 2*pi*fr;
A         = sqrt(2/(pi*gamma(abs(s)+1)));
%��������ı��ʽ
u0        = A*exp(-r.^2/w0^2).*(sqrt(2)*r/w0).^s.*exp(1i*s*theta);
u         = u0;
phi       = func_influence(lamda,z,Cn2);
for j=1:Numz
      u1 = fft2(exp(1i*phi).*u);
      u1 = fftshift(u1);
      u  = ifft2(ifftshift(exp(1i*k*dz)*exp(-1i*dz*kr.^2/(2*k)).*u1));
end
 
%ǿ�Ⱥ���λͼ
[X,Y] = meshgrid(m2*dx);
uz    = griddata(x,y,u,X,Y);

figure
subplot(121);
imshow(abs(uz).^2,[]);
axis square;
xlabel('x/mm');ylabel('y/mm')
subplot(122);
imshow(angle(uz),[]);
axis square;
xlabel('x/mm');ylabel('y/mm')

