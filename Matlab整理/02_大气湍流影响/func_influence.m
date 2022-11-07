function phi=func_influence(lamda,z,Cn2)

N      = 300; 
k      = 2*pi/lamda;%���������Ͳ���
Numz   = 20;%�Ѵ�������ֳ�Numz��
dz     = z/Numz;%ÿ�εľ���
L      = 70e-2;%���ڿ��
m1     = 1:N;
df     = 1/L;
%�ռ�Ƶ��
[fx,fy]= meshgrid(m1*df);
fr     = sqrt(fx.^2+fy.^2);
%�ռ䲨��
kx     = 2*pi*fx;
ky     = 2*pi*fy;
kr     = 2*pi*fr;
 
L0     = 20; 
l0     = 5e-3; 
kl     = 3.3/l0;
 
r0     = 0.185*(lamda^2/(dz*Cn2))^(3/5);
phi0   = 2*pi/L*0.0241*r0^(-5/6)*(fr.^2+1/L0^2).^(-11/12);
Gau    =(randn(N,N)+sqrt(-1)*randn(N,N))/sqrt(2); 
 
phi    =(fft2(Gau.*phi0));
phi    = abs(phi);
end 