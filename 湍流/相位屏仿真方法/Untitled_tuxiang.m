clc
clear;
lemda=1.064e-6;
Cn2=2e-12;
delta_z=0.3;
r0=r0(lemda,Cn2,delta_z)
D =2; % length of one side of square phase screen [m] 
%r0 = 0.1; % coherence diameter [m] 
N = 512; % number of grid points per side 
L0 = inf; % outer scale [m] 
%L0 = 100; % outer scale [m] 
l0 = 0.01;% inner scale [m] 
k0=2*pi/L0;
subh=2;% ����г�������������FFT�����ﲻ������
n=1;%����n����λ����Ȼ��ȡ��λ����ƽ������Ϊ������λ��
delta = D/N; 
x=(-N/2:N/2-1)*delta;
y=x;
phz=zeros(N,N);
for idxscr = 1 : 1 : n  
%          [phz_lo phz_hi]=ft_sh_phase_screen(r0,N,delta,L0,l0,subh);
%           phz = phz+phz_lo + phz_hi;
          phz =ft_phase_screen(r0, N, delta, L0, l0);%�����FFT�ģ��Ͱ������У�������һ��
     end

phz=phz/n;
figure(1)
pcolor(x,y,phz)%����λ���Ķ�άͼ��
axis([-1 1 -1 1])
xlabel('x[m]','FontSize',12)
ylabel('y[m]','FontSize',12)
colormap gray
shading interp
colorbar
set(gca,'FontSize',12);
title('����ǿ�ȣ�Cn2=2e-12','FontSize',14)

 figure(2)
surf(x,y,phz)%����λ������άάͼ��
%pcolor(x,y,phz)
xlabel('x[m]','FontSize',12)
ylabel('y[m]','FontSize',12)
zlabel('phz[rad]','FontSize',12)
colormap jet
shading interp
grid on
set(gca,'FontSize',12);
title('����ǿ�ȣ�Cn2=2e-12','fontsize',14)








