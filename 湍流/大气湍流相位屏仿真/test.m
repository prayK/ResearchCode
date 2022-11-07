%% example_ft_sh_phase_screen.m  ����г���������������λ��������
clc;
Dx = 2; % length of one side of square phase screen [m]����λ���߳�����λ��m
Dy = 2;
Nx = 256; % number of grid points per side��ȡ������Ϊ256
Ny = Nx;
N=Nx;
L0 = 10; % outer scale [m]�������⾶
l0 = 0.001;% inner scale [m]�������ھ�
Cn2 = 1e-12; % coherence diameter [m^(-2/3)].��Cn��Ϊ1e-14ʱ��1550nmʱ��r0=0.010255m��3800nmʱ��r0=0.03008m��r0�ļ������ʽ9.42
Lambda=1550*1e-9;  %��������λm
K_k = 2*pi / Lambda; % optical wavenumber [rad/m]���Ⲩ��
Delta_z=5e2;  %�ֶξ��룬��λm
subh=3;%%г������
 % SW and PW coherence diameters [m]�����沨��ƽ�沨�ĸ���뾶
 r0sw = (0.423 * K_k^2 * Cn2 * 3/8 * Delta_z)^(-3/5);   %%���沨
 r0pw = (0.423 * K_k^2 * Cn2 * Delta_z)^(-3/5);  %%ƽ�沨
 
delta_x=Dx/Nx;   %% grid step along x axis; x����������
delta_y=Dy/Ny;  %%�������
 
x = (-Nx/2 : Nx/2-1) * delta_x;  %spatial grid����ɢ���ڲ�ͬ�����㴦��Ӧ��x����ֵ
y = (-Ny/2:Ny/2-1)*delta_y;  %�������Ӧ��y����ֵ
r_val=sqrt(x((Nx/2+1):Nx).^2+y((Ny/2+1):Ny).^2);  %%ԭ��-��Ե��rֵ
[x y]=meshgrid(x, y);    %%���nx*ny�ľ���
 
 
Structure_ther=6.88*(r_val/r0pw).^(5/3);  %%���۽ṹ����ֵ
 mask = circ(x, y, 1);
% D_val_hill = structure_function_hill(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);  %%hill�׵�fft�ṹ����
% %figure,pcolor(x,y,D_val);
% figure,plot(1:(Ny/2),D_val_hill(Nx/2+1,(Ny/2+1):Ny))
% 
% D_val_hill_sub = structure_function_hill_sub(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh)    %%hill�׵�subharmonic�ṹ����
% subh1=10;
% D_val_hill_sub1 = structure_function_hill_sub(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh1)
% figure,plot(1:(Ny/2),D_val_hill_sub(Nx/2+1,(Ny/2+1):Ny),1:(Ny/2),D_val_hill_sub1(Nx/2+1,(Ny/2+1):Ny))
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�������ص�һ����������Ľ��
% wvl=Lambda;
% D=Dx;
% dz=Delta_z;
% N=Nx;
% CN=Cn2;
% 
% delta=D/N;
% x=(-N/2:N/2-1)*delta;
% y=x;
% [X Y]=meshgrid(x,y);
% del_f=1/(N*delta);
% fx=(-N/2:N/2-1)*del_f;
% [kx ky]=meshgrid(2*pi*fx);
% k=2*pi/wvl;
% [th ka]=cart2pol(kx,ky);
% km=5.92/l0;
% k0=2*pi/L0;
% % r0=0.185*(wvl^2/(dz*CN))^(3/5);
% PSD_phi=0.033*CN*exp(-(ka/km).^2)./(ka.^2+k0^2).^(11/6);
% PSD_phi(N/2+1,N/2+1)=0;
% cn=2*pi*k.^2*dz.*PSD_phi*(2*pi*del_f).^2;
% phz_hi=ifft2((randn(N)+1i*randn(N)).*sqrt(cn));%????matlab��??��??FFT???��??????????cn????????��??????????????????��del_f=1;
% phz_hi=real(phz_hi);
% % figure;imagesc(phz_hi);colorbar;
% %% ????????
% phz_lo=zeros(size(phz_hi));
% for p=1:3
%     del_fp=1/(3^p*D); 
%     fx1=(-1:1)*del_fp;
%     [kx1 ky1]=meshgrid(2*pi*fx1);
%     [th1 k1]=cart2pol(kx1,ky1);
%     km=5.92/l0;
%     k0=2*pi/L0;%outscale frequency
%     PSD_phi1=0.033*CN*exp(-(k1/km).^2)./(k1.^2+k0^2).^(11/6);
%     PSD_phi1(2,2)=0;
%     %random draws of Fourier coefficient
%     cn1=2*pi*k.^2*dz.*PSD_phi1*(2*pi*del_fp).^2;
%     cn1=(randn(3)+1i*randn(3)).*sqrt(cn1);
%     SH=zeros(N);
%     for ii=1:9
%         SH=SH+cn1(ii)*exp(1i*(kx1(ii)*X+ky1(ii)*Y));
%     end
%     phz_lo=phz_lo+SH;
% end
% phz_lo=real(phz_lo)-mean(real(phz_lo(:)));
% phz=phz_hi+phz_lo;
%  figure;imagesc(phz_hi);colorbar;
%  figure;imagesc(phz);colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
% generate a random draw of an atmospheric phase screen������λ��
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%��⽫��λ����Ƶ�����Ƿ���Ϊ0ʱ����������λ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ӳ�������λ����������������Ϊ0��ͼ����ɫ���������ڵĵط����ڣ���ͼ�ζ���һ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%������Ϊ0֮���ֵ�ر�С���ǲ��ǲ�̫�ԣ�������û�н�����Ϊ0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phz_central_zero phz_central_notzero] = ft_phase_screen_zero_or_notzero(r0, N, delta, L0, l0);  %%��⽫��Ƶ�����Ƿ����ó�Ϊ0ʱ�Ľ���Ƿ�һ��
% %[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0,subh);
% %phz=phz_hi;
% %phz=phz_lo;
% %phz = phz_lo + phz_hi;%%���յ������λ���ǵ�Ƶ���ֺ͸�Ƶ���ֵļӺ�
% figure % �½���ͼ���ڣ���������ǰͼ�Σ����»�һ���µ�
% subplot(2,1,1),pcolor(x,y,phz_central_zero);
% colormap gray
% colormap gray
% shading interp
% colorbar
% subplot(2,1,2),pcolor(x,y,phz_central_notzero);
% colormap gray
% shading interp
% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ƿ�ʹ��fftshiftʱ����������λ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ӳ�������λ��������������fftshift���Ͳ���������Ҫ��ͼ����������Ϊʲô
%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ƿ���ʹ��ifftshift��ֻ�Ƕ���������ͼ�������˶Խ��ߵ�Ų��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phz_1 phz_2] = ft_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);  
% %[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0,subh);
% %phz=phz_hi;
% %phz=phz_lo;
% %phz = phz_lo + phz_hi;%%���յ������λ���ǵ�Ƶ���ֺ͸�Ƶ���ֵļӺ�
% figure % �½���ͼ���ڣ���������ǰͼ�Σ����»�һ���µ�
% subplot(2,1,1),pcolor(x,y,phz_1);
% colormap gray
% shading interp
% colorbar
% subplot(2,1,2),pcolor(x,y,phz_2);
% colormap gray
% shading interp
% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%�Ҿ��ã����ֱ����fft�任������λ���Ļ�����Ӧ�ý�Ƶ���е�ԭ��λ�ã�Ӧ���Ǿ����е�(N/2+1,N/2+1)λ�ã���ֵ��Ϊ��
%%%��Ҫʹ��subharmonic����ʱ����Ҫ������λ����Ϊ�㣬��Ϊ��һ����ԭ��λ�ý����˸�ϸ�Ļ�����
 
 
[phz_lo phz_hi]= ft_sh_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh);%%phz_loֻ��г�������ĵ�Ƶ���֣�phz_hi�����ķ���Ϊ0ʱ��FFT
[phz_z phz_nz] = ft_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);
 phz_1=phz_hi+phz_lo; %%�ܵ�г����������λ��
 phz_2=phz_hi;  %%FFT����
 
[phz_lo_1 phz_hi_2]  = ft_sh_phase_screen(r0pw, Nx, delta_x, L0, l0,subh);
figure
pcolor(x,y,phz_lo_1+phz_hi_2);%imagesc(phz_2);
 colormap gray
 shading interp
colorbar
 
 figure % �½���ͼ���ڣ���������ǰͼ�Σ����»�һ���µ�
subplot(2,1,1),pcolor(x,y,phz_1);%imagesc(phz_1);
 colormap gray
 shading interp
colorbar
subplot(2,1,2),pcolor(x,y,phz_2);%imagesc(phz_2);
 colormap gray
 shading interp
colorbar
 
[cs1 cs2]= ft_sh_phase_screen(r0pw, Nx, delta_x, L0, l0,subh);
mask = ones(Nx);
CCC1=str_fcn2_ft(phz_hi_2, mask, delta_x);
CCC2=str_fcn2_ft(phz_lo_1+phz_hi_2, mask, delta_x);
CCC3=str_fcn2_ft(cs1+cs2, mask, delta_x);
CCC4=Structure_ther;
figure
plot(1:(Ny/2),CCC2(Nx/2+1,(Ny/2+1):Ny),'b',1:(Ny/2),CCC3(Nx/2+1,(Ny/2+1):Ny),'r',1:(Ny/2),CCC1(Nx/2+1,(Ny/2+1):Ny),'g')
%pcolor(x,y,CCC);%imagesc(phz_2);
 
 
 
 
phz = phz_lo + phz_hi;%%���յ������λ���ǵ�Ƶ���ֺ͸�Ƶ���ֵļӺ�
figure % �½���ͼ���ڣ���������ǰͼ�Σ����»�һ���µ�
pcolor(x,y,phz);
% min(phz(:))
% max(phz(:))
% imagesc(phz)
colormap gray(255)  %ͼ���ǻ�ɫ��
shading interp  %��ͼ����������ɫ.faceted:������Ϊ��ɫ��flat:�����߷ֿ���ɫ��interp����ɫ�Ĺ�˳�����
colorbar
 
 
mesh(x,y,phz);view(0,90); %%����άͼ��%��ά����ͼ��view�Ǹı��ӽ�
contourf(x,y,phz) %�ȸ���ͼ
pcolor(x,y,phz);shading interp%α��ɫͼ
surf(x,y,phz);
plot3(x,y,phz);%����ά����ͼ
stem3(x,y,phz);
contour3(x,y,phz);
imagesc(phz);