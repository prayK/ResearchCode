%% example_ft_sh_phase_screen.m  利用谐波方法产生随机相位屏的例子
clc;
Dx = 2; % length of one side of square phase screen [m]，相位屏边长，单位是m
Dy = 2;
Nx = 256; % number of grid points per side，取样点数为256
Ny = Nx;
N=Nx;
L0 = 10; % outer scale [m]，湍流外径
l0 = 0.001;% inner scale [m]，湍流内径
Cn2 = 1e-12; % coherence diameter [m^(-2/3)].当Cn设为1e-14时，1550nm时，r0=0.010255m；3800nm时，r0=0.03008m。r0的计算见公式9.42
Lambda=1550*1e-9;  %波长，单位m
K_k = 2*pi / Lambda; % optical wavenumber [rad/m]，光波数
Delta_z=5e2;  %分段距离，单位m
subh=3;%%谐波次数
 % SW and PW coherence diameters [m]，球面波和平面波的干涉半径
 r0sw = (0.423 * K_k^2 * Cn2 * 3/8 * Delta_z)^(-3/5);   %%球面波
 r0pw = (0.423 * K_k^2 * Cn2 * Delta_z)^(-3/5);  %%平面波
 
delta_x=Dx/Nx;   %% grid step along x axis; x方向网格间隔
delta_y=Dy/Ny;  %%采样间隔
 
x = (-Nx/2 : Nx/2-1) * delta_x;  %spatial grid，离散后在不同采样点处对应的x坐标值
y = (-Ny/2:Ny/2-1)*delta_y;  %采样点对应的y坐标值
r_val=sqrt(x((Nx/2+1):Nx).^2+y((Ny/2+1):Ny).^2);  %%原点-边缘的r值
[x y]=meshgrid(x, y);    %%变成nx*ny的矩阵
 
 
Structure_ther=6.88*(r_val/r0pw).^(5/3);  %%理论结构方程值
 mask = circ(x, y, 1);
% D_val_hill = structure_function_hill(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);  %%hill谱的fft结构方程
% %figure,pcolor(x,y,D_val);
% figure,plot(1:(Ny/2),D_val_hill(Nx/2+1,(Ny/2+1):Ny))
% 
% D_val_hill_sub = structure_function_hill_sub(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh)    %%hill谱的subharmonic结构方程
% subh1=10;
% D_val_hill_sub1 = structure_function_hill_sub(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh1)
% figure,plot(1:(Ny/2),D_val_hill_sub(Nx/2+1,(Ny/2+1):Ny),1:(Ny/2),D_val_hill_sub1(Nx/2+1,(Ny/2+1):Ny))
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%网上下载的一个程序里面的结果
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
% phz_hi=ifft2((randn(N)+1i*randn(N)).*sqrt(cn));%????matlab±??í??FFT???¨??????????cn????????±??????????????????±del_f=1;
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
 
 
 
% generate a random draw of an atmospheric phase screen生成相位屏
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%检测将相位屏零频分量是否设为0时，产生的相位屏的区别
%%%%%%%%%%%%%%%%%%%%%%%%%%%从产生的相位屏上来看，不设置为0的图形颜色更鲜明，黑的地方更黑，但图形都是一样的
%%%%%%%%%%%%%%%%%%%%%%%%%%%不过设为0之后的值特别小，是不是不太对？所以我没有将其设为0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phz_central_zero phz_central_notzero] = ft_phase_screen_zero_or_notzero(r0, N, delta, L0, l0);  %%检测将零频分量是否设置成为0时的结果是否一样
% %[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0,subh);
% %phz=phz_hi;
% %phz=phz_lo;
% %phz = phz_lo + phz_hi;%%最终的随机相位屏是低频部分和高频部分的加和
% figure % 新建画图窗口，即保存以前图形，重新画一个新的
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%检测是否使用fftshift时，产生的相位屏的区别
%%%%%%%%%%%%%%%%%%%%%%%%%%%从产生的相位屏上来看，先用fftshift，就产生不出想要的图样，不明白为什么
%%%%%%%%%%%%%%%%%%%%%%%%%%%最后是否再使用ifftshift，只是对最后产生的图样进行了对角线的挪移
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phz_1 phz_2] = ft_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);  
% %[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0,subh);
% %phz=phz_hi;
% %phz=phz_lo;
% %phz = phz_lo + phz_hi;%%最终的随机相位屏是低频部分和高频部分的加和
% figure % 新建画图窗口，即保存以前图形，重新画一个新的
% subplot(2,1,1),pcolor(x,y,phz_1);
% colormap gray
% shading interp
% colorbar
% subplot(2,1,2),pcolor(x,y,phz_2);
% colormap gray
% shading interp
% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%我觉得，如果直接用fft变换产生相位屏的话，不应该将频谱中的原点位置（应该是矩阵中的(N/2+1,N/2+1)位置）的值设为零
%%%而要使用subharmonic方法时，需要将中心位置设为零，因为进一步对原点位置进行了更细的划分了
 
 
[phz_lo phz_hi]= ft_sh_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z, subh);%%phz_lo只是谐波产生的低频部分，phz_hi是中心分量为0时的FFT
[phz_z phz_nz] = ft_phase_screen_use_cn2(Cn2, Nx, Ny, Dx, Dy, Lambda, L0, l0, Delta_z);
 phz_1=phz_hi+phz_lo; %%总的谐波产生的相位屏
 phz_2=phz_hi;  %%FFT方法
 
[phz_lo_1 phz_hi_2]  = ft_sh_phase_screen(r0pw, Nx, delta_x, L0, l0,subh);
figure
pcolor(x,y,phz_lo_1+phz_hi_2);%imagesc(phz_2);
 colormap gray
 shading interp
colorbar
 
 figure % 新建画图窗口，即保存以前图形，重新画一个新的
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
 
 
 
 
phz = phz_lo + phz_hi;%%最终的随机相位屏是低频部分和高频部分的加和
figure % 新建画图窗口，即保存以前图形，重新画一个新的
pcolor(x,y,phz);
% min(phz(:))
% max(phz(:))
% imagesc(phz)
colormap gray(255)  %图形是灰色的
shading interp  %对图形网格线着色.faceted:网格线为黑色；flat:网格线分块着色；interp：着色的光顺性最好
colorbar
 
 
mesh(x,y,phz);view(0,90); %%画三维图，%三维曲面图，view是改变视角
contourf(x,y,phz) %等高线图
pcolor(x,y,phz);shading interp%伪彩色图
surf(x,y,phz);
plot3(x,y,phz);%画三维曲线图
stem3(x,y,phz);
contour3(x,y,phz);
imagesc(phz);