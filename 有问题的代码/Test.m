
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