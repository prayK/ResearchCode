
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