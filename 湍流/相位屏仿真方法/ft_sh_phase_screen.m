function [phz_lo phz_hi] = ft_sh_phase_screen(r0,N,delta,L0,l0,subh) 
% function [phz_lo phz_hi] ... 
% = ft_sh_phase_screen(r0, N, delta, L0, l0,subh) 


D = N*delta; 
% high-frequency screen from FFT method 
phz_hi = ft_phase_screen(r0, N, delta, L0, l0); 
% spatial grid [m] 
[x y] = meshgrid((-N/2 : N/2-1) * delta); 
% initialize low-freq screen 
phz_lo = zeros(size(phz_hi)); 
% loop over frequency grids with spacing 1/(3^p*L) 
for p = 1:subh 
    % setup the PSD 
    del_f = 1 / (3^p*D); %frequency grid spacing [1/m] 
    fx = (-1 : 1) * del_f; 
    % frequency grid [1/m] 
    [fx fy] = meshgrid(fx); 
    [th f] = cart2pol(fx, fy); % polar grid 
    fm = 5.92/l0/(2*pi); % inner scale frequency [1/m] 
    f0 = 1/L0; % outer scale frequency [1/m] 
    % Kolmogorov phase PSD 
    PSD_phi = 0.023*r0^(-5/3) * f.^(-11/6);
    PSD_phi(2,2) = 0; 
    % random draws of Fourier coefficients 
    cn = (randn(3) + i*randn(3)) .* sqrt(PSD_phi)*del_f; 
    
    SH = zeros(N); 
    % loop over frequencies on this grid 
    for ii = 1:9 
        SH = SH + cn(ii)* exp(i*2*pi*(fx(ii)*x+fy(ii)*y)); 
        
    end
    phz_lo = phz_lo + SH; % accumulate subharmonics 
end
phz_lo = real(phz_lo) - mean(real(phz_lo(:)));

%  function [phz_lo phz_hi]= ft_sh_phase_screen(r0, N, delta, L0, l0,subh)
% %%subh是谐波次数
% %%这一个程序只是对应于lane方法中产生低频相位屏的部分。总的相位屏是低频和高频的总和。