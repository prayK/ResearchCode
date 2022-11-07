% Plots intensities in beams created with NEC
% Input files should be with 5 degrees angular. Please use 'L(l)_5.dat'
% resolution and theta =< 90
% degrees. Input file should have columns
% (theta, phi, vert.gain, hor.gain, tot.gain, AR, tilt angle, Ethetamagn, Ethetaphase, Ephimagn, Ephiphase)

clear;
c = 3e+8;
e0 = 8.854e-12;
% number of plots
korningar = input('number of different plots: ');
T40 = 0;      % needed?
RES = 0;      % needed
for korn = 1:korningar
    clear T40;
    clear RES;

    file = input('filnamn: ','s');     %inputfile
    RES = load(file);

    R = 1;      %may change
    T40 = RES;

    % change to radians and fix neg.
    theta = pi/180*-1*T40(:,1);
    phi = pi/180*T40(:,2);

    % extract E magnitudes and phases
    Ethetamagn = T40(:,8);
    Ethetaphase = pi/180*T40(:,9);
    Ephimagn = T40(:,10);
    Ephiphase = pi/180*T40(:,11);

    %size of E
    Esqr = Ethetamagn.^2+Ephimagn.^2;

    Int = c*e0/2*Esqr; %intesnity in each (theta,phi)

    % sort values for ploting
    for n = 0:72
        for m = 2:19
            INT(n*19+m) = Int((n+1)*19-(m-2));
        end
    end
    
    r = R*tan(theta);
 
    for n = 1:19
        rs(n) = r(n);
    end
    
    RS = sort(rs);

    for n = 1:18
        RS2(n)=RS(n);
    end
    
    for n = 1:16
        RS3(n) = RS2(n);
    end
    for m = 1:73
        for n = 1:16
            I0(n,m) = INT((m-1)*19+n);
        end
    end
    
    %create mesh
    [rr,tt] = meshgrid(RS3,linspace(-pi,pi,73));
    xx = rr.*cos(tt);
    yy = rr.*sin(tt);
    %intensity in dB
    I0dB = 10*log10(I0/(max(max(I0))));

    %plots all beams in one plot
    figure(4);
    subplot(korningar,2,korn*2-1);
    surf(xx,yy,I0'); view(2); axis([-4 4 -4 4]);
    shading interp; axis image; axis off;
    subplot(korningar,2,korn*2);
    surf(xx,yy,I0dB'); view(2); axis([-4 4 -4 4]);
    shading interp; axis image; axis off; colorbar;

    %plots all beams in individual plots
    figure(5+korn);
    surf(xx,yy,I0'); view(2); axis([-4 4 -4 4]);
    shading interp; axis image;
    axis off;
    figure(15+korn);
    surf(xx,yy,I0dB'); view(2); axis([-4 4 -4 4]);
    shading interp; axis image; axis off; colorbar;
end

