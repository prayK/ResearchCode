clc;
close all;
clear all;
warning off;
addpath 'func\'

%% 
c           = 1;
%OAM起始点
OAM_begin   = 0;
Wavelength  = 1.5; 
%幅度
A           = 0.5; 
NUM         = 100;  
Times2      = 10;  
Ranges      = 2*Wavelength;  
k           = 2*pi/Wavelength;
w           = k*c;
kvt         = [1 0 0]*k;
zh          = linspace(0,Ranges,NUM);  
zL          = length(zh);
r           = 0.5*ones(1,NUM);
steps       = Wavelength/(2*2*pi*c); 
Tvt         = 0:steps:Times2;
Zr          = Ranges/Wavelength;
steg        = zeros(1,Zr); 
t           = 1;
n           = 1;
l           = OAM_begin;
while t <= Times2
    while l <= OAM_begin+3
        if l==0
            clf;
            subplot(2,2,n);
            hold on;
            axis([0 Ranges -2 2 -2 2]);
            view(51,30);
            grid on;
            
            z_value_tmp       = zh;
            phvt              = linspace(0,2*pi,NUM);
            phL               = length(phvt);
            R                 = ones(1,phL);
            [x_value ,y_value]= pol2cart(phvt,R);
            temp              = ones(1,zL);
            zh                = temp*w*Tvt(t)/k;
            clear temp;
            for p=1:Zr
                zvp = zh+(p-1)*Wavelength;
                if zvp >= (steg(1,p)+1)*Ranges
                   steg(1,p)=steg(1,p)+1;
                end
                zvp = zvp-2*steg(p)*Wavelength;
                m   = 1;
                while m < zL  
                    X=[zvp(m) ;zvp(m) ;zvp(m+1);zvp(m+1) ];
                    Y=[0 ;y_value(m) ;y_value(m+1);0 ];
                    Z=[0 ;x_value(m) ;x_value(m+1);0 ];
                    fill3(X,Y,Z,'y');
                    m=m+1;
                end
            end
            zh=z_value_tmp;
            else
            subplot(2,2,n);
            hold on;
            axis([0 Ranges -2 2 -2 2]);
            view(51,30);
            grid on
            l_steg=0;
            while l_steg<l
                temp              = ones(1,zL);
                Phi_begin         = l_steg*2*pi/l*temp;
                phvt              = Phi_begin+k*zh/l-w*Tvt(t)/l;
                phL               = length(phvt);
                R                 = ones(1,phL);
                [x_value,y_value] = pol2cart(phvt,R);
                clear temp;
                m                 = 1;
                while m<zL
                    X=[zh(m) ;zh(m) ;zh(m+1);zh(m+1) ];
                    Y=[0 ;y_value(m) ;y_value(m+1);0 ];
                    Z=[0 ;x_value(m) ;x_value(m+1);0 ];
                    fill3(X,Y,Z,'y')
                    m=m+1;
                end
                l_steg=l_steg+1;
            end
        end
        l=l+1;
        n=n+1;
    end 
    hold off;
    n = 1; 
    t = t+1;
    pause(0.01);
    l = OAM_begin;
end


%%     
Gmin  = [4 4 4 4];      
Gmax  = [14 20 30 35];   
Gstep = 5;               
steps  = 4;                   
figure;
for m=1:4
    if m == 1
       Data = load('d1.mat');    
    end
    if m == 2
       Data = load('d2.mat');  
    end    
    if m == 3
       Data = load('d3.mat');  
    end    
    if m == 4
       Data = load('d4.mat');  
    end
    for grad = Gmin:Gstep:Gmax   
        grad
        [Ptheta,Pphi] = func_phase_get(Data.Data,grad,steps);
        d             = length(Ptheta);
        for k=1:d
            if Ptheta(k,2)>=180
                if Ptheta(k,4)>0
                   Ptheta(k,4)=Ptheta(k,4)-180;
                else
                   Ptheta(k,4)=Ptheta(k,4)+180;
                end
                Ptheta(k,3) = Ptheta(k,3)*(-1);
            end
            if Pphi(k,2)>=90 & Pphi(k,2)<270
                if Pphi(k,4)>0
                   Pphi(k,3)=Pphi(k,3)*(-1);
                   Pphi(k,4)=Pphi(k,4)-180;
                else
                   Pphi(k,3)=Pphi(k,3)*(-1);
                   Pphi(k,4)=Pphi(k,4)+180;
                end
            end
        end
        R          = 1;             
        r          = R.*tan(pi/180*grad).*ones(length(Ptheta),1);  
        [X,Y]      = pol2cart(pi/180.*Ptheta(:,2),r);
        [X_th,Y_th]= pol2cart(pi/180.*Ptheta(:,4),1);
        [X_ph,Y_ph]= pol2cart(pi/180.*Pphi(:,4),1);
        subplot(2,2,m);
        quiver(X,Y,X_ph,Y_ph,0.5);
        hold on;
        xlabel('X');
        ylabel('Y');
    end
end
 

