% Plots the phase planes for four consecutive l values, in subfigures.
% LG beam.

clear;

%----------------constants-------------------

c=1;

%-------------things to change---------------

l_start=0;                   %start OAM
lam=1.5;                     %Wavelength
A=0.5;                       %Amplitude
nog=100;
t_stop=10;                   %Amount of timesteps
z_range=2*lam;               %Length of z values

%--------------------------------------------

k=2*pi/lam;
w=k*c;
k_vector=[1 0 0]*k;

%---------vector x--------

z_value=linspace(0,z_range,nog);     %z on x-axis in plot
L_z=length(z_value);

r=0.5*ones(1,nog);

%----- --------- -------

t_step=lam/(2*2*pi*c);               %Not to big timesteps
t_vector=0:t_step:t_stop;

%------------------------

z_antal=z_range/lam;
steg=zeros(1,z_antal);               

t=1;
n=1; 
l=l_start; 
while t<=t_stop
      while l<=l_start+3
           
      if l==0 %Plane wave.
         clf;
         subplot(2,2,n);
         hold on;
   
         axis([0 z_range -2 2 -2 2]);
         view(51,30);
    
         axis fill;
         axis square;
         grid on;
         %light
    
         z_value_tmp=z_value;
         PHI_vector=linspace(0,2*pi,nog);
         L_PHI=length(PHI_vector);
         R=ones(1,L_PHI);
         [x_value y_value]=pol2cart(PHI_vector,R);
         temp=ones(1,L_z);
         z_value=temp*w*t_vector(t)/k;
         clear temp;

         for p=1:z_antal 
             z_value_p=z_value+(p-1)*lam;
             if z_value_p>=(steg(1,p)+1)*z_range
                 steg(1,p)=steg(1,p)+1;
             end
             
             z_value_p=z_value_p-2*steg(p)*lam;
             m=1;
             while m<L_z             %Plottar Phasplanes
                 X=[z_value_p(m) ;z_value_p(m) ;z_value_p(m+1) ;z_value_p(m+1)];
                 Y=[0            ;;y_value(m)  ;y_value(m+1)   ;0             ];
                 Z=[0            ;;x_value(m)  ;x_value(m+1)   ;0             ];
                 fill3(X,Y,Z,'y');
                 m=m+1;
                 alpha(1);
             end
             
         end
         %for loop
         
         z_value=z_value_tmp;
      
      else
          subplot(2,2,n);
          hold on;
          axis([0 z_range -2 2 -2 2]);
          view(51,30);
          
          axis fill;
          axis square;
          grid on;
          
          xlabel('x');
          ylabel('y');
          zlabel('z');
          
          light;
          material shiny;
          l_steg=0;
          while l_steg<l
              
              temp=ones(1,L_z);
              PHI_start=l_steg*2*pi/l*temp;
              PHI_vector=PHI_start+k*z_value/l-w*t_vector(t)/l;
              L_PHI=length(PHI_vector);
              R=ones(1,L_PHI);
              
              [x_value y_value]=pol2cart(PHI_vector,R);
              clear temp;
              
              m=1;
              while m<L_z                 %%Plottar Phasplanes
                  X=[z_value(m) ;z_value(m) ;z_value(m+1) ;z_value(m+1)];
                  Y=[0          ;y_value(m) ;y_value(m+1) ;0           ];
                  Z=[0          ;x_value(m) ;x_value(m+1) ;0           ];
                  fill3(X,Y,Z,'y')
                  m=m+1;
              
              end
              l_steg=l_steg+1;
          end   %while l_steg...
      end   %else satsen
      l=l+1;
      n=n+1;
      
      end %while l.

      %Film(:,t)=getframe(gcf);
      hold off;
      n=1;            %?terst?llning av r?knare
      t=t+1;          %Tidsstegning
      pause(0.1);
      l=l_start;
end
%movie2avi(Film,'phaseplanes_4.avi');
%movieview('phaseplanes_4.avi')