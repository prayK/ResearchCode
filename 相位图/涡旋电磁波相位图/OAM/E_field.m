%
% This program plots what the electric field looks like from 4nec2 data.


% Source files from 4nec, modified. "LINEAR", "LEFT", "RIGHT" replaced
  %by appropriate amount of space. Radiation patterns values are kept.
  %Columns as follows:
   %( 1 - Theta               2 - Phi                        ) Angles                        
   %( 3 - Vert. DB            4 - Hor. DB       5 - Tot. DB  ) Power gains 
   %( 6 - Axial Ratio         7 - Tilt deg                   ) Polarization     
   %( 8 - E_theta magnitude.  9 - Phase_theta                ) E(Theta)
   %(10 - E_phi magnitude.   11 - Phase_phi                  ) E(Phi)
       
clear;
close all;
       
names_of_files=['L0.dat';     %Names of 4nec2 files
                'L1.dat';
                'L2.dat';
                'L4.dat'];
       
grad_min=[4 4 4 4];        % Min meassurement theta angle.
grad_max=[14 20 30 35];    % Max Meassurement theta angle (ex mainlobe)
grad_steg=1;               % Resolution in 4nec2 file.
step=2;                    % Every step'th value to be plotted

for m=1:length(names_of_files(:,1))
    Data=load(names_of_files(m,:));         %Data from file
    for grad=grad_min:grad_steg:grad_max    %Degrees for data.
        [Phase_theta_10 Phase_phi_10]=nec(Data,grad,step);
                     % nec_org can also be used here, making it general.
                     % nec obtains the data from the 4nec2 output file.
    % E_phase:
        E_theta=pi/180.*Phase_theta_10(:,3).*(cos(pi/180.*Phase_theta_10(:,4)));
        E_phi=pi/180.*Phase_phi_10(:,3).*(cos(pi/180.*Phase_phi_10(:,4)));
    % Distance between each ring:
        R=1;
        z=0*ones(length(E_theta),1);
        rho=R.*tan(pi/180*grad).*ones(length(E_theta),1);
    % Placement:
        X=rho.*cos(pi/180.*Phase_theta_10(:,2));
        Y=rho.*sin(pi/180.*Phase_theta_10(:,2));
        Z=z;
        
        grad_theta=-pi/180.*Phase_theta_10(:,1);
        grad_phi=pi/180.*Phase_theta_10(:,2);
    % Change of coordinate system:
        EX=E_theta.*cos(grad_phi).*cos(grad_theta)-E_phi.*sin(grad_phi);
        EY=E_theta.*cos(grad_theta).*sin(grad_phi)+E_phi.*cos(grad_phi);
        EZ=-E_theta.*sin(grad_theta);
    % PLOT:
        figure(m)
        quiver3 (X,Y,Z,EX,EY,EZ,0.3);        % Vector arrows
        hold on;
        ax_values=(R*tan(pi/180*grad)+0.1);  % Axes values
        axis ([-ax_values ax_values -ax_values ax_values -ax_values ax_values]);
        axis square;
        xlabel('X','fontsize',24);
        ylabel('Y','fontsize',24);
        zlabel('Z','fontsize',24);
    end %GRAD 
end

