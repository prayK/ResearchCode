%
% This program will plot what phase arrows from 4nec2 data.
%
%Source files from 4nec, modified. "LINEAR", "LEFT", "RIGHT" replaced
  

  %by appropriate amount of space. Radiation patterns values are kept.
  %Columns as follows:
        %( 1 - Theta               2 - Phi                        ) Angles                        
        %( 3 - Vert. DB            4 - Hor. DB       5 - Tot. DB  ) Power gains
        %( 6 - Axial Ratio         7 - Tilt deg                   ) Polarization     
        %( 8 - E_theta magnitude.  9 - Phase_theta                ) E(Theta)
        %(10 - E_phi magnitude.   11 - Phase_phi                  ) E(Phi)

clear;

close all;

names_of_files=['L0.dat';      %Names of 4nec2 files
                'L1.dat';
                'L2.dat';
                'L4.dat'];
        
grad_min=[4 4 4 4];        % Min meassurement theta angle.
grad_max=[14 20 30 35];    % Max Meassurement theta angle (ex mainlobe)
grad_steg=5;               % Resolution in 4nec2 file.
step=4;                    % Every step'th value.

for m=1:length(names_of_files(:,1))
    Data=load(names_of_files(m,:));             %Data from file
    for grad=grad_min:grad_steg:grad_max        %Degrees for data.
        [Phase_theta_10 Phase_phi_10]=nec(Data,grad,step);
                        % nec obtains the data from the 4nec2 output file.
        % Ecos(theta)=-Ecos(theta+-180) same for phi. (For looks only)
        % No real change in the values of the E-field just easier to visualize the phase.
        d=length(Phase_theta_10);
        for k=1:d
            if Phase_theta_10(k,2)>=180
                if Phase_theta_10(k,4)>0
                    Phase_theta_10(k,4)=Phase_theta_10(k,4)-180;
                else
                    Phase_theta_10(k,4)=Phase_theta_10(k,4)+180;
                end
                Phase_theta_10(k,3)=Phase_theta_10(k,3)*(-1);
            end
            if Phase_phi_10(k,2)>=90 & Phase_phi_10(k,2)<270
                if Phase_phi_10(k,4)>0
                    Phase_phi_10(k,3)=Phase_phi_10(k,3)*(-1);
                    Phase_phi_10(k,4)=Phase_phi_10(k,4)-180;
                else
                    Phase_phi_10(k,3)=Phase_phi_10(k,3)*(-1);
                    Phase_phi_10(k,4)=Phase_phi_10(k,4)+180;
                end
            end
        end


% Distance between each ring:
     R=1;                %Virtual distance just to get a feel for it.
     r=R.*tan(pi/180*grad).*ones(length(Phase_theta_10),1); %Placements
% Placement:
     [X Y]=pol2cart(pi/180.*Phase_theta_10(:,2),r);
         %Since phi=theta, in our case we can use theta or phi
         %for the phase.CAREFUL!!!!!
% Change of coordinates:
     [X_th Y_th]=pol2cart(pi/180.*Phase_theta_10(:,4),1);
     [X_ph Y_ph]=pol2cart(pi/180.*Phase_phi_10(:,4),1);
% PLOT:
     figure(m)
        quiver (X,Y,X_ph,Y_ph,0.5);
        hold on;
        ax_values=(R*tan(pi/180*grad)+0.1);
        axis ([-ax_values ax_values -ax_values ax_values]);
        axis square;
        xlabel('X','fontsize',24);
        ylabel('Y','fontsize',24);
    end %GRAD
end %name_of_files