%
% Numerical Electromagnetic Code
%
% Returns the information, Theta, Phi, E_magn and phase) from 4nec2 data.
%
% The Data input file is obtained from the 4nec2 outputfile with modifications:
%   Only the values from "Radiation patterns" are kept.
%   "LINEAR", "RIGHT" and "LEFT" are replaced by appropriate amount of spaces.

function [Phase_theta_10 Phase_phi_10]=nec(L,grad,step)


% L               % Data file from 4nec2
% grad_min        % Min meassurement theta angle.
% grad_max        % Max meassurement theta angle (16 24 32 45 mainlobe)
% gradsteg        % Resolution in 4nec2.
% step            % Every step'th value.


p=1;o=1;
d=length(L);      
for l=1:d               % All values in L to be analyzed 
    if L(l,1)==-grad    % If the theta angle equals grad
        if p==1
            Phase_theta1=[L(l,1), L(l,2), L(l,8), L(l,9)];
            p=2;
        else
            Phase_theta1=[Phase_theta1;L(l,1), L(l,2), L(l,8), L(l,9)];
        end
        if o==1
            Phase_phi1=[L(l,1), L(l,2), L(l,10), L(l,11)];
            o=2;
        else
            Phase_phi1=[Phase_phi1; L(l,1), L(l,2), L(l,10), L(l,11)];
        end
    end
end

Phase_theta(:,:)=[Phase_theta1];    %[Theta         Phi
                                    %E_magn(theta)  E_phase(theta)]

Phase_phi(:,:)=[Phase_phi1];        %[Theta         Phi
                                    %E_magn(theta)  E_phase(theta)]


%___________________________ Every "step" value ___________________________
l_ph_theta=length(Phase_theta(:,1));
for ten=1:step:l_ph_theta
    if ten==1;
        Phase_theta_10(1,:)=Phase_theta(ten,:);
        Phase_phi_10(1,:)=Phase_phi(ten,:);
    else
        Phase_theta_10((ten-1)/step+1,:)=Phase_theta(ten,:);
        Phase_phi_10((ten-1)/step+1,:)=Phase_phi(ten,:);
    end
end