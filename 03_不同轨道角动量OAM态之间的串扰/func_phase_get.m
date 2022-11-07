function [Ptheta,Pphi] = func_phase_get(L,grad,step)
p = 1;
o = 1;
d = length(L);    

for l=1:d              
    if L(l,1) == -1*double(grad)
        if p==1
           Phase_theta1 = [L(l,1), L(l,2), L(l,8), L(l,9)];
           p            = 2;
        else
           Phase_theta1 = [Phase_theta1;L(l,1), L(l,2), L(l,8), L(l,9)];
        end
        if o==1
           Phase_phi1   = [L(l,1), L(l,2), L(l,10), L(l,11)];
           o            = 2;
        else
           Phase_phi1   = [Phase_phi1; L(l,1), L(l,2), L(l,10), L(l,11)];
        end
    end
end
Phase_theta(:,:) = Phase_theta1; 
Phase_phi(:,:)   = Phase_phi1; 
l_ph_theta       = length(Phase_theta(:,1));
for ten=1:step:l_ph_theta
    if ten==1;
        Ptheta(1,:)=Phase_theta(ten,:);
        Pphi(1,:)=Phase_phi(ten,:);
    else
        Ptheta((ten-1)/step+1,:)=Phase_theta(ten,:);
        Pphi((ten-1)/step+1,:)=Phase_phi(ten,:);
    end
end