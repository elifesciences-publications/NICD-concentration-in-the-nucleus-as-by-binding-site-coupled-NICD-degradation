%this function calculated the total NICD with no CoR in for the degradation
%coupled transcripitin model. There is WITH effect of rebinding of phosphorylated NICD
function [NICD_tot, NICD_p, NICD_up]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites)
NICD_tot=zeros(1, length(number_of_sites));
NICD_p=zeros(1, length(number_of_sites));
NICD_up=zeros(1, length(number_of_sites));

site1=number_of_sites(1);

% for N= 0:(length(number_of_sites)-1)     %the number of binding sites
for N = number_of_sites     %the number of binding sites
A= -Gamma_p.*Gamma_up; %the term of NICD_tot^3 in the cubic equation
B= P_NICD.*Gamma_p-N.*k_p.*Gamma_p - 2.*k_alpha.*Gamma_p.*Gamma_up + 2.*N.*Gamma_p.*Gamma_up; %the term of NICD tot^2
C= N.*P_NICD.*k_p - N.*Gamma_p.*k_alpha.*k_p + N.*N.*k_p.*Gamma_p - Gamma_p.*Gamma_up.*k_alpha.*k_alpha +2.*N.*Gamma_p.*Gamma_up.*k_alpha - N.*N.*Gamma_p.*Gamma_up + 2.*P_NICD.*Gamma_p.*k_alpha - N.*P_NICD.*Gamma_p;
D= N.*P_NICD.*k_alpha.*k_p + P_NICD.*Gamma_p.*k_alpha.*k_alpha - N.*P_NICD.*Gamma_p.*k_alpha; %the term of NICD_tot^0

p=[A B C D]; %defining the polynomial

r=roots(p); %findnig the roots of the polynomial
NICD_tot(N-site1+1)=max(r); %taking only the max of the 3 solutions (the rest are negative or don't fit)
NICD_up(N-site1+1)= P_NICD.*(k_alpha + NICD_tot(N-site1+1))./(N.*k_p +Gamma_up.*(k_alpha + NICD_tot(N-site1+1)- N-site1));
NICD_p(N-site1+1)=NICD_tot(N-site1+1)-NICD_up(N-site1+1);
end
end