function dx=degrader2d(t,x,params)

dx=zeros(2,1);

% x(1) is NICD_up
% x(2) is NICD_p

P_NICD=params.P_NICD;
Gamma_up=params.Gamma_up;
Gamma_p=params.Gamma_p;
k_p=params.k_p;
k_a=params.k_a;
num_sites=params.num_sites;

a_up = x(1)./k_a;
a_p = x(2)./k_a;

NICD_pb=num_sites.*(a_p./(1+a_p+a_up));
NICD_upb=num_sites.*(a_up./(1+a_p+a_up));


dx(1) = P_NICD - Gamma_up.*x(1) - NICD_upb.*(k_p-Gamma_up);
dx(2) = NICD_upb.*k_p - Gamma_p.*x(2) + NICD_pb.*Gamma_p;

end