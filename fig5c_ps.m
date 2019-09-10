function [t1 x1_tot t2 x2_tot t3 x3_tot]=fig5c_ps(params);

% close all

if nargin<1
    
    NICD0=2000;
    slope_wt=-0.0294;
    slope_skd=-0.0156678;
    a_skd=1.20694;
    a_wt=1.05899;
    %
    % slope_wt=-0.05;
    % slope_skd=slope_wt/2;
    
    Ne=(a_skd-a_wt)/(slope_skd-slope_wt);
    kpp= NICD0/(-1/slope_wt+Ne);
    Gamma_up=1/120;
    Gamma_p=1/8;
    k_p=kpp*Gamma_up/(1-Gamma_up/Gamma_p);
    P_NICD=NICD0*Gamma_up;
    
    params.P_NICD=P_NICD;
    params.Gamma_up=Gamma_up;
    params.Gamma_p=Gamma_p;
    params.k_p=k_p;
    params.k_a=10;
    params.Ne=Ne;
    params.tend=300;
end

tspan=[0 params.tend];
IC=[0 0];

%% 24x case
params.num_sites=params.Ne+24;
%steady state 24x
[t1 x1]=ode23(@degrader2d,tspan,IC,[],params);
x1_tot=x1(:,1)+x1(:,2);

%% case wildtype

% steady state 2x
params.num_sites=params.Ne;

[t2 x2]=ode23(@degrader2d,tspan,IC,[],params);
x2_tot=x2(:,1)+x2(:,2);

%% case Notch Het

% steady state 2x
params.num_sites=params.Ne;
params.P_NICD=params.P_NICD/2;
[t3 x3]=ode23(@degrader2d,tspan,IC,[],params);
x3_tot=x3(:,1)+x3(:,2);



