%% setting parameters

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

%% running simulations

[t1 x1_tot t2 x2_tot t3 x3_tot]=fig5c_ps(params);


%% ploting curves

figure;
%plot all variables
% t2=t2*Gamma_up;
plot(t2,x2_tot,'linewidth',3);
hold on
% plot(t2p,x2p_tot,'linewidth',3);
% t1=t1*Gamma_up;
plot(t1,x1_tot,'-','linewidth',3);
% plot(t1p,x1p_tot,'-.','linewidth',3);

% t3=t3*Gamma_up;
plot(t3,x3_tot,'-','linewidth',3);
% plot(t3p,x3p_tot,'--','linewidth',3);

xlabel('Time');
ylabel('NICD Concentration [#NICD/nucleus]')
% legend('WT steady state','WT pulse','WT+24xSPS steady state','WT+24xSPS pulse', 'N^{+/-} steady state','N^{+/-} pulse')
legend('WT','2 x G6S-LacZ','N^{+/-}')
set(gca,'fontsize',16)

