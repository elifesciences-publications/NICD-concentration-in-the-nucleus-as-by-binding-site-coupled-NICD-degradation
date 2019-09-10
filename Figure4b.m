%The degradation model without coreprresor - with binding unbinding of
%phosphorylated notch and only degradation of NICD_ub via Gamma_p and
%Gamma_up
clear all
%The parameters

NICD0=2000;
slope_wt=-0.0294;
Ne=11;
kpp= NICD0/(-1/slope_wt+Ne);
Gamma_p=1/8; %the degradation rate of phosphorylated NICD
Gamma_up=1/120; %the degradaion rate of unphosphorilated NICD
P_NICD=NICD0*Gamma_up; %the rate of incoming NICD into the nucleus
% kp0=0.415;     %the phosphorilation rate of bound NICD
% kp0=kpp*Gamma_up*(1-Gamma_up/Gamma_p);
kp0=1; 
k_alpha=12; %the affinity of NICD to the binding sites
number_of_sites=0:50;
number_of_sites2=0:25;

number_of_SPSsites=number_of_sites/2;

%% Plot of NICD levels vs #SPS sites for WT and CDK+/- 
% Option I : # of SPS sites includes both endogenous and for the degrader
figure
k_p=kp0; %WT
NICD_tot=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
plot(number_of_SPSsites, NICD_tot, 'm', 'linewidth', 3);
k_p=kp0/2; %CDK8 het
hold on
NICD_tot=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
plot(number_of_SPSsites, NICD_tot, 'r', 'linewidth', 3);
k_p=kp0/4; %CDK8 het
hold on
NICD_tot=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
plot(number_of_SPSsites, NICD_tot, 'b', 'linewidth', 3);
legend('k_{p} = 1 min^{-1}','k_{p} = 0.5 min^{-1}','k_{p} = 0.25 min^{-1}')
legend('boxoff')
axis([0 25 0 1.2*NICD0])
set(gca,'fontsize',16)

% Option II : # of SPS sites includes only those of the degrader
% figure
% numWT=Ne;
% number_of_sites=numWT:numWT+50;
% k_p=kp0; %WT
% NICD_tot=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% Sim.WT=NICD_tot;
% Sim.Ne=number_of_sites-numWT;
% plot((number_of_sites-numWT), NICD_tot, 'm', 'linewidth', 3);
% k_p=kp0/2; %CDK8 het
% hold on
% NICD_tot=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% Sim.CDK=NICD_tot;
% plot((number_of_sites-numWT), NICD_tot, 'r', 'linewidth', 3);
% legend('WT','skd^{+/-}')
% legend('boxoff')
% axis([0 50 0 1.2*NICD0])
% set(gca,'fontsize',16)

numWT=Ne;
Sim.WT=NICD_tot;
Sim.Ne=number_of_sites-numWT;
Sim.CDK=NICD_tot;
save('D:\My Documents\Dropbox (TAU Storage)\Brian_Rafi_David\Model files\Degradation_only_CoA_with_rebinding\Sim.mat','Sim')

%% plot phosphorylated vs unphosphorylated
% [NICD_tot1, NICD_p1, NICD_up1]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% [NICD_tot1lin, NICD_p1lin, NICD_up1lin]=degrader2_linear(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites2);
% 
% k_p=18;
% 
% 
% [NICD_tot2, NICD_p2, NICD_up2]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% [NICD_tot2lin, NICD_p2lin, NICD_up2lin]=degrader2_linear(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites2);
% 
% figure
% plot(number_of_sites, NICD_tot1, 'b', 'linewidth', 3);
% hold on
% plot(number_of_sites2, NICD_tot1lin, 'b--', 'linewidth', 3);
% % plot(number_of_sites, NICD_tot2, 'r', 'linewidth', 3);
% % plot(number_of_sites2, NICD_tot2lin, 'r--', 'linewidth', 3);
% % 
% plot(number_of_sites, NICD_p1, 'g.', 'linewidth', 1.5);
% plot(number_of_sites, NICD_up1, 'r.', 'linewidth', 1.5);
% 


% k_alpha=5;
% [NICD_tot2, NICD_p2, NICD_up2]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% k_alpha=10;
% k_p=2.5;
% [NICD_tot3, NICD_p3, NICD_up3]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% k_p=5.0;
% P_NICD=250;
% [NICD_tot4, NICD_p4, NICD_up4]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% k_p=2.5;
% P_NICD=250;
% [NICD_tot2, NICD_p2, NICD_up2]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% plot(number_of_sites, NICD_p1./NICD_up1, 'b', 'linewidth', 3);
% hold on
% plot(number_of_sites, NICD_p3./NICD_up3,'k', 'linewidth', 3 );
% plot(number_of_sites, NICD_p4./NICD_up4,'r', 'linewidth', 3 );
% plot(number_of_sites, NICD_p2./NICD_up2,'--g', 'linewidth', 3 );


% plot(number_of_sites, NICD_up, 'b', 'linewidth', 3);
% plot(number_of_sites, NICD_p, 'r', 'linewidth', 3);
% P_NICD=250;
% NICD_tot1=degrader(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% [NICD_tot2, NICD_p, NICD_up]=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% plot(number_of_sites, NICD_tot1, 'k', 'linewidth', 3);
% plot(number_of_sites, NICD_tot2, '--k', 'linewidth', 3);

