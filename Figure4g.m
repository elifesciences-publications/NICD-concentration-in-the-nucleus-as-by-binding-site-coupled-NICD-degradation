%The degradation model without coreprresor - with binding unbinding of
%phosphorylated notch and only degradation of NICD_ub via Gamma_p and
%Gamma_up
clear all
%The parameters

NICD0=2000;
slope_wt=-0.0294;
Ne=11;
Gamma_p=1/8; %the degradation rate of phosphorylated NICD
Gamma_up=1/120; %the degradaion rate of unphosphorilated NICD
kpp= NICD0/(-1/slope_wt+Ne);
P_NICD=NICD0*Gamma_up; %the rate of incoming NICD into the nucleus
% kp0=0.415;     %the phosphorilation rate of bound NICD
kp0=kpp*Gamma_up/(1-Gamma_up/Gamma_p);
k_alpha=12; %the affinity of NICD to the binding sites
number_of_sites=0:50;
number_of_sites2=0:25;

%% Plot of NICD levels vs #SPS sites for WT and CDK+/- 
% Option I : # of SPS sites includes both endogenous and for the degrader
% figure
% k_p=kp0; %WT
% NICD_tot_wt=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% plot(number_of_sites, NICD_tot_wt, 'm', 'linewidth', 3);
% k_p=kp0/2; %CDK8 het
% NICD_tot_0=NICD_tot_wt(1);
% hold on
% NICD_tot_skd=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
% plot(number_of_sites, NICD_tot_skd, 'r', 'linewidth', 3);
% legend('WT','CDK8^{+/-}')
% legend('boxoff')
% axis([0 50 0 1.2*NICD0])
% set(gca,'fontsize',16)

% Option II : # of SPS sites includes only those of the degrader

numWT=Ne;
number_of_sites=numWT:numWT+50;
k_p=kp0; %WT
NICD_tot_wt=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);
NICD_tot_0=NICD_tot_wt(1);
k_p=kp0/2; %CDK8 het
NICD_tot_skd=degrader2_rebinding(P_NICD, Gamma_p, Gamma_up, k_p, k_alpha, number_of_sites);


figure
plot((number_of_sites-numWT), NICD_tot_wt, 'm', 'linewidth', 3);
hold on
plot((number_of_sites-numWT), NICD_tot_skd, 'r', 'linewidth', 3);
legend('WT','skd^{+/-}')
legend('boxoff')
axis([0 50 0 1.2*NICD0])
set(gca,'fontsize',16)

Sim.WT=NICD_tot_wt;
Sim.Ne=number_of_sites-numWT;
Sim.CDK=NICD_tot_skd;
save('D:\My Documents\Dropbox (TAU Storage)\Brian_Rafi_David\Model files\Degradation_only_CoA_with_rebinding\Sim.mat','Sim')

%% linear fit to data 

%fit results for y=a1+a2x
%YW (fit only to the first 3 points without the 36 sps)
%a1=1.05899 ± 0.09705319
%a2 = -0.02940908 ± 0.00296674

%Sdk
%a1=1.20694 ± 0.07668673
%a2=-0.0156678 ± 0.001887674

%the data
%YW
x_YW(1:10)=0;
x_YW(11:21)=12;
x_YW(22:34)=24;
x_YW(35:42)=36;
% rand_x_YW= - 1 + 0.5.*(rand(1, length(x_YW)) - rand(1, length(x_YW)));
% x_YW_rand=x_YW + rand_x_YW;
x_YW_rand=-1+x_YW;
y_YW=[1.09583406 1.225732159 0.918434366 1.038126397 0.980958086 0.808458838 0.779967492 1.44343489 0.807801359 0.901252355 0.879244773 0.854524044 1.13663605 1.046587277 0.677410107 0.652400578 0.782833343 0.884525115 0.317499109 0.899489915 1.124034163 0.345639198 0.19399373 0.302055212 0.388858194 0.343890164 0.283858473 0.193737221 0.330279113 0.495715317 0.259680328 0.210794856 0.563106295 0.490583683 0.195853604 0.312580853 0.374239507 0.227269407 0.346150229 0.121552368 0.164824787 0.284364866 ];

%Skd
x_Skd(1:27)=0 ;
x_Skd(28:41)=12 ;
x_Skd(42:60)=24 ;
x_Skd(61:72)=36 ;
% rand_x_Skd=  1 + 0.5.*(rand(1, length(x_Skd)) - rand(1, length(x_Skd)));
% x_Skd_rand=x_Skd + rand_x_Skd;
x_Skd_rand= 1 + x_Skd;
y_Skd=[0.811493413 1.018353451 1.317878366 1.494930337 0.945871887 1.104723397 1.066965067 1.235363455 1.161144727 1.43764649 0.959268912 1.402704078 1.158712265 1.056582984 1.161611014 1.073310263 0.821086218 1.608339094 1.2758261 1.792728154 1.319049747 0.985665092 0.5846247 1.162404855 0.977717987 1.179409818 1.074748724 0.800509431 1.211135433 1.17499132 0.886982855 1.006809894 0.997829273 1.047117705 0.979607146 1.020876295 1.127503487 1.101942299 0.842328245 1.241132995 1.254806336 0.565595065 0.651731086 0.997118067 1.138908543 0.82954887 1.126463152 1.30258213 1.36160105 1.461405629 1.273042506 1.161881645 1.066240997 0.902263877 0.597238798 0.696510427 0.63126025 0.565944165 0.74124442 0.743613208 0.781595216 0.408004069 0.913478742 0.769292738 0.461253347 0.620261944 0.691616512 0.639558617 0.404389409 0.370223573 0.512050692 0.424173047];

%the fits
x_fit=0:40;
y_fit_YW=1.05899 -0.02940908.*x_fit;
y_fit_Skd=1.20694 -0.0156678.*x_fit;

%the fits negative
x_fit_neg=-11:0;
y_fit_YW_neg=1.05899 -0.02940908.*x_fit_neg;
y_fit_Skd_neg=1.20694 -0.0156678.*x_fit_neg;


%average data
GFP_YW=([1 0.84 0.34 0.253])';
dGFP_YW=([0.0948280340332811 0.0562281329109105 0.12617033611488 0.0907790131760148])';

GFP_Skd=([1.16 1.05 0.94 0.58])';
dGFP_Skd=([0.0312166418074514 0.0542437700725964 0.0549447784308477 0.11205396414483])';


%plotting the data
figure
plot(x_YW_rand, y_YW, 'db','markersize',4)
hold on
plot(x_Skd_rand, y_Skd, 'or','markersize',4)
plot(x_fit, y_fit_YW, 'b','linewidth',3)
plot(x_fit, y_fit_Skd, 'r','linewidth',3)

norm=NICD_tot_0/1.05899;
plot((number_of_sites-numWT), NICD_tot_wt/norm, 'b-.', 'linewidth', 2);
plot((number_of_sites-numWT), NICD_tot_skd/norm, 'r-.', 'linewidth', 2);

plot(x_fit_neg, y_fit_Skd_neg, 'r--','linewidth',3)
plot(x_fit_neg, y_fit_YW_neg, 'b--','linewidth',3)
errorbar(([0 12 24 36] -1), GFP_YW, dGFP_YW, 'b.' ,'linewidth',2)
errorbar(([0 12 24 36]+1), GFP_Skd, dGFP_Skd, 'r.' ,'linewidth',2)
xticks([0 12 24 36])
ylim([0 2])
xlim([-11 40])
set(gca,'FontSize',16)
lgd=legend('wt','skd^{+/-}','fit wt','fit skd^{+/-}','model wt','model skd^{+/-}');
lgd.FontSize=12;
% legend('boxoff')

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

