

NICD0=5000;
slope_wt=-0.0294;
err_slope=0.0030;
Ne=11;
Gamma_p=1/8; %the degradation rate of phosphorylated NICD
Gamma_up=1/120; %the degradaion rate of unphosphorilated NICD
P_NICD=NICD0*Gamma_up; %the rate of incoming NICD into the nucleus
% kp0=0.415;     %the phosphorilation rate of bound NICD

k_alpha=12; %the affinity of NICD to the binding sites


NICD0=100:100:10000;
kpp= NICD0/(-1/slope_wt+Ne);

Gamma_up=1/120;

kp0=kpp*Gamma_up/(1-Gamma_up/Gamma_p);

kpp_errp= NICD0/(-1/(slope_wt+err_slope)+Ne);
kp0_errp=kpp_errp*Gamma_up/(1-Gamma_up/Gamma_p);

kpp_errm= NICD0/(-1/(slope_wt-err_slope)+Ne);
kp0_errm=kpp_errm*Gamma_up/(1-Gamma_up/Gamma_p);


figure

plot(NICD0,kp0,'k-','linewidth',2)
hold on
plot(NICD0,kp0_errp,'k--','linewidth',1,'color',[0.5 0.5 0.5])
plot(NICD0,kp0_errm,'k--','linewidth',1,'color',[0.5 0.5 0.5])


set(gca,'fontsize',16)

figure

loglog(NICD0,kp0,'k-','linewidth',2)
hold on
loglog(NICD0,kp0_errp,'--','linewidth',1,'color',[0.5 0.5 0.5])
loglog(NICD0,kp0_errm,'--','linewidth',1,'color',[0.5 0.5 0.5])


figure
Gamma_up_vals=[1/30 1/60 1/120 1/240 1/480];

for i=1:5
    
    Gamma_up=Gamma_up_vals(i);
    
    kp0=kpp*Gamma_up/(1-Gamma_up/Gamma_p);
    loglog(NICD0,kp0,'-','linewidth',2)
    hold on
    % loglog(NICD0,kp1_errp,'--','linewidth',1,'color',[0.5 0.5 0.5])
    % loglog(NICD0,kp1_errm,'--','linewidth',1,'color',[0.5 0.5 0.5])
    axis([100 10000 10^-2 10^1])
end
legend({'\Gamma_{up} = 1/30 min^{-1}','\Gamma_{up} = 1/60 min^{-1}','\Gamma_{up} = 1/120 min^{-1}','\Gamma_{up} = 1/240 min^{-1}','\Gamma_{up} = 1/480 min^{-1}'},'Location','northwest','FontSize',11);
set(gca,'fontsize',16)



figure
Gamma_up_vals=logspace(-1.5,-3,100);
NICD0=logspace(2,4,100);
kpp = NICD0/(-1/slope_wt+Ne);
kp0 = zeros(length(NICD0),length(Gamma_up_vals)); 

for i=1:length(Gamma_up_vals)
    
    Gamma_up=Gamma_up_vals(i);
    
    kp0(:,i)=kpp*Gamma_up/(1-Gamma_up/Gamma_p);
end

imagesc(NICD0,Gamma_up_vals,log10(kp0))
colorbar('Ticks',[-2,-1,0,1],...
         'TickLabels',{'10^{-2}','10^{-1}','10^{0}','10^{1}'})
set(gca,'fontsize',16)

figure

surf(NICD0,Gamma_up_vals,log10(kp0),'EdgeColor','none')
colorbar('Ticks',[-2,-1,0,1],...
         'TickLabels',{'10^{-2}','10^{-1}','10^{0}','10^{1}'})
set(gca,'fontsize',16)
