%% part two: fractional equilibration with interannual noise
clear all;
color1 = [1 0.5 0.1];
color2 = [0.2 0.6 1];
color3 = 0.2*[1 1 1];
load('Cascades_tau12_schematic','tplot','noisy_retreat_SNR1percentury','Leq','bp','bramp')
noisy12 = noisy_retreat_SNR1percentury;
Leq12 = Leq;
lw = 1;
tau1 = 12;
dt = 1;
eps = 1/sqrt(3);
kap = 1 - dt*(eps.*tau1).^(-1);
sigb = 1; % m/yr
bdot = 0.01; % m/yr /yr
psi = tau1.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
t = 1:170;
sigfrac = t.^(-1) * (psi/tau1) * (sigb/bdot);
frac3_12 = 1 - (3*eps*tau1./t).*(1-exp(-t/(eps*tau1))) + exp(-t/(eps*tau1)).*(t/(2*eps*tau1) + 2);

tplot2 = 1880+t;

subplot(3,3,4); hold on
patch([tplot2,flip(tplot2)],[frac3_12+sigfrac,flip(frac3_12-sigfrac)],color1,'edgecolor','none'); alpha(0.2);
plot(tplot2,frac3_12,':','color',color1,'linewidth',lw)
plot(tplot,(noisy12)./(Leq12-Leq12(1)),'color',color1,'linewidth',lw)

%for another glacier
load('Cascades_tau48_schematic','noisy_retreat_SNR1percentury','Leq')
noisy48 = noisy_retreat_SNR1percentury;
Leq48 = Leq;
tau2 = 48;
kap = 1 - dt*(eps.*tau2).^(-1);
psi = tau2.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
t = 1:170;
sigfrac = t.^(-1) * (psi/tau2) * (sigb/bdot);
frac3_48 = 1 - (3*eps*tau2./t).*(1-exp(-t/(eps*tau2))) + exp(-t/(eps*tau2)).*(t/(2*eps*tau2) + 2);
patch([tplot2,flip(tplot2)],[frac3_48+sigfrac,flip(frac3_48-sigfrac)],color2,'edgecolor','none'); alpha(0.2);
plot(tplot2,frac3_48,':','color',color2,'linewidth',lw)
plot(tplot,(noisy48)./(Leq48-Leq48(1)),'color',color2,'linewidth',lw)

ylabel('Fractional equilibration')
xlim([1880 2020])
xlabel('Year CE')
ylim([0 1])
set(gca,'fontsize',9,'xtick',[1880:40:2020])
%text(1890,0.9,'SNR_b = 1 century^{-1}','FontSize',9)
%title('a)')

axes('position',[0.19 0.65 0.1 0.08]); hold on; box on
%subplot(3,3,1); hold on
patch([tplot,flip(tplot)],[bramp+sigb;flip(bramp-sigb)],[1 0.6 0.6],'edgecolor','none'); alpha(0.5);
plot(tplot,bp,'color',[1 0.5 0.5]);
plot(tplot,bramp,':r','linewidth',lw);
title('Trend + noise (SNR = 1 century^{-1})')
ylim([-3 2]);
set(gca,'xtick',[1880:60:2020],'fontsize',9,'ytick',[])
xlim([1880 2020])
ylabel('b''')
%%
load('Cascades_tau12_hiatus');
tplot = t - yrs_pre + date_start; % adjust to year trend starts

subplot(3,3,5); hold on; 
plot(tplot,((L-L(1))./(Leq-Leq(1))),'color',color1,'linewidth',lw)    
plot(tplot2,frac3_12,':','color',color1,'linewidth',lw)

load('Cascades_tau48_hiatus');

plot(tplot,((L-L(1))./(Leq-Leq(1))),'color',color2,'linewidth',lw)  
plot(tplot2,frac3_48,':','color',color2,'linewidth',lw)

ylabel('Fractional equilibration')
xlim(plotlims)
xlabel('Year CE')
ylim([0 1])
set(gca,'fontsize',9,'xtick',[1880:40:2020])

%subplot(3,3,2)
axes('position',[0.47 0.65 0.1 0.08]); 
plot(tplot,-mu*(Tramp-Tramp(1)),'r','linewidth',lw);
title('30-yr hiatus')
xlim([1880 2020])
set(gca,'fontsize',9,'xtick',[1880:60:2000],'ytick',[]); grid on
ylabel('b''')

% now compare with uncertainty in tau...

tau = tau1*[0.5 1 1.5]'; tau = repmat(tau,1,170);
kap = 1 - dt*(eps.*tau).^(-1);
psi = tau.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
t = repmat([1:170],3,1);
frac3 = 1 - (3*eps*tau./t).*(1-exp(-t./(eps*tau))) + exp(-t./(eps*tau)).*(t./(2*eps*tau) + 2);

subplot 336; hold on
%subplot('position',[0.35 0.1 0.21 0.21]); hold on
f12=patch([tplot2,flip(tplot2)],[frac3(1,:),flip(frac3(3,:))],color1,'edgecolor','none'); alpha(0.2);
plot(tplot2,frac3(2,:),'color',color1);

tau = tau2*[0.5 1 1.5]'; tau = repmat(tau,1,170);
kap = 1 - dt*(eps.*tau).^(-1);
psi = tau.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
t = repmat([1:170],3,1);
frac3 = 1 - (3*eps*tau./t).*(1-exp(-t./(eps*tau))) + exp(-t./(eps*tau)).*(t./(2*eps*tau) + 2);

f48=patch([tplot2,flip(tplot2)],[frac3(1,:),flip(frac3(3,:))],color2,'edgecolor','none'); alpha(0.2);
plot(tplot2,frac3(2,:),'color',color2);
text(1890,0.9,'\sigma_{\tau} = 50%','FontSize',10,'fontweight','bold')

ylabel('Fractional equilibration')
xlim([1880 2020])
xlabel('Year CE')
ylim([0 1])
set(gca,'fontsize',9,'xtick',[1880:40:2020])
legend([f12 f48],'\tau = 12 yr','\tau = 48 yr','location','northoutside')
set(legend,'fontsize',12)