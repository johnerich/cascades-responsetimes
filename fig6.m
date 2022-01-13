%%% plots specifc and total melt flux changes for idealized glaciers under
%%% warming scenarios. (Fig. 6 of Christian et al., 2021)
% Output comes from 3-stage model (Roe and Baker 2014) with additional melt
% calculation bolted on, which simply integrates melt over the evolving
% glacier surface as the climate warms and the glacier retreats. 
% JEC 2021

clearvars;
figure(1); clf
color1 = [1 0.5 0.1];
color2 = [0.2 0.6 1];
date_start = 1880; % start of trend
date_end = 2020; % end of trend
plotlims = [date_start-20 date_end];
lw = 1;

%%% load output from glacier w/ 12 yr response time, linear trend
load('Cascades_tau12');
tplot = t - yrs_pre + date_start; % adjust to year trend starts
Qs = Qs*w; % multiply by nonminal width to get volume
Qs_eq = Qs_eq*w;

% plot total melt flux
axes('position',[0.33 0.45 0.28 0.25]); hold on; grid on
Q12=plot(tplot,Qs,'color',color1,'linewidth',lw-0.5);    
plot(tplot,Qs_eq,'-.','color',color1);

%%% load output from glacier w/ 48 yr response time
load('Cascades_tau48');
Qs = Qs*w; % multiply by nonminal width to get volume
Qs_eq = Qs_eq*w;

Q48=plot(tplot,Qs,'color',color2,'linewidth',lw)    
plot(tplot,Qs_eq,'-.','color',color2);
eq=plot(0,0,'-.','color','k');
title('   Total melt (M_g)','FontWeight','normal')
ylabel('m^3 / yr')
xlim(plotlims)
ylim([1e7 2.2e7])
xlabel('Year CE')
leg=legend([Q12 Q48 eq],'\tau = 12 yr','\tau = 48 yr','eq. response','location','best');
leg.FontSize = 6;
set(gca,'xtick',[1880:40:2000])
box on

%%% specific melt anomaly (e.g., increase in melt rate per area)
axes('position',[0.33 0.81 0.28 0.08]); 
plot(tplot,mu*(Tramp-Tramp(1)),'r','linewidth',lw);
title('specific melt anom.','FontWeight','normal')
ylabel('m / yr')
xlim([1860 2020])
ylim([0 1])
set(gca,'fontsize',8,'xtick',[1880:40:2000],'XTicklabel',[]); grid on

%%% repeat for scenario with pause in warming
load('Cascades_tau12_hiatus');
tplot = t - yrs_pre + date_start; % adjust to year trend starts
Qs = Qs*w; % multiply by nonminal width to get volume
Qs_eq = Qs_eq*w;
% subplot 221; hold on; 
% center12 = plot(tplot,L/1e3,'color',color1,'linewidth',lw);
% eq = plot(tplot,Leq/1e3,'-.k');

%subplot('position',[0.2 0.5 0.3 0.4]); hold on; 
axes('position',[0.7 0.45 0.28 0.28]); hold on; grid on
Q12=plot(tplot,Qs,'color',color1,'linewidth',lw-0.5);    
plot(tplot,Qs_eq,'-.','color',color1);

load('Cascades_tau48_hiatus');
Qs = Qs*w; % multiply by nonminal width to get volume
Qs_eq = Qs_eq*w;

Q48=plot(tplot,Qs,'color',color2,'linewidth',lw)    
plot(tplot,Qs_eq,'-.','color',color2);
eq=plot(0,0,'-.','color','k');

title('   Total melt (M_g)','FontWeight','normal')
ylabel('m^3 / yr')
xlim(plotlims)
ylim([1e7 2.2e7])
xlabel('Year CE')
set(gca,'xtick',[1880:40:2000])
box on

axes('position',[0.7 0.81 0.28 0.08]); 
plot(tplot,mu*(Tramp-Tramp(1)),'r','linewidth',lw);
title('specific melt anom.','FontWeight','normal')
ylabel('m / yr')
xlim([1860 2020])
ylim([0 1])
set(gca,'fontsize',8,'xtick',[1880:40:2000],'XTicklabel',[]); grid on
