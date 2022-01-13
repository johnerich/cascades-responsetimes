% script to compare length disequilibrium and length variability 
clear all;

sigb = 1; % m/yr

color1 = [1 0.5 0.1];
color2 = [0.2 0.6 1];
color3 = 0.2*[1 1 1];
load('Cascades_tau12_schematic')
% panel A: climate:
%axes('position',[0.21 0.88 0.15 0.12]); hold on; grid on
subplot 231; hold on;
plot(tplot,bp,'color',[1 0.6 0.6]);
plot(tplot,bramp,'--k','linewidth',lw);
plot(tplot,bpsmooth,'r','linewidth',lw);
legend('trend + noise','forced trend','30-yr mean','location','best')
ylim([-4 3]);
set(gca,'xtick',[1860:40:2020],'fontsize',9)
ylabel('b'' (m yr^{-1})')
xlabel('Year CE')
xlim([1850 2050])
%title('a)')
text(1860,-3.2,'SNR_b = 1 century^{-1}','FontSize',9)

% Panel B: glacier response
subplot 232; hold on
plot([2020 2020],[0 5],'-','color',0.9*[1 1 1])

noisy12 = noisy_retreat_SNR1percentury+L(1);
Leq12 = Leq;
eps = 1/sqrt(3);
kap = 1 - dt*(eps.*tau).^(-1);
psi = tau.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
sigL12 = beta.*psi*sigb;
D12 = L(241)-Leq(241);

tplot = t - yrs_pre + date_start; % adjust to year trend starts
patch([tplot,flip(tplot)],[(L+sigL12)/1e3;flip(L-sigL12)/1e3],color1,'edgecolor','none');alpha(0.2);
eq = plot(tplot,Leq/1e3,'--','color',color3);
center = plot(tplot,L/1e3,':','color',color1);
L12=plot(tplot,noisy12/1e3,'color',color1,'linewidth',lw);
plot(2020,L(241)/1e3,'*','color',color1);

load('Cascades_tau48_schematic')
noisy48 = noisy_retreat_SNR1percentury+L(1);
Leq48 = Leq;
kap = 1 - dt*(eps.*tau).^(-1);
psi = tau.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);
sigL48 = beta.*psi*sigb;
D48 = L(241)-Leq(241);

tplot = t - yrs_pre + date_start; % adjust to year trend starts
patch([tplot,flip(tplot)],[(L+sigL48)/1e3;flip(L-sigL48)/1e3],color2,'edgecolor','none');alpha(0.2);
eq = plot(tplot,Leq/1e3,'--','color',color3);
center48 = plot(tplot,L/1e3,':','color',color2);
L48=plot(tplot,noisy48/1e3,'color',color2,'linewidth',lw);
plot(2020,L(241)/1e3,'*','color',color2)
legend([L12 L48 eq],'\tau = 12 yrs','\tau = 48 yrs','Equilibrium length','Location','best')
plotlims = [1850 2050];
ylabel('Glacier Length (km)')
xlabel('Year CE')
xlim(plotlims)
set(gca,'ytick',[0:0.5:4.5],'xtick',[1860:40:2020],'fontsize',9)
ylim([1 4.8])
text(1860,2.1,'SNR_b = 1 century^{-1}','FontSize',9)
%title('b)')

%%
% Panel B: sigmaL vs. Disequilibrium by response time
% glacier parameters:
beta = 100; % this can be arbitrary - it cancels in the ratio
tau = [5:100]; 
dt = 1;


eps = 1/sqrt(3);
kap = 1 - dt*(eps.*tau).^(-1);
psi = tau.*(((1-kap).*(1+4.*kap.^2 + kap.^4))./(1 + kap).^5).^(0.5);


% climate parameters:
sigb = 1; % m/yr
bdot = 0.01; % m/yr /yr
t = 140*ones(size(tau)); % yrs of trend
snr_century = bdot*100/sigb;
diseq = tau.*beta*bdot.*t .* (3*eps*tau./t.*(1 - exp(-t./(eps*tau))) - exp(-t./(eps*tau)).*(t./(2*eps*tau) + 2));
diseq_lim = 3*eps*tau.^2.*beta.*bdot; % limit as t >> tau
sigL = beta.*psi*sigb;

ratio = diseq./sigL;
ratio_lim = diseq_lim./sigL;
ratio_lim_plus = diseq_lim./(sigL*2);
ratio_lim_minus = diseq_lim./(sigL*0.5);

subplot 233; hold on; grid on;
%plot(tau,ratio,'linewidth',2);
center=plot(tau,ratio_lim,'linewidth',2,'color',0.8*[0.9 0.9 1]);
plus=plot(tau,ratio_lim_plus,'-.','linewidth',1,'color',0.8*[0.9 0.9 1]);
min=plot(tau,ratio_lim_minus,'--','linewidth',1,'color',0.8*[0.9 0.9 1]);
plot(12,D12/sigL12,'*','color',color1)
plot(48,D48/sigL48,'*','color',color2)

xlabel('Response time (yrs)');
ylabel('Disequilibrium / \sigma_L')
ylim([0 12]);
xlim([0 60])
set(gca,'ytick',[1:2:15],'xtick',[0:10:100],'fontsize',9)
legend([min center plus],'SNR_b = 2 century^{-1}','SNR_b = 1 century^{-1}','SNR_b = 0.5 century^{-1}','location','northoutside')
set(legend,'fontsize',7)
%title('c)')

%% part two: forcing with hiatus

% panel C:
color1 = [0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];

load('./TCascades_1800-2018.mat','yrs','BE_aprsep','BE_octmar')
load('./PCascades_1900-2017','years','octmar_cut')

winsize = 30;

T = BE_aprsep(81:end);yrsT = yrs(81:end);
Tsmooth = conv(T,ones(1,winsize)/winsize,'same');
P = octmar_cut*6;
Psmooth = conv(P,ones(1,winsize)/winsize,'same');

subplot 223; hold on

[hax line1 line2] = plotyy(years,P,yrsT,T);

ylabel(hax(2),'Apr-Sep Temp. (\circC)')
ylim(hax(2),[-4 2]);
set(hax(2),'ytick',[-1:3])
%line2.Color = [1 0.6 0.6];

ylabel(hax(1),'Oct-Mar Precip. (cm)')
ylim(hax(1),[80 300]);
set(hax(1),'ytick',[100:25:200])
%line1.Color = [0.6 0.6 1];

plot(hax(1),years(winsize/2:end-winsize/2),Psmooth(winsize/2:end-winsize/2),'color',color1,'linewidth',1);
hold(hax(2))
plot(hax(2),yrsT(winsize/2:end-winsize/2),Tsmooth(winsize/2:end-winsize/2),'color',color2,'linewidth',1);

xlim(hax(1),[1880 2020])
xlim(hax(2),[1880 2020])

set(hax(1),'FontSize',9,'xtick',[1860:40:2020])
set(hax(2),'FontSize',9,'xtick',[1860:40:2020])

xlabel('Year CE')
%title('d)')
%%

% panel D:
load('Cascades_tau12_hiatus');
tplot = t - yrs_pre + date_start; % adjust to year trend starts

subplot 224; hold on; 
center12 = plot(tplot,L/1e3,'color',color1,'linewidth',lw);
eq = plot(tplot,Leq/1e3,'--','color',color3);

load('Cascades_tau48_hiatus');

center48 = plot(tplot,L/1e3,'color',color2,'linewidth',lw);
eq = plot(tplot,Leq/1e3,'--','color',color3);
%legend([eq center12 center48],'equilib. response','\tau = 12 yr','\tau = 48 yr','location','southoutside','orientation','horizontal')
ylabel('Glacier Length (km)')
xlabel('Year CE')
xlim([1850 2020])
ylim([2.3 4.6]);
set(gca,'ytick',[0:0.5:4.5],'xtick',[1860:40:2020],'fontsize',9)
%title('e)')

legend([center12 center48 eq],'\tau = 12 yrs','\tau = 48 yrs','Equilibrium length','Location','best')


axes('position',[0.82 0.4 0.08 0.07]); 
plot(tplot,-mu*(Tramp-Tramp(1)),'r','linewidth',lw);
ylabel('b'' (m yr^{-1})')
xlim([1880 2020])
%ylim([0 1.1])
set(gca,'fontsize',7,'xtick',[1880:60:2000],'ytick',-1:0.5:0); grid on


