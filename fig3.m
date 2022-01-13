clear all; close all

% Choose output to plot - reponse times based on horizontal or vertical MB
% gradients:
load('final_scaling_estimates_dbdx')
%load('final_scaling_estimates_dbdz')
names2index;

% panel A: plot fractional equilibration for whole collection of glaciers
figure(1); 

eps = 1/sqrt(3);
trend_yrs = 1880:1:2018;        % assumed duration of ramp mass balance trend
t = trend_yrs - trend_yrs(1);
t = repmat(t,length(idx_big),1);
tau_big = repmat(tau(idx_big),1,length(trend_yrs));

frac_eq = 1 - (3*eps*tau_big./t).*(1-exp(-t./(eps*tau_big))) + exp(-t./(eps*tau_big)).*(t./(2*eps*tau_big) + 2);

subplot('position',[0.08 0.4 0.3 0.4]); box off; hold on

gray=plot(trend_yrs,frac_eq,'color',0.7*[1 1 1],'linewidth',0.5);
ylabel('Fractional equilibration (L''/L''_{eq})')
xlabel('Year')
xlim([1880 2020]); 
set(gca,'xtick',[1880:40:2000],'fontsize',10)



% panel b: plot distribution of current fractional equilibration
frac_today = frac_eq(:,end);   
bins_frac = [0:0.05:1];
centers_frac = [(bins_frac(1)+bins_frac(2))/2:bins_frac(2)-bins_frac(1):(bins_frac(end)+bins_frac(end-1))/2];
figure(2)
F = histogram(frac_today,bins_frac);F = F.Values;
close(figure(2))

figure(1);subplot('position',[0.44 0.64 0.2 0.16]); hold on
bar(centers_frac,100*F/length(idx_big),'facecolor',0.3*[1 1 1])
ylabel('% of glaciers')
xlabel('Fractional equilibration')
set(gca,'xtick',bins_frac(1:5:end),'fontsize',8)
ylim([0 15])

% panel c: Nisqually and South Cascade
% Load length data for Nisqually and S. Cascade glaciers, originally published in: 
% Leclercq PW, Oerlemans J, Basagic HJ, Bushueva I, Cook A and Le Bris R (2014) 
% A data set of worldwide glacier fluctuations. The Cryosphere, 8, 659?672 
% (doi: 10.5194/tc-8-659-2014)
load ('Lengths_Nis_sc.mat')

trend_yrs = 1880:2018;
t = 1:length(trend_yrs);
eps = 1/sqrt(3);

% nisqually
tau_nis = tau(Nisqually);
tau_nis_min = 48/9;
tau_nis_max = 122/9;

frac3_nis = 1 - (3*eps*tau_nis./t).*(1-exp(-t/(eps*tau_nis))) + exp(-t/(eps*tau_nis)).*(t/(2*eps*tau_nis) + 2);
frac3_nis_min = 1 - (3*eps*tau_nis_min./t).*(1-exp(-t/(eps*tau_nis_min))) + exp(-t/(eps*tau_nis_min)).*(t/(2*eps*tau_nis_min) + 2);
frac3_nis_max = 1 - (3*eps*tau_nis_max./t).*(1-exp(-t/(eps*tau_nis_max))) + exp(-t/(eps*tau_nis_max)).*(t/(2*eps*tau_nis_max) + 2);
 
% south cascade
tau_sc = tau(SouthCascade);
tau_sc_min = 99/5;
tau_sc_max = 203/5;

frac3_sc = 1 - (3*eps*tau_sc./t).*(1-exp(-t/(eps*tau_sc))) + exp(-t/(eps*tau_sc)).*(t/(2*eps*tau_sc) + 2);
frac3_sc_min = 1 - (3*eps*tau_sc_min./t).*(1-exp(-t/(eps*tau_sc_min))) + exp(-t/(eps*tau_sc_min)).*(t/(2*eps*tau_sc_min) + 2);
frac3_sc_max = 1 - (3*eps*tau_sc_max./t).*(1-exp(-t/(eps*tau_sc_max))) + exp(-t/(eps*tau_sc_max)).*(t/(2*eps*tau_sc_max) + 2);


  
subplot('position',[0.46 0.4 0.16 0.16]); hold on
color1 = [0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];

nis=plot(trend_yrs,frac3_nis,'color',color1,'linewidth',1); 
patch([trend_yrs, flip(trend_yrs)],[frac3_nis_min, flip(frac3_nis_max)],color1,'edgecolor','none'); alpha(0.2)
patch([trend_yrs, flip(trend_yrs)],[frac3_sc_min, flip(frac3_sc_max)],color2,'edgecolor','none'); alpha(0.2)
sc=plot(trend_yrs,frac3_sc,'color',color2,'linewidth',1); alpha(0.5)
ylabel('L''/L''_{eq}')
legend([nis sc],'Nisqually, \tau \sim 5-14 yr','South Cascade, \tau \sim 19-41 yr')
xlabel('Year')
xlim([1880 2020]); 
set(gca,'xtick',[1880:40:2000],'fontsize',8)


%%% panel d: length plot:
% nisqually

dL_nis = 2174/1e3; % total length change over ~20th century period
L_Nisq(:,2) = L_Nisq(:,2)/1e3;
years_nis = 1880:2001;
t = years_nis(end) - years_nis(1)+1;
frac3_nis = 1 - (3*eps*tau_nis./t).*(1-exp(-t/(eps*tau_nis))) + exp(-t/(eps*tau_nis)).*(t/(2*eps*tau_nis) + 2);
frac3_nis_min = 1 - (3*eps*tau_nis_min./t).*(1-exp(-t/(eps*tau_nis_min))) + exp(-t/(eps*tau_nis_min)).*(t/(2*eps*tau_nis_min) + 2);
frac3_nis_max = 1 - (3*eps*tau_nis_max./t).*(1-exp(-t/(eps*tau_nis_max))) + exp(-t/(eps*tau_nis_max)).*(t/(2*eps*tau_nis_max) + 2);
 
Lcom_nis = dL_nis/frac3_nis(end) - dL_nis;
Lcom_nis_bounds = dL_nis./[frac3_nis_min frac3_nis_max]-dL_nis;


% south cascade

dL_sc = 1802/1e3; % total length change over ~20th century period
years_sc = 1890:2007;
L_Scg(:,2) = L_Scg(:,2)/1e3;
t = years_sc(end) - years_sc(1)+1;

frac3_sc = 1 - (3*eps*tau_sc./t).*(1-exp(-t/(eps*tau_sc))) + exp(-t/(eps*tau_sc)).*(t/(2*eps*tau_sc) + 2);
frac3_sc_min = 1 - (3*eps*tau_sc_min./t).*(1-exp(-t/(eps*tau_sc_min))) + exp(-t/(eps*tau_sc_min)).*(t/(2*eps*tau_sc_min) + 2);
frac3_sc_max = 1 - (3*eps*tau_sc_max./t).*(1-exp(-t/(eps*tau_sc_max))) + exp(-t/(eps*tau_sc_max)).*(t/(2*eps*tau_sc_max) + 2);
    
Lcom_sc = dL_sc/frac3_sc - dL_sc;
Lcom_sc_bounds = dL_sc./[frac3_sc_min frac3_sc_max]-dL_sc;

subplot('position',[0.72 0.4 0.25 0.4]); hold on
sc=plot(L_Scg(:,1),L_Scg(:,2),'-o','color',color2,'markerfacecolor',color2,'markersize',3);
nis=plot(L_Nisq(:,1),L_Nisq(:,2),'-o','color',color1,'markerfacecolor',color1,'markersize',3);

plot(L_Nisq(end,1),L_Nisq(end,2)-Lcom_nis,'+','color',color1)
plot(L_Scg(end,1),L_Scg(end,2)-Lcom_sc,'+','color',color2)
committed=plot(0,0,':','color',0.5*[1 1 1],'linewidth',2);
scaling=plot(0,0,'+','color',0.5*[1 1 1]);

plot([L_Nisq(end,1) L_Nisq(end,1)],L_Nisq(end,2)-Lcom_nis_bounds,':','color',color1,'linewidth',2)
plot([L_Scg(end,1) L_Scg(end,1)],L_Scg(end,2)-Lcom_sc_bounds,':','color',color2,'linewidth',2)
ylabel('Length anomaly (km)')
xlabel('Year CE')
legend([sc nis committed scaling],'South Cascade','Nisqually','committed retreat (upper/lower bounds)','scaling estimate')
xlim([1880 2020]);
ylim([-4200 100]/1e3)
set(gca,'ytick',[-5000:1000:0]/1e3)
set(gca,'xtick',[1880:40:2000],'fontsize',10)
grid on
