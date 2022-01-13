clear all; %close all;

% Fig 1: Schematic, adapted from Christian et al., J. Glac., 2018
% timeseries are ouput from 3stage model, forced with simple ramp
% warming

load('idealized_retreat_3stage')
dod_3s = L_3s - L_eq; % absolute disequilibrium
frac3s = L_3s./L_eq; frac3s(isnan(frac3s))=0; % fractional equilb.

years = [0:10:690] - 100;
yrlim = 120;

colorL = [173,216,230]/255;
colordod = [46,139,87]/255;
colorfrac = [144,238,144]/255;
coloreq = [72,61,139]/255;
fontsz = 9;
figure(1); clf

% length plot:
subplot('position',[0.1 0.45 0.45 0.5]); hold on; grid off
plot(years,L_3s,'color',colorL,'linewidth',3);
plot(years,L_eq,':','color',coloreq,'linewidth',2);
legend('Transient response','Equilibrium response')
xlim([-10 yrlim]); xlabel('Time (a)')
ylim([-2 0.1]); ylabel('Length anomaly (km)')
set(gca,'fontsize',fontsz,'ytick',[-2:0.5:0])

%%% panel b, glacier diseq:
subplot('position',[0.1 0.1 0.18 0.25]); hold on; 
plot(years,dod_3s,'color',colordod,'linewidth',3)
xlim([-10 yrlim]); xlabel('Time (a)')
ylim([-0.02 0.7]); ylabel('Disequilibrium (km)')
set(gca,'fontsize',fontsz,'ytick',[0:0.2:1])

% fractional equilibration
subplot('position',[0.37 0.1 0.18 0.25]); hold on; 
plot(years,frac3s,'color',colorfrac,'linewidth',3);
ylim([0 1.01]); ylabel('Fractional equilib.')
xlim([0 yrlim]); xlabel('Time (a)')
set(gca,'fontsize',fontsz)
