
% The code implements the three-stage glacier model of Roe and Baker,
% JGlac, 2014 (RB14). You provide a time series of climate anomalies (either accumlation
% and melt-season temperature *or* annual-mean mass balance anomalies. The
% code then calculates the time series of the glacier response for the
% specified glacier parameters.


% Modified by John Erich Christian to include glacier melt season runoff contribution in response to a
% warming trend, for analyses in Christian et al. 2021, J. glac. 
% The given set of parameters give a small glacier with 12-year response
% time, and model is set up to run simple linear/step-wise warming trends.
% Variability and Precip anomalies can be added, too. 


clear all; %close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define ramp warming trend:
Ti = 0.5; % initial melt season temperature at maximum elevation (deg C)
%Tf = 2.5; % final temp 
dT1 = 1/0.8; % total temp change

yrs_pre = 100; % years of initial climate
yrs_trend = 060; % years of changing climate
yrs_post = 100; % years of new climate
Tramp1 = dT1*[zeros(yrs_pre,1);linspace(0,1,yrs_trend)';ones(yrs_post,1)] + Ti;

%2nd trend... just treated as an additional anomaly, same magnitude
hiatus = 30; % years of no warming
yrs_pre2 = yrs_pre + yrs_trend + hiatus;
yrs_trend2 = 40; % years of 2nd stage warming
yrs_post2 = length(Tramp1) - yrs_trend2 - yrs_pre2; % years of new climate
dT2 = 1/0.8;

% ramp forcing to be used
Tramp = Tramp1 + dT2*[zeros(yrs_pre2,1);linspace(0,1,yrs_trend2)';ones(yrs_post2,1)];
winsize = 1; % years of running mean
Tramp_smooth= conv(Tramp,ones(winsize,1)/winsize,'same');
Tramp(50:end-50) = Tramp_smooth(50:end-50);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define time array.
tf = length(Tramp)-1;          % length of integration [yrs]
ts = 0;             % starting year
dt = 1;             % time step [keep at 1yr always]
nts = (tf-ts)/dt+1; % number of time steps
t = ts:dt:tf;       % array of times [yr]

    %% glacier model parameters
    Zmax = 2800; % top of glacier
    Pbar = 4.5; % accumulation rate, m yr^-1
    mu = 0.8;      % melt factor [m yr^-1 K^-1]
    gamma = 6.5e-3; % assumed surface lapse rate [K m^-1] 
    tanphi = 0.35;   % assumed basal slope [no units]
    dbdz = mu*gamma;
    
    ela0 = Zmax - (Pbar/mu - Ti)/gamma;
    L0 = -2*(ela0 - Zmax)/tanphi;
    L0 = L0;
    w = 500;        % characteristic width of the glacier tongue [m]. 
    Atot = L0*w;  % total area of the glacier [m^2]
    ATgt0 = Atot;  % area of the glacier where some melting occurs [m^2]
    Aabl = Atot/2;  % ablation area [m^2] 
    dt = 1;         % incremental time step [yr]
    hbar = 50; % characteristic ice thickness near the terminus [m]


    % natural climate variability - for temperature and precipitation forcing
    sigP = 0;     % std. dev. of accumulation variability [m yr^-1]
    sigT = 0;     % std. dev. of melt-season temperature variability [m yr^-1]
    % natural climate variability - for mass balance forcing
    %sigb = 1;     % std. dev. of annual-mean mass balance [m yr^-1]

    % linear model coefficients, combined from above parameters
    % play with their values by choosing different numbers...
    alphaa = mu*ATgt0*dt/(w*hbar);
    beta = Atot*dt/(w*hbar);

    % glacier memory [ys]
    % this is the glacier response time  (i.e., memory) based on the above glacier geometry
    % if you like, just pick a different time scale to see what happens. 
    % Or also, use the simple, tau = hbar/b_term, if you know the terminus
    % balance rate from, e.g., observations
    tau = w*hbar/(mu*gamma*tanphi*Aabl);

    % coefficient needed in the model integrations
    % keep fixed - they are intrinsic to 3-stage model
    eps = 1/sqrt(3);
    phi = 1-dt/(eps*tau);

    %% note at this point you could create your own climate time series, using
    % random forcing, trends, oscillations etc.
    % Tp = array of melt-season temperature anomalise
    % Pp = array of accumlation anomalies
    % bp = array of mass balance anomalies
   
        % add noise, and re-define temp as anomaly, assuming glacier geometry corresponds to initial Temp
        Tp = sigT*randn(nts,1);
        Tp = Tp + Tramp-Tramp(1); 
        
        Pp = sigP*randn(nts,1);
    %   bp = sigb*randn(nts,1);
        %load('pnw_synthetic'); Pp = Pp_pnw; Tp = Tp_pnw + Tramp-Tramp(1);

        %% integrate the 3 stage model equations forward in time
        L3s = zeros(size(Tramp));
        for i = 4:nts
            L3s(i) = 3*phi*L3s(i-1)-3*phi^2*L3s(i-2)+1*phi^3*L3s(i-3)...
             + dt^3*tau/(eps*tau)^3 * (beta*Pp(i) - alphaa*Tp(i));
        % if you want to use mass balance anomalies instead comment out the 2 lines
        % above, and uncomment the 2 lines below
        % L3s(i) = 3*phi*L3s(i-1)-3*phi^2*L3s(i-2)+1*phi^3*L3s(i-3)...
        %        + dt^3*tau/(eps*tau)^3 * (beta*bp);
        end

 %% Calculate instantaneous equilibrium length and degree of disequilibrium (dod) - JEC
 % The instantaneous eq. length is length is where the glacier would be
 % were it allowed to equilibrate with the instantaneous equilibrium.
 % Disequilibrium is the difference between that and the transient length. 
 % note this will be super noisy if variability in T or P is included.
 % recommend not plotting in that case (or could plot a running mean)

T = Tramp(1) + Tp; % add back in initial temp at zmax
P = Pbar + Pp;
ela = Zmax - (P/mu - T)/gamma;
Leq = -2*(ela - Zmax)/tanphi;
L = L3s + L0;
dod = L - Leq;

% mass balance anomaly for plotting:
% point MB:
bp = Pp - mu*Tp;
bramp = -mu*Tp;
winsize = 30;
bpsmooth = conv(bp,ones(winsize,1)/winsize,'same'); % running mean

% 3-stage model net mass balance calculation (see Roe et al., 2021, TC):
% this is the traditional, glacier-averaged mass balance that incorporates
% evolving glacier geometry:
bnet = bp - L3s/(beta*tau);

%% Calculate melt-season runoff
% Qs is simply the summer melt (mu*T) integrated over the evolving glacier
% surface. With the model assumptions here (linear mass balance gradient,
% constant width), this simplifies to the avg. between melt at the terminus
% and melt at the head of the glacier, multiplied by the instantaneous
% glacier area. One could, in principle, numerically integrate over an
% arbitrary geometry, however. The point here is to capture the competeing
% tendecies of more melt per area vs. less area to melt. See Christian et
% al., 2021 for discussion. 

Qs = w*0.5*mu*L.*(2*T + gamma*tanphi*L);


% melt flux that would occur over equilibrium geometry:
Qs_eq = w*0.5*mu*Leq.*(2*T + gamma*tanphi*Leq);

%% Plot it up:

plotlims = [-50 300];
tplot = t - 100; % adjust to year trend starts
date_start = 1880; % start of trend
date_end = 2020; % end of trend
tplot = t - yrs_pre + date_start; % adjust to year trend starts
plotlims = [date_start-20 date_end];
lw = 1;

color1 = [1 0.5 0.1];
color2 = [0.2 0.6 1];
color3 = [0.5 0.8 1];

figure(1);
subplot 221; hold on; grid on
plot(tplot,bp,'r','linewidth',lw-0.5);
plot(tplot,bpsmooth,'r','linewidth',lw+0.5);
plot(tplot,bnet,'--','color',[1 0.5 0.5],'linewidth',lw+0.5);
legend('Point MB','Point (smoothed)','Net MB','Location','best')
ylabel('m yr^{-1}')
xlabel('Time (yrs)')
xlim(plotlims)
title('mass balance forcing')

subplot 222; hold on; grid on
center = plot(tplot,L/1e3,'color',color1,'linewidth',lw);
if sigP == 0 && sigT ==0
    eq = plot(tplot,Leq/1e3,'-.k');
    legend([eq center],'equilibrium','3-stage model')
end
ylabel('Glacier Length (km)')
xlabel('Time (yrs)')
xlim(plotlims)
title('Glacier response')
%%
subplot 223; hold on; grid on
plot(tplot,dod/1e3,'color',color1,'linewidth',lw);
ylabel('Diseq. (km)')
xlabel('Time (yrs)')
xlim(plotlims)
title('Length disequilibrium')

%%
subplot 224; hold on; grid on
plot(tplot,Qs,'color',color3,'linewidth',lw-0.5)
plot(tplot,Qs_eq,'-.k')
title('Melt-season integrated melt flux')
legend('Transient','Equilib.')
ylabel('m^3 / yr')
xlabel('Time (yrs)')



 

