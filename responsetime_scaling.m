%%%
% script to estimate glacier response timescales from Haeberli and Hoelzle
% (1995) scaling, applied to RGI 6.0 data. Algorithm could be adjusted a
% bit based on assumptions about particular glacier geometry (e.g., is ELA
% really at median elevation), but key free parameters are really the basal
% shear stress (for thickness) and mass balance gradient (for bt).
% Parameters reflect those used for the Cascades analysis in Christian et
% al., 2021

% script also plots scatterplot and histograms of response times

% RGI data citation:
% RGI Consortium (2017) Randolph glacier inventory?a dataset of global glacier 
% outlines: Version 6.0. Technical report, Global Land Ice Measurements from 
% Space, Colorado, USA (doi: https://doi.org/10.7265/N5-RGI-60)

clear all; close all;
% load data. This file is a subsection of the RGI Western Canada/US region,
% already cut to encompass the WA cascades (46-49 N, 120-122.5 W)
 load('rgi6_WAcas_all.mat')
% Get names for well-known glaciers and groups:
names2index;

%% Haeberli and Hoezle Scaling:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants:
rho = 900;                  % ice density
g = 9.81;
dbdz = 6e-3;              % mass balance gradient, mwe per m; HH95 use 7.5e-3
dbdz = dbdz*(1000/rho);     % convert to ice equivalent
dbdx = 0.003;               % alternative horizontal mass balance gradient

f = 0.8;                    % shape factor
sb = 1.5e5*ones(size(Lmax));% constant basal shear stress:

Zmean = (Zmax + Zmin)/2;
delZ = Zmax - Zmin;

% threshold to ignore tiny glaciers:
Athresh = 0.1;                % area threshold, in km^2
dZthresh = 250;             % vertical span threshold, in m
% index for focusing on big glaciers [e.g., "histogram(tau(idx_big))" ]
idx_big = find((Area>Athresh).*(delZ>dZthresh));


% thickness part:
slopetot = atand((Zmax - Zmin)./Lmax);
hf = sb./(f*rho*g*sind(slopetot));  % avg flowline thickness
have = pi*hf/4;                     % avg total glacier thickness (following HH95) - note these assume a valley galcier geometry
hmax = 2.5*hf;                      % max thickness (following HH95)
Vtot = Area.*have/1e3;              % rough volume estimate in km^3

% assume ELA in middle of the glacier, then extrapolate mb gradient to terminus:
% bt = -dbdz*(Zmed-Zmin);    % for vertical gradient
bt = -dbdx*(Lmax/2);        % for horizontal gradient
tau = -hf./bt;              % response time

%% scatter H vs. bt
figure(1); 
color1 = [0    0.4470    0.7410];
subplot('position',[0.25 0.45 0.55 0.5]); hold on ; 
x = -linspace(0,1.1*max(abs(bt(idx_big))),100); x2 = repmat(x,100,1);
y = linspace(0,1.1*max(hf(idx_big)),100); y2 = repmat(y,100,1)';
tau_contours = -y2./x2;
ascale = 40;

[c h] = contour(x,y,tau_contours,[0:5:100]);
clabel(c,h,[0:10:100],'labelspacing',300);
scatter(bt(idx_big),hf(idx_big),abs(log(Area(idx_big)*ascale)),'k','filled')
scatter(bt(idx_big),hf(idx_big),Area(idx_big)*ascale,color1,'filled')

 %%% "key" for marker size
yy = 110;
xx = -14;
text(xx-1,yy+10,'Area (km^2)','fontweight','bold')
key0p5 = scatter(xx,yy-20,0.5*ascale,'k','filled'); text(xx+0.5,yy-20,'0.5')
key0p1 = scatter(xx,yy-10,2*ascale,'k','filled'); text(xx+0.5,yy-10,'2')
key0p02 = scatter(xx,yy-0,8*ascale,'k','filled'); text(xx+0.5,yy-0,'8')

alpha(0.7)
colormap(gray)
axis([x(end) -0.75 20 140])
ylabel('Ice thickness (m)')
xlabel('Terminus balance rate (m/yr)')
title('Individual glacier estimates (countours give \tau = -H/b_t (yrs))')

%%% uncomment to plot your favorite glacier!! see names2index script for
%%% options...
% scatter(bt(SouthCascade),hf(SouthCascade),Area(SouthCascade)*ascale,'r','filled')

%% very inelegant method for creating area-weighted distributions...

% bin EDGES (for making histograms)
bins_tau = [0:5:80];

% bin CENTERS (for plotting in non histogram mode)
centers_tau = [(bins_tau(1)+bins_tau(2))/2:bins_tau(2)-bins_tau(1):(bins_tau(end)+bins_tau(end-1))/2];

[tau_sort,tau_idx] = sort(tau(idx_big));

figure(2); % keep this in, Matlab for some reason needs a separate figure handle to retain histogram data. 
T = histogram(tau_sort,bins_tau);T = T.Values;

figure(1); 
subplot('position',[0.15 0.1 0.3 0.25]); hold on
percentlims = [0 35];
bar(centers_tau,100*T/length(idx_big),1,'facecolor',0.3*[1 1 1])
set(gca,'xtick',bins_tau(1:2:end))
ylabel('% of glaciers')
xlabel('\tau (years)')
ylim(percentlims)

close(figure(2))

A = Area(idx_big);          % original order of Volume vector
Atau(:) = A(tau_idx);       % sort Volume by increasing tau
A_at_tau = zeros(length(bins_tau)-1,1);

%%%% find total area in each bin range. T.Values is a vector describing
%%%% the number of glaciers in each bin. Since they are sorted by
%%%% increasing tau, we can then find total area in each
%%%% bin. 
% jj is index for A (ordered list of glacier volumes
% ii is bin index

% for response time:
jj = 1;
for ii = 1:(length(bins_tau) - 1);
    if T(ii) > 0
        A_at_tau(ii) = sum(Atau(jj:(jj+T(ii)-1)));
    else
        A_at_tau(ii) = 0;
    end
    jj = jj + T(ii);
end
jj = 1;


% plot "area-weighted" distributions
subplot('position',[0.6 0.1 0.3 0.25]); hold on ; 
bar(centers_tau,100*A_at_tau/sum(A),1,'facecolor',0.3*[1 1 1])
set(gca,'xtick',bins_tau(1:2:end))
ylabel('% of total Area')
xlabel('\tau (years)')
ylim(percentlims)

