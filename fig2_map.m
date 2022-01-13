clear all
% Choose output to plot - reponse times based on horizontal or vertical MB
% gradients:
load('final_scaling_estimates_dbdx')
%load('final_scaling_estimates_dbdz')

% DEM basemap from SRTM
% citation: Jarvis A, Reuter H, Nelson A and Guevara E (2008) Hole-filled 
% SRTM for the globe Version 4, available from the CGIAR-CSI SRTM 90m Database 
% (http://srtm.csi.cgiar.org). Accessed 29 Oct, 2020
% load('srtm_12_03.mat')

figure(1); clf; hold on;

ax1 = axes('position',[0.2 0.15 0.45 0.75]);
% imagesc(lon,lat,flipud(dem)); axis xy
xlim([-122 -120.5])
ylim([46 49]);

ax2 = axes('position',[0.2 0.15 0.45 0.75]);
ascale = 30;
scatter(CenLon(idx_big),CenLat(idx_big),Area(idx_big)*ascale,tau(idx_big),'filled');alpha(0.7)
xlim([-123 -120])
linkaxes([ax1,ax2]);

ax2.YTickLabel = [];
ax2.XTickLabel = [];
ax2.Visible = 'off'

colormap(ax1,'pink');caxis(ax1,[-2000 2000]);
cb1 = colorbar(ax1,'position',[.15 .9 .75 .02],'orientation','horizontal');
colormap(ax2,1*flip(winter));caxis(ax2,[5 60]);
cb2 = colorbar(ax2,'Position',[.15 .05 .75 .02],'orientation','horizontal');
title(cb1,'Elevation (m)');
title(cb2,'Glacier response time (years)','Fontsize',14)
set(ax2,'YTickLabel',[],'XTickLabel',[]);
set([ax1,ax2],'Position',[.15 .14 .75 .75]);
axis equal
xlabel(ax1,'Longitude (^{\circ}E)')
ylabel(ax1,'Latitude (^{\circ}N)')

%Area key:
 yy = 47;
 xx = -121;
hold on
kcolor = [87 7 6]/255; %deep maroon
text(xx,yy,'Area ','fontsize',16)
key0p5 = scatter(xx,yy-0.2,8*ascale,'markeredgecolor',kcolor); text(xx+0.1,yy-0.2,'8 km^2','color',kcolor)
key0p1 = scatter(xx,yy-0.4,2*ascale,'markeredgecolor',kcolor); text(xx+0.1,yy-0.4,'2 km^2','color',kcolor)
key0p02 = scatter(xx,yy-0.6,0.5*ascale,'markeredgecolor',kcolor); text(xx+0.1,yy-0.6,'0.5 km^2','color',kcolor)
