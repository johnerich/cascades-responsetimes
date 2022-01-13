%%% script to cut a region from full RGI western CAN/US. Assumes you have 
%%% first loaded an RGI excel file into matlab

% for data, see:
% RGI Consortium (2017) Randolph glacier inventory?a dataset of global 
% glacier outlines: Version 6.0. Technical report, Global Land Ice Measurements 
% from Space, Colorado, USA (doi: https://doi.org/10.7265/N5-RGI-60)

clear all;
 
Nlim = 49;
Slim = 46; % between mt. Adams and Hood

Elim = -120;
Wlim = -122.5;

idx = find((CenLat>Slim).*(CenLat<Nlim).*(CenLon>Wlim).*(CenLon<Elim));

Area = Area(idx);
Aspect = Aspect(idx);
BgnDate = BgnDate(idx);
CenLat = CenLat(idx);
CenLon = CenLon(idx);
Connect = Connect(idx);
EndDate = EndDate(idx);
GLIMSId = GLIMSId(idx);
Lmax = Lmax(idx);
O1Region = O1Region(idx);
O2Region = O2Region(idx);
RGIId = RGIId(idx);
Slope = Slope(idx);
Status = Status(idx);
Zmax = Zmax(idx);
Zmed = Zmed(idx);
Zmin = Zmin(idx);
BgnDate = (BgnDate - 9999)/1e4;
