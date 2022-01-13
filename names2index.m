% % glacier ID's for some well known glaciers. Returns indices so that
% attributes can be queried for specific glaciers by name [e.g., for
% response time, tau(SouthCascade) ]. Also creates some groups of indices,
% e.g. for the Volcanos. 
% This script is called in other scripts, so assumes that a set of
% glaciers has been loaded and 'RGIId' is an existing vector. 
% Therefore, the indices will depend on which subset of glaciers has been loaded. 

% If this isn't the case, load some set of glaciers:
%  load('rgi6_WAcas_all.mat')

% Olympics:
Blue =          find(strcmp(RGIId,'RGI60-02.14017'));
Hoh =           find(strcmp(RGIId,'RGI60-02.18763'));
%%% North Cascades:
SouthCascade =  find(strcmp(RGIId,'RGI60-02.18778'));
LeConte =       find(strcmp(RGIId,'RGI60-02.17902'));
NorthKlawatti = find(strcmp(RGIId,'RGI60-02.18077'));
Noisy =         find(strcmp(RGIId,'RGI60-02.17369'));
Sandalee =      find(strcmp(RGIId,'RGI60-02.17863'));
Silver =        find(strcmp(RGIId,'RGI60-02.17612'));
Boston =        find(strcmp(RGIId,'RGI60-02.18094'));
Inspiration =   find(strcmp(RGIId,'RGI60-02.18280'));
Lewis =         find(strcmp(RGIId,'RGI60-02.18200'));
Colonial =      find(strcmp(RGIId,'RGI60-02.17124'));
columbia =      find(strcmp(RGIId,'RGI60-02.18415'));
Daniels =       find(strcmp(RGIId,'RGI60-02.18487'));
Foss =          find(strcmp(RGIId,'RGI60-02.18337'));
honeycomb =     find(strcmp(RGIId,'RGI60-02.18486'))';
Iceworm =       find(strcmp(RGIId,'RGI60-02.18346'));
lowercurtis =   find(strcmp(RGIId,'RGI60-02.17307'));
Lyman =         find(strcmp(RGIId,'RGI60-02.18005'));
lynch =         find(strcmp(RGIId,'RGI60-02.18816'));
neve =          find(strcmp(RGIId,'RGI60-02.17731'));
yawning =       find(strcmp(RGIId,'RGI60-02.17806'));


%%% Mt. Baker:
Rainbow =       find(strcmp(RGIId,'RGI60-02.17733'));
Park =          find(strcmp(RGIId,'RGI60-02.17734'));
Roosevelt =     find(strcmp(RGIId,'RGI60-02.17735'));
Coleman =       find(strcmp(RGIId,'RGI60-02.17736'));
Thunder =       find(strcmp(RGIId,'RGI60-02.17737'));
Deming =        find(strcmp(RGIId,'RGI60-02.17738'));
Easton =        find(strcmp(RGIId,'RGI60-02.17739'));
Squak =         find(strcmp(RGIId,'RGI60-02.17740'));
Boulder =       find(strcmp(RGIId,'RGI60-02.17741'));
Talum =         find(strcmp(RGIId,'RGI60-02.17742'));
idx_baker = [Rainbow Park Roosevelt Coleman Thunder Deming Easton Squak Boulder Talum];
%%% Mt. Rainier:
Carbon =        find(strcmp(RGIId,'RGI60-02.14256'));
Winthrop =      find(strcmp(RGIId,'RGI60-02.14259'));
Emmons =        find(strcmp(RGIId,'RGI60-02.14297'));
Fryingpan =     find(strcmp(RGIId,'RGI60-02.14315'));
Whitman =       find(strcmp(RGIId,'RGI60-02.14335'));
Ingraham =      find(strcmp(RGIId,'RGI60-02.14333'));
Cowlitz =       find(strcmp(RGIId,'RGI60-02.18817'));
Nisqually =     find(strcmp(RGIId,'RGI60-02.14336'));
Kautz =         find(strcmp(RGIId,'RGI60-02.14341'));
SouthTahoma =   find(strcmp(RGIId,'RGI60-02.14352'));
Tahoma =        find(strcmp(RGIId,'RGI60-02.14323'));
Puyallup =      find(strcmp(RGIId,'RGI60-02.14329'));
SouthMowich =   find(strcmp(RGIId,'RGI60-02.14308'));
Edmunds =       find(strcmp(RGIId,'RGI60-02.14307'));
NorthMowich =   find(strcmp(RGIId,'RGI60-02.14293'));
Russell =       find(strcmp(RGIId,'RGI60-02.14277'));
idx_rainier = [Carbon Winthrop Emmons Fryingpan Whitman Ingraham Cowlitz Nisqually...
    Kautz SouthTahoma Tahoma Puyallup SouthMowich Edmunds NorthMowich Russell];
%%% Glacier Peak:
GP1 =           find(strcmp(RGIId,'RGI60-02.18282'));
Sitkum =        find(strcmp(RGIId,'RGI60-02.18369'));
Scimitar =      find(strcmp(RGIId,'RGI60-02.18475'));
Ptarmigan =     find(strcmp(RGIId,'RGI60-02.18476'));
Kennedy =       find(strcmp(RGIId,'RGI60-02.18477'));
Vista =         find(strcmp(RGIId,'RGI60-02.18478'));
Ermine =        find(strcmp(RGIId,'RGI60-02.18479'));
Dusty =         find(strcmp(RGIId,'RGI60-02.18480'));
NorthGuardian = find(strcmp(RGIId,'RGI60-02.18481'));
Chocolate =     find(strcmp(RGIId,'RGI60-02.18482'));
Cool =          find(strcmp(RGIId,'RGI60-02.18483'));
idx_GP = [GP1 Sitkum Scimitar Ptarmigan Kennedy Vista Ermine Dusty NorthGuardian Chocolate Cool];

%%% Mt. Adams:
Lava =          find(strcmp(RGIId,'RGI60-02.14462'));
Lyman =         find(strcmp(RGIId,'RGI60-02.17100'));
Wilson =        find(strcmp(RGIId,'RGI60-02.17101'));
Rusk =          find(strcmp(RGIId,'RGI60-02.17102'));
Klikitat =      find(strcmp(RGIId,'RGI60-02.17099'));
Mazama =        find(strcmp(RGIId,'RGI60-02.17103'));
WhiteSalmon =   find(strcmp(RGIId,'RGI60-02.14467'));
Pinnacle =      find(strcmp(RGIId,'RGI60-02.17098'));
Adams =         find(strcmp(RGIId,'RGI60-02.14600'));
idx_adams = [Lava Lyman Wilson Rusk Klikitat Mazama WhiteSalmon Pinnacle Adams];

%%% four major volcanos:
idx_vol = [idx_baker idx_GP idx_rainier idx_adams];
%%% BC:
Klinaklini =    find(strcmp(RGIId,'RGI60-02.05157'));
