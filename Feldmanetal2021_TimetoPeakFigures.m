clc
clear

load coast   
clat = lat;                                                        
clon = long ; 

Ncr = 200;
Ncb = 200;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [blu;flipud(red)];

Regime4Flag_Africa_Full = [];
TTPSeas_Full = [];
dVWC_Full = [];
TTPThreshold_Full = [];
RootzoneDiff_Full = [];
upperBoundMat_Full = [];
ZerosDiffMat_Full = [];
TTPZeros_Full = [];
DrydownCharMat_Full = [];
ANOVAAntRZSM_Full = [];
ANOVAAntSM_Full = [];
ANOVAAntVOD_Full = [];
ANOVASMPulse_Full = [];
TTPOverpassNumMat_Full = [];
TTPLongShortFrac_Full = [];
TTPCorrSMAnt_Full = [];
TTPCorrSMAnt0_Full = [];
TTPCorrSMPulse_Full = [];
TTPCorrSMPulse0_Full = [];
TTPCorrSMRootAnt_Full = [];
TTPCorrSMRootAnt0_Full = [];
TTPCorrVODAnt_Full = [];
TTPCorrVODAnt0_Full = [];
VWCAntpulsefraction_Full = [];

cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy_TauEfolding/Figures_ToDistribute')
dr = cat(1,dir('*TimeToPeak_Globe_V18*.mat'));
load('GlobeAncillaryHalfDegree')
load('GlobeAncillaryHalfDegreeTwoYearGPM')
PrecipMat=PrecipMat./2;

for j = 1:length(dr)
load(dr(j).name)    

Regime4Flag_Africa_Full = ...
    cat(4,Regime4Flag_Africa_Full,Regime4Flag_Africa);
TTPZeros_Full = ...
    cat(4,TTPZeros_Full,TTPZeros);
upperBoundMat_Full = ...
    cat(4,upperBoundMat_Full,upperBoundMat);
ZerosDiffMat_Full = ...
    cat(4,ZerosDiffMat_Full,ZerosDiffMat);
DrydownCharMat_Full = ...
    cat(4,DrydownCharMat_Full,DrydownCharMat);
TTPOverpassNumMat_Full = ...
    cat(4,TTPOverpassNumMat_Full,TTPOverpassNumMat);
TTPLongShortFrac_Full = ...
    cat(4,TTPLongShortFrac_Full,TTPLongShortFrac);

TTPCorrSMAnt_Full = ...
    cat(4,TTPCorrSMAnt_Full,TTPCorrSMAnt);
TTPCorrSMAnt0_Full = ...
    cat(4,TTPCorrSMAnt0_Full,TTPCorrSMAnt0);
TTPCorrSMPulse_Full = ...
    cat(4,TTPCorrSMPulse_Full,TTPCorrSMPulse);
TTPCorrSMPulse0_Full = ...
    cat(4,TTPCorrSMPulse0_Full,TTPCorrSMPulse0);
TTPCorrSMRootAnt_Full = ...
    cat(4,TTPCorrSMRootAnt_Full,TTPCorrSMRootAnt);
TTPCorrSMRootAnt0_Full = ...
    cat(4,TTPCorrSMRootAnt0_Full,TTPCorrSMRootAnt0);
TTPCorrVODAnt_Full = ...
    cat(4,TTPCorrVODAnt_Full,TTPCorrVODAnt);
TTPCorrVODAnt0_Full = ...
    cat(4,TTPCorrVODAnt0_Full,TTPCorrVODAnt0);

ANOVAAntRZSM_Full = ...
    cat(4,ANOVAAntRZSM_Full,ANOVAAntRZSM);
ANOVAAntVOD_Full = ...
    cat(4,ANOVAAntVOD_Full,ANOVAAntVOD);
ANOVAAntSM_Full = ...
    cat(4,ANOVAAntSM_Full,ANOVAAntSM);
ANOVASMPulse_Full = ...
    cat(4,ANOVASMPulse_Full,ANOVASMPulse);

VWCAntpulsefraction_Full = ...
    cat(4,VWCAntpulsefraction_Full,VWCAntpulsefraction);
end

Regime4Flag_Africa = nanmean(Regime4Flag_Africa_Full,4);
TTPZeros = nanmean(TTPZeros_Full,4);
upperBoundMat = nanmean(upperBoundMat_Full,4);
ZerosDiffMat = nanmean(ZerosDiffMat_Full,4);
DrydownCharMat = nanmean(DrydownCharMat_Full,4);
TTPOverpassNumMat = nanmean(TTPOverpassNumMat_Full,4);
TTPLongShortFrac = nanmean(TTPLongShortFrac_Full,4);
VWCAntpulsefraction = nanmean(VWCAntpulsefraction_Full,4);

TTPCorrSMAnt = nanmean(TTPCorrSMAnt_Full,4);
TTPCorrSMAnt0 = nanmean(TTPCorrSMAnt0_Full,4);
TTPCorrSMPulse = nanmean(TTPCorrSMPulse_Full,4);
TTPCorrSMPulse0 = nanmean(TTPCorrSMPulse0_Full,4);
TTPCorrSMRootAnt = nanmean(TTPCorrSMRootAnt_Full,4);
TTPCorrSMRootAnt0 = nanmean(TTPCorrSMRootAnt0_Full,4);
TTPCorrVODAnt = nanmean(TTPCorrVODAnt_Full,4);
TTPCorrVODAnt0 = nanmean(TTPCorrVODAnt0_Full,4);

ANOVAAntRZSM = nanmean(ANOVAAntRZSM_Full,4);
ANOVAAntVOD = nanmean(ANOVAAntVOD_Full,4);
ANOVAAntSM = nanmean(ANOVAAntSM_Full,4);
ANOVASMPulse = nanmean(ANOVASMPulse_Full,4);

dll = 0.5;   %aggregate to dll degrees 
alat = ( -60 : dll : 60 );                              
alon = (-179 : dll : 179 );                             
nlat = length(alat);                                    
nlon = length(alon); 
alat = flipud(transpose(repmat(alat,nlon,1))); %lat matrix         
alon = repmat(alon,nlat,1); 

[~,~,~,~,...
   TopLeftCornerLat,~,~,~] = InterpEdgesCorners(alat);
TopLeftCornerLat(:,1) = TopLeftCornerLat(:,2);
TopLeftCornerLat(1,:) = TopLeftCornerLat(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLon,~,~,~] = InterpEdgesCorners(alon);
TopLeftCornerLon(1,:) = TopLeftCornerLon(2,:);
TopLeftCornerLon(:,1) = TopLeftCornerLon(:,2)-0.5;

Regime4Flag = Regime4Flag_Africa(:,:,1);
Regime4Flag(Regime4Flag>0.03)=nan;
Regime4Flag(Regime4Flag<=0.03)=1;

Regime4FlagA = Regime4Flag_Africa(:,:,1);
Regime4FlagA(Regime4FlagA>0.03)=1;
Regime4FlagA(Regime4FlagA<=0.03)=nan;

SMStDFlag = Regime4Flag_Africa(:,:,1);
SMStDFlag(SMStDFlag>0.03)=1;
SMStDFlag(SMStDFlag<=0.03)=nan;

%%%%% Bare and forest flags %%%%%%
Bareflag = IGBPMat(:,:,16);
Bareflag(Bareflag>0.75)=nan;
Bareflag(isfinite(Bareflag))=1;

Frozenflag = IGBPMat(:,:,15);
Frozenflag(Frozenflag>0.75)=nan;
Frozenflag(isfinite(Frozenflag))=1;

Regime4Flag_Africa = Regime4Flag_Africa.*Bareflag.*Frozenflag;
ZerosDiffMat = ZerosDiffMat.*Bareflag.*Frozenflag;
upperBoundMat = upperBoundMat.*Bareflag.*Frozenflag;
DrydownCharMat = DrydownCharMat.*Bareflag.*Frozenflag;
TTPZeros = TTPZeros.*Bareflag.*Frozenflag;
TTPOverpassNumMat = TTPOverpassNumMat.*Bareflag.*Frozenflag;
TTPLongShortFrac = TTPLongShortFrac.*Bareflag.*Frozenflag;
VWCAntpulsefraction = VWCAntpulsefraction.*Bareflag.*Frozenflag;

ANOVAAntRZSM = ANOVAAntRZSM.*Bareflag.*Frozenflag;
ANOVAAntVOD = ANOVAAntVOD.*Bareflag.*Frozenflag;
ANOVAAntSM = ANOVAAntSM.*Bareflag.*Frozenflag;
ANOVASMPulse = ANOVASMPulse.*Bareflag.*Frozenflag;

TTPFlag = TTPZeros(:,:,7);
TTPFlag(TTPFlag>0.5)=1;
TTPFlag(TTPFlag<=0.5)=nan;

% Calculate percentages for bins 2 and 3
AntRZSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
AntVODCount = nan(size(ANOVAAntRZSM(:,:,7)));
AntSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
SMPulseCount = nan(size(ANOVAAntRZSM(:,:,7)));

AntRZSMCount(TTPZeros(:,:,7)>0.5 & ANOVAAntRZSM(:,:,7)<0.05)=1;
AntVODCount(TTPZeros(:,:,7)>0.5 & ANOVAAntVOD(:,:,7)<0.05)=1;
AntSMCount(TTPZeros(:,:,7)>0.5 & ANOVAAntSM(:,:,7)<0.05)=1;
SMPulseCount(TTPZeros(:,:,7)>0.5 & ANOVASMPulse(:,:,7)<0.05)=1;

PerRZSM = nansum(AntRZSMCount(:))/nansum(TTPFlag(:));
PerAntVOD = nansum(AntVODCount(:))/nansum(TTPFlag(:));
PerAntSM = nansum(AntSMCount(:))/nansum(TTPFlag(:));
PerSMPulse = nansum(SMPulseCount(:))/nansum(TTPFlag(:));

AntRZSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
AntVODCount = nan(size(ANOVAAntRZSM(:,:,7)));
AntSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
SMPulseCount = nan(size(ANOVAAntRZSM(:,:,7)));

AntRZSMCount(TTPZeros(:,:,7)>0.5 & ANOVAAntRZSM(:,:,7)<0.05 ...
    & -1*ANOVAAntRZSM(:,:,4)>0)=1;
AntVODCount(TTPZeros(:,:,7)>0.5 & ANOVAAntVOD(:,:,7)<0.05 ...
    & (-1*ANOVAAntVOD(:,:,4))>0)=1;
AntSMCount(TTPZeros(:,:,7)>0.5 & ANOVAAntSM(:,:,7)<0.05 ...
    & -1*ANOVAAntSM(:,:,4)>0)=1;
SMPulseCount(TTPZeros(:,:,7)>0.5 & ANOVASMPulse(:,:,7)<0.05 ...
    & -1*ANOVASMPulse(:,:,4)>0)=1;

PerRZSMPos = nansum(AntRZSMCount(:))/nansum(TTPFlag(:));
PerAntVODPos = nansum(AntVODCount(:))/nansum(TTPFlag(:));
PerAntSMPos = nansum(AntSMCount(:))/nansum(TTPFlag(:));
PerSMPulsePos = nansum(SMPulseCount(:))/nansum(TTPFlag(:));

clayAf = ClayMat;
sandAf = SandMat;

SatMatSand = -12.1;
SatMatClay = -40.5;
SatMatLoam = -47.8;
bSand = 4.05;
bClay = 11.4;
bLoam = 5.39;
nSand = 0.395;
nClay = 0.482;
nLoam = 0.451;
KSSand = 63.36;
KSClay = 0.46;
KSLoam = 2.5;
loamAf = 1-clayAf-sandAf;
beffMat = bSand.*sandAf+bClay.*clayAf+bLoam.*loamAf;
neffMat = nSand.*sandAf+nClay.*clayAf+nLoam.*loamAf;

AfFlag = ones(size(Bareflag));
AfFlag(50:95,446:500)=nan;
AfFlag(50:77,434:500)=nan;
AfFlag(50:72,429:500)=nan;
AfFlag(50:88,438:500)=nan;
AfFlag(50:91,443:500)=nan;
AfFlag(50:91,443:500)=nan;
AfFlag(1:53,384:500)=nan;
AfFlag(45:48,341:357)=nan;
AfFlag(:,500:end)=nan;
AfFlag(:,1:310)=nan;
AfFlag(1:47,:)=nan;

percZero = TTPZeros(:,:,7);
percZero(percZero<0.5)=nan;
percZero(isfinite(percZero))=1;
percZeroTotal = percZero;
percZeroAfOnly = percZeroTotal.*AfFlag;
nansum(percZeroAfOnly(:))./nansum(percZeroTotal(:))

ForestFlagMat = TreeCoverMat;
ForestFlagMat(ForestFlagMat>40)=nan;
ForestFlagMat(isfinite(ForestFlagMat))=1;

%% Median Time to Peak

TTPZeroMean = TTPZeros(:,:,1).*ForestFlagMat; %median
TTPZeroMeanNoZeros = TTPZeroMean;
TTPZeroMeanNoZeros(TTPZeroMeanNoZeros==0)=nan;
nanmedian(TTPZeroMeanNoZeros(:))

Ncr = 150;
Ncb = 400;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBuMod = flipud([blu;flipud(red(Ncr*(1/3):Ncr,:))]);


figure
iplotC = subplot(2,4,[1:4]);
pcolor(TopLeftCornerLon,TopLeftCornerLat,...
    TTPZeroMean) ; shading flat
hold on
geoshow(clat,clon,'LineWidth',1,'Color','k')
colormap(cmapRedBuMod)
set(gca,'color',[0.8 0.8 0.8])
ylim([-60.01 60.01])
caxis([0 5])
set(gca,'fontsize',18)
colorbar
set(gcf,'color','white')
set(gcf,'position',[100 100 1000 400])

labeltxt = text(208,27,'Median {\itt_p} (Days)','fontsize',20,'color','k');
set(labeltxt,'rotation',270)

iplotA = subplot(2,4,[5 6]);
A = TTPZeroMean;
Xticks = [1 2 3 4 5 6 7];
Xticks1 = {'0' '0-250' '250-500' '500-1000' '1000-2000' '2000-5000'};
edge = [0 250 500 1000 2000 5000];
[n,bin] = histc(PrecipMat(:),edge);
boxplot(A(:),bin,'notch','on');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'y','FaceAlpha',0.2);
end
hold on
set(findobj(gca,'type','line'),'linew',1.5)
h=findobj(gca,'tag','Outliers');
xlim([1.5 6.5])
delete(h) 
ylabel('Median {\itt_p} (Days)')
grid on
xlabel('Mean Annual Precipitation (mm)')
set(gca,'fontsize',15)
set(gca,'fontname','arial')
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
ylim([-0.5 5])
xtickangle(10)

iplotB = subplot(2,4,[7 8]);
A = TTPZeroMean;
Xticks = [1 2 3 4 5 6];
Xticks1 = {'0' '0-5' '5-10' '10-20' '20-30' '30-40'};
edge = [0 5 10 20 30 40];
[n,bin] = histc(TreeCoverMat(:),edge);
boxplot(A(:),bin,'notch','on');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'y','FaceAlpha',0.2);
end
hold on
% plot([0 8],[0.50 0.50],'--k','linewidth',1.5)
set(findobj(gca,'type','line'),'linew',1.5)
h=findobj(gca,'tag','Outliers');
xlim([1.5 6.5])
delete(h) 
grid on
xlabel('Tree Cover (%)')
set(gca,'fontsize',15)
set(gca,'fontname','arial')
ylim([-0.5 5])
ylabel('Median {\itt_p} (Days)')
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
% ylim([0 1])
% ylim([0.25 0.75])
xtickangle(10)

set(gcf,'position',[10 10 1000 600])

AMSA = get(iplotA,'Position'); % top right fig
set(iplotA,'Position',[AMSA(1)-0.01 AMSA(2)-0.02 AMSA(3)-0.02 AMSA(4)-0.05]);
AMSA = get(iplotA,'Position'); % top right fig

AMSB = get(iplotB,'Position'); % top right fig
set(iplotB,'Position',[AMSB(1) AMSA(2) AMSA(3) AMSA(4)]);
AMSB = get(iplotB,'Position'); % top right fig

AMSC = get(iplotC,'Position'); % top right fig
set(iplotC,'Position',[AMSC(1)-0.01 AMSC(2)-0.10 AMSC(3)+0.05 AMSC(4)+0.13]);
AMSC = get(iplotC,'Position'); % top right fig

Atext = ['B'];
annotation('textbox',[0.05 0.32 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.49 0.32 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['A'];
annotation('textbox',[0.05 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')

%% Time Series Examples
for iP = 1:4
 if iP == 1
     load('dVWC_Pulse_AfricanGrassland.mat')
    iplotD1 = subplot(2,2,iP);
 elseif iP ==2
load('dVWC_Pulse_CONUSShrubland.mat')
    iplotD2 = subplot(2,2,iP);
 elseif iP == 3
load('dVWC_Pulse_AfricanWSavanna.mat')
    iplotD3 = subplot(2,2,iP);
 elseif iP ==4
load('dVWC_Pulse_BrazilSavanna.mat')
    iplotD4 = subplot(2,2,iP);
 end

dVWCAllN = nan(1,20);
for i = 1:size(dVWCVecAll,2)
    dVWCA = dVWCVecAll(:,i);
    dVWCA(isnan(dVWCA))=[];
    dVWCAllN(i) = length(dVWCA);
end
dVWCAllN = dVWCAllN./size(dVWCVecAll,1);


tPulse = (0:size(dVWCVecAll,2)-1)-0.5;
dVWCVec = prctile(dVWCVecAll,[25 50 75],1);
missing = dVWCAllN<0.05;
dVWCVec(:,missing)=nan;
Notmissing = ~isnan(dVWCVec(1,:));

patch([-1 0 0 -1],[-0.1 -0.1 0.1 0.1],[0.5 0.5 0.5],...
    'facealpha',0.5,'linestyle','none')
hold on
binwidth = 0.125;
for ip1 = 1:length(tPulse)
    plot([tPulse(ip1)-binwidth tPulse(ip1)+binwidth],[dVWCVec(3,ip1) dVWCVec(3,ip1)],'-k')   
    hold on
    plot([tPulse(ip1)-binwidth tPulse(ip1)+binwidth],[dVWCVec(1,ip1) dVWCVec(1,ip1)],'-k')  
    hold on
    plot([tPulse(ip1)-binwidth tPulse(ip1)-binwidth],[dVWCVec(1,ip1) dVWCVec(3,ip1)],'-k')     
    hold on
    plot([tPulse(ip1)+binwidth tPulse(ip1)+binwidth],[dVWCVec(1,ip1) dVWCVec(3,ip1)],'-k') 
    hold on
    patch([tPulse(ip1)-binwidth tPulse(ip1)-binwidth tPulse(ip1)+binwidth tPulse(ip1)+binwidth],...
    [dVWCVec(1,ip1) dVWCVec(3,ip1) dVWCVec(3,ip1) dVWCVec(1,ip1)],...
    [0 0 0.8],'facealpha',0.1,'linestyle','none')
end
hold on
p2=plot(tPulse(Notmissing),dVWCVec([2],Notmissing),'.-k','linewidth',2,'markersize',30);
hold on
plot([-1 15],[0 0],'--k','linewidth',1)
set(gca,'xtick',[-1:1:8])
xlim([tPulse(1)-0.5 8])

if iP ==1|| iP==3
ylabel([{'Normalized'};{'dVOD/dt (1/Day)'}])
end
grid on
set(gca,'fontsize',14)
hold on
xlabel('Time After Pulse (Days)')

 if iP == 1
title('Sahelian Grasslands')
 ylim([-0.1 0.1])
 text(3.5,0.07,[{'Water'}; {'Storage'}],'fontsize',...
     14,'color','k','horizontalalignment','center')
text(3.5,-0.07,[{'Water'}; {'Loss'}],'fontsize',...
    14,'color','k','horizontalalignment','center')
annotation('textarrow',[0.3 0.3],[0.79 0.85])
annotation('textarrow',[0.3 0.3],[0.74 0.68])
 elseif iP ==2
title('Southwest US Shrublands')
ylim([-0.1 0.1])
 elseif iP == 3
title('African Woody Savannas') 
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.025:0.05])
 elseif iP ==4
title('Brazilian Savannas')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.025:0.05])
end
 
if iP == 1
axes('position',[0.343 0.835 0.12 0.11])
elseif iP == 2
axes('position',[0.785 0.835 0.12 0.11])    
elseif iP == 3
axes('position',[0.343 0.36 0.12 0.11])
elseif iP == 4  
axes('position',[0.785 0.36 0.12 0.11])        
end

pcolor(TopLeftCornerLon,TopLeftCornerLat,...
    TTPZeroMean) ; shading flat
hold on
geoshow(clat,clon,'LineWidth',1,'Color','k')
colormap(cmapRedBuMod)
set(gca,'color',[0.8 0.8 0.8])
caxis([0 5])
set(gca,'fontsize',18)
set(gca,'xticklabel','')
set(gca,'yticklabel','')
set(gcf,'color','white')
set(gcf,'position',[100 100 1000 400])


if iP == 1
ylim([-10 40]) 
xlim([-30 60])
 hold on
rectangle('Position',[-10  10 46 7],'linewidth',2)
elseif iP == 2
ylim([10 60]) 
xlim([-145 -55])
hold on
rectangle('Position',[-120 25 20 15],'linewidth',2)
elseif iP == 3
ylim([-35 15]) 
xlim([-30 60]) 
hold on
rectangle('Position',[ 9 -15  21 14 ],'linewidth',2)
elseif iP == 4 
ylim([-35 15]) 
xlim([-110 -20]) 
hold on
rectangle('Position',[-60 -28 20  13],'linewidth',2)
end

end


AMSD1 = get(iplotD1,'Position'); % top right fig
set(iplotD1,'Position',[AMSD1(1) AMSD1(2) AMSD1(3) AMSD1(4)+0.02]);
AMSD1 = get(iplotD1,'Position'); % top right fig

AMSD2 = get(iplotD2,'Position'); % top right fig
set(iplotD2,'Position',[AMSD2(1) AMSD2(2) AMSD1(3) AMSD1(4)]);
AMSD2 = get(iplotD2,'Position'); % top right fig

AMSD3 = get(iplotD3,'Position'); % top right fig
set(iplotD3,'Position',[AMSD3(1) AMSD3(2) AMSD1(3) AMSD1(4)]);
AMSD3 = get(iplotD3,'Position'); % top right fig

AMSD4 = get(iplotD4,'Position'); % top right fig
set(iplotD4,'Position',[AMSD4(1) AMSD4(2) AMSD1(3) AMSD1(4)]);
AMSD4 = get(iplotD4,'Position'); % top right fig

Atext = ['A'];
annotation('textbox',[0.05 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.5 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.05 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['D'];
annotation('textbox',[0.5 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')

set(gcf,'position',[10 10 800 600])

%% Pulse Characteristics

rTTPCorrSMAnt = TTPCorrSMAnt(:,:,1);
rTTPCorrSMAnt0 = TTPCorrSMAnt0(:,:,1);
rTTPCorrSMPulse = TTPCorrSMPulse(:,:,1);
rTTPCorrSMPulse0 = TTPCorrSMPulse0(:,:,1);
rTTPCorrSMRootAnt = TTPCorrSMRootAnt(:,:,1);
rTTPCorrSMRootAnt0 = TTPCorrSMRootAnt0(:,:,1);
rTTPCorrVODAnt = TTPCorrVODAnt(:,:,1);
rTTPCorrVODAnt0 = TTPCorrVODAnt0(:,:,1);

percZero = TTPZeros(:,:,7);
percZeroFlag = percZero;
percZeroFlag(percZero<0.5)=nan;
percZeroFlag(percZeroFlag>0)=1;

rTTPCorrSMAnt(percZero<0.5) = nan;
rTTPCorrSMAnt0(percZero<0.5) = nan;
rTTPCorrSMPulse(percZero<0.5) = nan;
rTTPCorrSMPulse0(percZero<0.5) = nan;
rTTPCorrSMRootAnt(percZero<0.5) = nan;
rTTPCorrSMRootAnt0(percZero<0.5) = nan;
rTTPCorrVODAnt(percZero<0.5) = nan;
rTTPCorrVODAnt0(percZero<0.5) = nan;

xbox3 = [1:1:4];

rTTPCorrSMPulsePosSig = rTTPCorrSMPulse;
rTTPCorrSMPulsePosSig(rTTPCorrSMPulsePosSig<0 | TTPCorrSMPulse(:,:,2)>0.05)=nan;
rTTPCorrSMPulsePosSig(isfinite(rTTPCorrSMPulsePosSig))=1;
rTTPCorrSMPulseNegSig = rTTPCorrSMPulse;
rTTPCorrSMPulseNegSig(rTTPCorrSMPulseNegSig>0 | TTPCorrSMPulse(:,:,2)>0.05)=nan;
rTTPCorrSMPulseNegSig(isfinite(rTTPCorrSMPulseNegSig))=1;
rTTPCorrSMPulseTot = rTTPCorrSMPulse;
rTTPCorrSMPulseTot(isfinite(rTTPCorrSMPulseTot))=1;
SMPulsePosSig = nansum(rTTPCorrSMPulsePosSig(:))/nansum(rTTPCorrSMPulseTot(:));
SMPulseNegSig = nansum(rTTPCorrSMPulseNegSig(:))/nansum(rTTPCorrSMPulseTot(:));

rTTPCorrSMAntPosSig = rTTPCorrSMAnt;
rTTPCorrSMAntPosSig(rTTPCorrSMAntPosSig<0 | TTPCorrSMAnt(:,:,2)>0.05)=nan;
rTTPCorrSMAntPosSig(isfinite(rTTPCorrSMAntPosSig))=1;
rTTPCorrSMAntNegSig = rTTPCorrSMAnt;
rTTPCorrSMAntNegSig(rTTPCorrSMAntNegSig>0 | TTPCorrSMAnt(:,:,2)>0.05)=nan;
rTTPCorrSMAntNegSig(isfinite(rTTPCorrSMAntNegSig))=1;
rTTPCorrSMAntTot = rTTPCorrSMAnt;
rTTPCorrSMAntTot(isfinite(rTTPCorrSMAntTot))=1;
SMAntPosSig = nansum(rTTPCorrSMAntPosSig(:))/nansum(rTTPCorrSMAntTot(:));
SMAntNegSig = nansum(rTTPCorrSMAntNegSig(:))/nansum(rTTPCorrSMAntTot(:));

rTTPCorrSMRootAntPosSig = rTTPCorrSMRootAnt;
rTTPCorrSMRootAntPosSig(rTTPCorrSMRootAntPosSig<0 | TTPCorrSMRootAnt(:,:,2)>0.05)=nan;
rTTPCorrSMRootAntPosSig(isfinite(rTTPCorrSMRootAntPosSig))=1;
rTTPCorrSMRootAntNegSig = rTTPCorrSMRootAnt;
rTTPCorrSMRootAntNegSig(rTTPCorrSMRootAntNegSig>0 | TTPCorrSMRootAnt(:,:,2)>0.05)=nan;
rTTPCorrSMRootAntNegSig(isfinite(rTTPCorrSMRootAntNegSig))=1;
rTTPCorrSMRootAntTot = rTTPCorrSMRootAnt;
rTTPCorrSMRootAntTot(isfinite(rTTPCorrSMRootAntTot))=1;
SMRootAntPosSig = nansum(rTTPCorrSMRootAntPosSig(:))/nansum(rTTPCorrSMRootAntTot(:));
SMRootAntNegSig = nansum(rTTPCorrSMRootAntNegSig(:))/nansum(rTTPCorrSMRootAntTot(:));

rTTPCorrVODAntPosSig = rTTPCorrVODAnt;
rTTPCorrVODAntPosSig(rTTPCorrVODAntPosSig<0 | TTPCorrVODAnt(:,:,2)>0.05)=nan;
rTTPCorrVODAntPosSig(isfinite(rTTPCorrVODAntPosSig))=1;
rTTPCorrVODAntNegSig = rTTPCorrVODAnt;
rTTPCorrVODAntNegSig(rTTPCorrVODAntNegSig>0 | TTPCorrVODAnt(:,:,2)>0.05)=nan;
rTTPCorrVODAntNegSig(isfinite(rTTPCorrVODAntNegSig))=1;
rTTPCorrVODAntTot = rTTPCorrVODAnt;
rTTPCorrVODAntTot(isfinite(rTTPCorrVODAntTot))=1;
VODAntPosSig = nansum(rTTPCorrVODAntPosSig(:))/nansum(rTTPCorrVODAntTot(:));
VODAntNegSig = nansum(rTTPCorrVODAntNegSig(:))/nansum(rTTPCorrVODAntTot(:));

percZeroFlag = percZero;
percZeroFlag(percZeroFlag<0.5)=nan;
percZeroFlag(isfinite(percZeroFlag))=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Antecedent Soil Moisture %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PulseTTPZero = ANOVAAntSM(:,:,8);
PulseTTPShort = ANOVAAntSM(:,:,9);
PulseTTPLong = ANOVAAntSM(:,:,10);
pAnovaAntSM = nan(4,1);
AntSMIGBP = nan(3,1);

PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
PulseGroups = [PulseTTPZeroit(:)' PulseTTPShortit(:)' PulseTTPLongit(:)'];
groupAnova = [1*ones(1,length(PulseTTPZeroit(:)))...
    2*ones(1,length(PulseTTPShortit(:))) ...
    3*ones(1,length(PulseTTPLongit(:)))];
groupAnova(isnan(PulseGroups))=[];
PulseGroups(isnan(PulseGroups))=[];
[A tbl stats1] = anovan(PulseGroups,groupAnova','display','off');
cAntSM = multcompare(stats1,'displayopt','off');    
pAnovaAntSM(1) = A;
pAnovaAntSM(2:4) = cAntSM(:,end); 

pAnovaPulseVec = pAnovaAntSM;
pAnovaPulseVec(pAnovaPulseVec<0.05)=0;
pAnovaPulseVec(pAnovaPulseVec~=0)=1;
pAnovaPulseVec = nansum(pAnovaPulseVec,1);
pAnovaPulseVec23 = pAnovaAntSM(4,:);
pAnovaPulseVec23(pAnovaPulseVec23<0.05)=1;
pAnovaPulseVec(pAnovaPulseVec==1 & pAnovaPulseVec23==1)=2;

AntSMIGBP(1) = nanmean(PulseTTPZeroit(:));
AntSMIGBP(2) = nanmean(PulseTTPShortit(:));
AntSMIGBP(3) = nanmean(PulseTTPLongit(:)); 

[h0KWPulse p0KWAntSM] = jbtest(PulseTTPZeroit(:));
[h1KWPulse p1KWAntSM] = jbtest(PulseTTPShortit(:));
[h2KWPulse p2KWAntSM] = jbtest(PulseTTPLongit(:));
[pKWPulse,tbl,stats] = kruskalwallis(PulseGroups',groupAnova','off');
chi_AntSM = tbl{2,5};
[p12] = ranksum(PulseTTPZeroit(:)',PulseTTPShortit(:)')
[p23] = ranksum(PulseTTPShortit(:)',PulseTTPLongit(:)')
[p13] = ranksum(PulseTTPZeroit(:)',PulseTTPLongit(:)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Antecedent Rootzone Soil Moisture %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseTTPZero = ANOVAAntRZSM(:,:,8).*neffMat;
PulseTTPShort = ANOVAAntRZSM(:,:,9).*neffMat;
PulseTTPLong = ANOVAAntRZSM(:,:,10).*neffMat;
pAnovaAntRZSM = nan(4,1);
AntRZSMIGBP = nan(3,1);

PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
PulseGroups = [PulseTTPZeroit(:)' PulseTTPShortit(:)' PulseTTPLongit(:)'];
groupAnova = [1*ones(1,length(PulseTTPZeroit(:)))...
    2*ones(1,length(PulseTTPShortit(:))) ...
    3*ones(1,length(PulseTTPLongit(:)))];
groupAnova(isnan(PulseGroups))=[];
PulseGroups(isnan(PulseGroups))=[];
[A tbl stats1] = anovan(PulseGroups,groupAnova','display','off');
cAntRZSM = multcompare(stats1,'displayopt','off');    
pAnovaAntRZSM(1) = A;
pAnovaAntRZSM(2:4) = cAntRZSM(:,end); 

pAnovaPulseVec = pAnovaAntRZSM;
pAnovaPulseVec(pAnovaPulseVec<0.05)=0;
pAnovaPulseVec(pAnovaPulseVec~=0)=1;
pAnovaPulseVec = nansum(pAnovaPulseVec,1);
pAnovaPulseVec23 = pAnovaAntRZSM(4,:);
pAnovaPulseVec23(pAnovaPulseVec23<0.05)=1;
pAnovaPulseVec(pAnovaPulseVec==1 & pAnovaPulseVec23==1)=2;

AntRZSMIGBP(1) = nanmean(PulseTTPZeroit(:));
AntRZSMIGBP(2) = nanmean(PulseTTPShortit(:));
AntRZSMIGBP(3) = nanmean(PulseTTPLongit(:)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Pulse Soil Moisture Magnitude %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseTTPZero = ANOVASMPulse(:,:,8);
PulseTTPShort = ANOVASMPulse(:,:,9);
PulseTTPLong = ANOVASMPulse(:,:,10);
pAnovaSMPulse = nan(4,1);
SMPulseIGBP = nan(3,1);

PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
PulseGroups = [PulseTTPZeroit(:)' PulseTTPShortit(:)' PulseTTPLongit(:)'];
groupAnova = [1*ones(1,length(PulseTTPZeroit(:)))...
    2*ones(1,length(PulseTTPShortit(:))) ...
    3*ones(1,length(PulseTTPLongit(:)))];
groupAnova(isnan(PulseGroups))=[];
PulseGroups(isnan(PulseGroups))=[];
[A tbl stats1] = anovan(PulseGroups,groupAnova','display','off');
cPulseSize = multcompare(stats1,'displayopt','off');    
pAnovaSMPulse(1) = A;
pAnovaSMPulse(2:4) = cPulseSize(:,end); 

pAnovaPulseVec = pAnovaSMPulse;
pAnovaPulseVec(pAnovaPulseVec<0.05)=0;
pAnovaPulseVec(pAnovaPulseVec~=0)=1;
pAnovaPulseVec = nansum(pAnovaPulseVec,1);
pAnovaPulseVec23 = pAnovaSMPulse(4,:);
pAnovaPulseVec23(pAnovaPulseVec23<0.05)=1;
pAnovaPulseVec(pAnovaPulseVec==1 & pAnovaPulseVec23==1)=2;

SMPulseIGBP(1) = nanmean(PulseTTPZeroit(:));
SMPulseIGBP(2) = nanmean(PulseTTPShortit(:));
SMPulseIGBP(3) = nanmean(PulseTTPLongit(:)); 

[h0KWPulse p0KWPulse] = jbtest(PulseTTPZeroit(:));
[h1KWPulse p1KWPulse] = jbtest(PulseTTPShortit(:));
[h2KWPulse p2KWPulse] = jbtest(PulseTTPLongit(:));
[pKWPulse,tbl,stats] = kruskalwallis(PulseGroups',groupAnova','off');
chi_PulseSize = tbl{2,5};
[p12] = ranksum(PulseTTPZeroit(:)',PulseTTPShortit(:)')
[p23] = ranksum(PulseTTPShortit(:)',PulseTTPLongit(:)')
[p13] = ranksum(PulseTTPZeroit(:)',PulseTTPLongit(:)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Antecedent VOD %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseTTPZero = ANOVAAntVOD(:,:,8);
PulseTTPShort = ANOVAAntVOD(:,:,9);
PulseTTPLong = ANOVAAntVOD(:,:,10);
pAnovaAntVOD = nan(4,1);
AntVODIGBP = nan(3,1);

PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
PulseGroups = [PulseTTPZeroit(:)' PulseTTPShortit(:)' PulseTTPLongit(:)'];
groupAnova = [1*ones(1,length(PulseTTPZeroit(:)))...
    2*ones(1,length(PulseTTPShortit(:))) ...
    3*ones(1,length(PulseTTPLongit(:)))];
groupAnova(isnan(PulseGroups))=[];
PulseGroups(isnan(PulseGroups))=[];
[A tbl stats1] = anovan(PulseGroups,groupAnova','display','off');
cAntVOD = multcompare(stats1,'displayopt','off');    
pAnovaAntVOD(1) = A;
pAnovaAntVOD(2:4) = cAntVOD(:,end); 

pAnovaPulseVec = pAnovaAntVOD;
pAnovaPulseVec(pAnovaPulseVec<0.05)=0;
pAnovaPulseVec(pAnovaPulseVec~=0)=1;
pAnovaPulseVec = nansum(pAnovaPulseVec,1);
pAnovaPulseVec23 = pAnovaAntVOD(4,:);
pAnovaPulseVec23(pAnovaPulseVec23<0.05)=1;
pAnovaPulseVec(pAnovaPulseVec==1 & pAnovaPulseVec23==1)=2;

AntVODIGBP(1) = nanmean(PulseTTPZeroit(:));
AntVODIGBP(2) = nanmean(PulseTTPShortit(:));
AntVODIGBP(3) = nanmean(PulseTTPLongit(:)); 

[h0KWPulse p0KWAntVOD] = jbtest(PulseTTPZeroit(:));
[h1KWPulse p1KWAntVOD] = jbtest(PulseTTPShortit(:));
[h2KWPulse p2KWAntVOD] = jbtest(PulseTTPLongit(:));
[pKWPulse,tbl,stats] = kruskalwallis(PulseGroups',groupAnova','off');
chi_AntVOD = tbl{2,5};
[p12] = ranksum(PulseTTPZeroit(:)',PulseTTPShortit(:)')
[p23] = ranksum(PulseTTPShortit(:)',PulseTTPLongit(:)')
[p13] = ranksum(PulseTTPZeroit(:)',PulseTTPLongit(:)')


figure

PulseTTPZero = ANOVASMPulse(:,:,8);
PulseTTPShort = ANOVASMPulse(:,:,9);
PulseTTPLong = ANOVASMPulse(:,:,10);
PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag;

iplotC = subplot(2,2,2);
abox = boxplot([PulseTTPZeroit(:) PulseTTPShortit(:) PulseTTPLongit(:)],...
'positions',[1:3],'Colors', [0 0 0; 1 0 0; 0 0 0.8],'Notch','on');

h = findobj(gca,'Tag','Box');
colors = flipud([0 0 0; 1 0 0; 0 0 1]);
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.1);
end
set(gca,'xticklabel','')
set(findobj(gca,'type','line'),'linew',1.5)
ylabel([{'Mean Soil Moisture'}; {'Pulse Magnitude [m^3 m^{-3}]'}])
grid on
set(gca,'fontsize',16)
h=findobj(gca,'tag','Outliers');
delete(h) 
ylim([0 0.2])
text(1,-0.023,[{'{\itt_p}=0'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.023,[{'1?{\itt_p}?3'}],'fontsize',16,'horizontalalignment','center')
text(3,-0.023,[{'{\itt_p}>3'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.05,[{'Time to Peak Plant Water Content'}],'fontsize',16,'horizontalalignment','center')

PulseTTPZero = ANOVAAntSM(:,:,8);
PulseTTPShort = ANOVAAntSM(:,:,9);
PulseTTPLong = ANOVAAntSM(:,:,10);
PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 

iplotA = subplot(2,2,1);
boxplot([PulseTTPZeroit(:) PulseTTPShortit(:) PulseTTPLongit(:)],...
'positions',[1:3],'Colors', [0 0 0; 1 0 0; 0 0 0.8],'Notch','on')

h = findobj(gca,'Tag','Box');
colors = flipud([0 0 0; 1 0 0; 0 0 1]);
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.1);
end
set(gca,'xticklabel','')
grid on
ylabel([{'Mean Antecedent Soil'}; {'Moisture [m^3 m^{-3}]'}])
set(gca,'fontsize',16)
h=findobj(gca,'tag','Outliers');
text(1,-0.023,[{'{\itt_p}=0'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.023,[{'1?{\itt_p}?3'}],'fontsize',16,'horizontalalignment','center')
text(3,-0.023,[{'{\itt_p}>3'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.06,[{'Time to Peak Plant Water Content'}],'fontsize',16,'horizontalalignment','center')
delete(h) 
ylim([0 0.3])
set(findobj(gca,'type','line'),'linew',1.5)

PulseTTPZero = ANOVAAntRZSM(:,:,8).*neffMat;
PulseTTPShort = ANOVAAntRZSM(:,:,9).*neffMat;
PulseTTPLong = ANOVAAntRZSM(:,:,10).*neffMat;
PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
diffPulse = PulseTTPLongit-PulseTTPShortit;

PulseTTPZero = ANOVAAntVOD(:,:,8);
PulseTTPShort = ANOVAAntVOD(:,:,9);
PulseTTPLong = ANOVAAntVOD(:,:,10);
PulseTTPZeroit = PulseTTPZero.*percZeroFlag;    
PulseTTPShortit = PulseTTPShort.*percZeroFlag;    
PulseTTPLongit = PulseTTPLong.*percZeroFlag; 
diffPulse = PulseTTPLongit-PulseTTPShortit;

iplotD = subplot(2,2,3);
boxplot([PulseTTPZeroit(:) PulseTTPShortit(:) PulseTTPLongit(:)],...
'positions',[1:3],'Colors', [0 0 0; 1 0 0; 0 0 0.8],'Notch','on')

h = findobj(gca,'Tag','Box');
colors = flipud([0 0 0; 1 0 0; 0 0 1]);
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.1);
end

set(gca,'xticklabel','')
grid on
ylabel([{'Mean Antecedent VOD'}])
set(gca,'fontsize',16)
ylim([0 0.6])
set(findobj(gca,'type','line'),'linew',1.5)
h=findobj(gca,'tag','Outliers');
delete(h) 
text(1,-0.023*2,[{'{\itt_p}=0'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.023*2,[{'1?{\itt_p}?3'}],'fontsize',16,'horizontalalignment','center')
text(3,-0.023*2,[{'{\itt_p}>3'}],'fontsize',16,'horizontalalignment','center')
text(2,-0.05*2.5,[{'Time to Peak Plant Water Content'}],'fontsize',16,'horizontalalignment','center')

AMSA = get(iplotA,'Position'); % top right fig
set(iplotA,'Position',[AMSA(1)-0.08 AMSA(2)+0.03 AMSA(3)+0.08 AMSA(4)]);
AMSA = get(iplotA,'Position'); % top right fig

AMSC = get(iplotC,'Position'); % top right fig
set(iplotC,'Position',[AMSC(1)-0.06 AMSA(2) AMSA(3) AMSA(4)]);
AMSC = get(iplotC,'Position'); % top right fig

AMSD = get(iplotD,'Position'); % top right fig
set(iplotD,'Position',[AMSD(1)+0.17 AMSD(2)-0.01 AMSA(3) AMSA(4)]);
AMSD = get(iplotD,'Position'); % top right fig

set(gcf,'position',[100 100 800 600])

Atext = ['A'];
annotation('textbox',[0.02 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.49 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.26 0.4 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
