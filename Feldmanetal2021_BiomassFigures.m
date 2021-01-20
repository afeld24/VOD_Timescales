clc
clear

load coast   
clat = lat;                                                        
clon = long ; 

cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy_TauEfolding/Figures_ToDistribute')
dr = cat(1,dir('*TimeToPeak_Africa_ANOVA_V25_N*.mat'),...
    dir('*TimeToPeak_Africa_ANOVA_V25_S*.mat'));
load('GlobeAncillaryHalfDegree')

Ncr = 200;
Ncb = 200;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [blu;flipud(red)];

Regime4Flag_Africa_Full = [];
R2Mat_Full = [];
STDMat_Full = [];
betaMat_Full = [];
upperBoundMat_Full = [];
ZerosDiffMat_Full = [];
TTPZeros_Full = [];
DrydownCharMat_Full = [];
ANOVAAntRZSM_Full = [];
ANOVAAntSM_Full = [];
ANOVAAntVOD_Full = [];
ANOVALAI3_Full = [];
ANOVAFAPAR3_Full = [];
ANOVALAI3Median_Full = [];
ANOVALAI5_Full = [];
ANOVASMPulse_Full = [];
RZInfTimeMean_Full = [];
RZInfTimeMedian_Full = [];
RZInfTimePer_Full = [];

TTPBins_Full = []; 
TTPGrowth_Full = [];
AntSMGrowth_Full = [];
pulsesizeGrowth_Full = [];

ANOVALAICorr_Full = [];
ANOVALAICorr0_Full = [];

SMPulseSeasMean_Full = [];
SMPulseSeasNShort_Full = [];
SMPulseSeasNLong_Full = [];
SMPulseSeasNRapid_Full = [];
SMPulseSeasMeanLAI_Full = [];
SMPulseSeasNShortLAI_Full = [];
SMPulseSeasNLongLAI_Full = [];
SMPulseSeasNRapidLAI_Full = [];
SMPulseSeasMeanLAIMODIS_Full = [];
SMPulseSeasNShortLAIMODIS_Full = [];
SMPulseSeasNLongLAIMODIS_Full = [];
SMPulseSeasNRapidLAIMODIS_Full = [];

ANOVALAI3_MODIS_Full = [];

for j = 1:length(dr)
load(dr(j).name)    

TTPBins_Full = ...
    cat(4,TTPBins_Full,TTPBins);
TTPGrowth_Full = ...
    cat(4,TTPGrowth_Full,TTPGrowth);
AntSMGrowth_Full = ...
    cat(4,AntSMGrowth_Full,AntSMGrowth);
pulsesizeGrowth_Full = ...
    cat(4,pulsesizeGrowth_Full,pulsesizeGrowth);

Regime4Flag_Africa_Full = ...
    cat(4,Regime4Flag_Africa_Full,Regime4Flag_Africa);
R2Mat_Full = ...
    cat(4,R2Mat_Full,R2Mat);
STDMat_Full = ...
    cat(4,STDMat_Full,STDMat);
betaMat_Full = ...
    cat(4,betaMat_Full,betaMat);
TTPZeros_Full = ...
    cat(4,TTPZeros_Full,TTPZeros);
upperBoundMat_Full = ...
    cat(4,upperBoundMat_Full,upperBoundMat);
DrydownCharMat_Full = ...
    cat(4,DrydownCharMat_Full,DrydownCharMat);

ANOVAFAPAR3_Full = ...
    cat(4,ANOVAFAPAR3_Full,ANOVAFAPAR3);

ANOVAAntRZSM_Full = ...
    cat(4,ANOVAAntRZSM_Full,ANOVAAntRZSM);
ANOVALAI3_Full = ...
    cat(4,ANOVALAI3_Full,ANOVALAI3);
ANOVALAI3Median_Full = ...
    cat(4,ANOVALAI3Median_Full,ANOVALAI3Median);
ANOVALAI5_Full = ...
    cat(4,ANOVALAI5_Full,ANOVALAI5);
ANOVAAntVOD_Full = ...
    cat(4,ANOVAAntVOD_Full,ANOVAAntVOD);
ANOVAAntSM_Full = ...
    cat(4,ANOVAAntSM_Full,ANOVAAntSM);
ANOVASMPulse_Full = ...
    cat(4,ANOVASMPulse_Full,ANOVASMPulse);

ANOVALAI3_MODIS_Full = ...
    cat(4,ANOVALAI3_MODIS_Full,ANOVALAI3_MODIS);

ANOVALAICorr_Full = ...
    cat(4,ANOVALAICorr_Full,ANOVALAICorr);
ANOVALAICorr0_Full = ...
    cat(4,ANOVALAICorr0_Full,ANOVALAICorr0);

RZInfTimeMean_Full = ...
    cat(4,RZInfTimeMean_Full,RZInfTimeMean);
RZInfTimeMedian_Full = ...
    cat(4,RZInfTimeMedian_Full,RZInfTimeMedian);
RZInfTimePer_Full = ...
    cat(4,RZInfTimePer_Full,RZInfTimePer);

SMPulseSeasMean_Full = ...
    cat(4,SMPulseSeasMean_Full,SMPulseSeasMean);
SMPulseSeasNShort_Full = ...
    cat(4,SMPulseSeasNShort_Full,SMPulseSeasNShort);
SMPulseSeasNLong_Full = ...
    cat(4,SMPulseSeasNLong_Full,SMPulseSeasNLong);
SMPulseSeasNRapid_Full = ...
    cat(4,SMPulseSeasNRapid_Full,SMPulseSeasNRapid);

SMPulseSeasMeanLAI_Full = ...
    cat(4,SMPulseSeasMeanLAI_Full,SMPulseSeasMeanLAI);
SMPulseSeasNShortLAI_Full = ...
    cat(4,SMPulseSeasNShortLAI_Full,SMPulseSeasNShortLAI);
SMPulseSeasNLongLAI_Full = ...
    cat(4,SMPulseSeasNLongLAI_Full,SMPulseSeasNLongLAI);
SMPulseSeasNRapidLAI_Full = ...
    cat(4,SMPulseSeasNRapidLAI_Full,SMPulseSeasNRapidLAI);

SMPulseSeasMeanLAIMODIS_Full = ...
    cat(4,SMPulseSeasMeanLAIMODIS_Full,SMPulseSeasMeanLAIMODIS);
SMPulseSeasNShortLAIMODIS_Full = ...
    cat(4,SMPulseSeasNShortLAIMODIS_Full,SMPulseSeasNShortLAIMODIS);
SMPulseSeasNLongLAIMODIS_Full = ...
    cat(4,SMPulseSeasNLongLAIMODIS_Full,SMPulseSeasNLongLAIMODIS);
SMPulseSeasNRapidLAIMODIS_Full = ...
    cat(4,SMPulseSeasNRapidLAIMODIS_Full,SMPulseSeasNRapidLAIMODIS);
end

Regime4Flag_Africa = nanmean(Regime4Flag_Africa_Full,4);
R2Mat = nanmean(R2Mat_Full,4);
STDMat = nanmean(STDMat_Full,4);
betaMat = nanmean(betaMat_Full,4);
TTPZeros = nanmean(TTPZeros_Full,4);
upperBoundMat = nanmean(upperBoundMat_Full,4);
ZerosDiffMat = nanmean(ZerosDiffMat_Full,4);
DrydownCharMat = nanmean(DrydownCharMat_Full,4);

TTPBins = nanmean(TTPBins_Full,4);
TTPGrowth = nanmean(TTPGrowth_Full,4);
AntSMGrowth = nanmean(AntSMGrowth_Full,4);
pulsesizeGrowth = nanmean(pulsesizeGrowth_Full,4);

RZInfTimeMean = nanmean(RZInfTimeMean_Full,4);
RZInfTimeMedian = nanmean(RZInfTimeMedian_Full,4);
RZInfTimePer = nanmean(RZInfTimePer_Full,4);

SMPulseSeasMean = nanmean(SMPulseSeasMean_Full,4);
SMPulseSeasNShort = nanmean(SMPulseSeasNShort_Full,4);
SMPulseSeasNLong = nanmean(SMPulseSeasNLong_Full,4);
SMPulseSeasNRapid = nanmean(SMPulseSeasNRapid_Full,4);

SMPulseSeasMeanLAI = nanmean(SMPulseSeasMeanLAI_Full,4);
SMPulseSeasNShortLAI = nanmean(SMPulseSeasNShortLAI_Full,4);
SMPulseSeasNLongLAI = nanmean(SMPulseSeasNLongLAI_Full,4);
SMPulseSeasNRapidLAI = nanmean(SMPulseSeasNRapidLAI_Full,4);

SMPulseSeasMeanLAIMODIS = nanmean(SMPulseSeasMeanLAIMODIS_Full,4);
SMPulseSeasNShortLAIMODIS = nanmean(SMPulseSeasNShortLAIMODIS_Full,4);
SMPulseSeasNLongLAIMODIS = nanmean(SMPulseSeasNLongLAIMODIS_Full,4);
SMPulseSeasNRapidLAIMODIS = nanmean(SMPulseSeasNRapidLAIMODIS_Full,4);

ANOVAAntRZSM = nanmean(ANOVAAntRZSM_Full,4);
ANOVALAI = nanmean(ANOVALAI3_Full,4);
ANOVAFAPAR = nanmean(ANOVAFAPAR3_Full,4);
ANOVALAIMedian = nanmean(ANOVALAI3Median_Full,4);
ANOVALAI5 = nanmean(ANOVALAI5_Full,4);
ANOVAAntVOD = nanmean(ANOVAAntVOD_Full,4);
ANOVAAntSM = nanmean(ANOVAAntSM_Full,4);
ANOVASMPulse = nanmean(ANOVASMPulse_Full,4);

ANOVALAI3_MODIS = nanmean(ANOVALAI3_MODIS_Full,4);
ANOVALAICorr = nanmean(ANOVALAICorr_Full,4);
ANOVALAI0Corr = nanmean(ANOVALAICorr0_Full,4);

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

Bareflag(50:95,446:500)=nan;
Bareflag(50:77,434:500)=nan;
Bareflag(50:72,429:500)=nan;
Bareflag(50:88,438:500)=nan;
Bareflag(50:91,443:500)=nan;
Bareflag(50:91,443:500)=nan;
Bareflag(1:53,384:500)=nan;
Bareflag(45:48,341:357)=nan;


Regime4Flag_Africa = Regime4Flag_Africa.*Bareflag;
R2Mat = R2Mat.*Bareflag;
STDMat = STDMat.*Bareflag;
betaMat = betaMat.*Bareflag;
upperBoundMat = upperBoundMat.*Bareflag;
DrydownCharMat = DrydownCharMat.*Bareflag;
TTPZeros = TTPZeros.*Bareflag;

betaMat1 = betaMat.*STDMat;

ANOVAAntRZSM = ANOVAAntRZSM.*Bareflag;
ANOVALAI = ANOVALAI.*Bareflag;
ANOVAFAPAR = ANOVAFAPAR.*Bareflag;
ANOVALAIMedian = ANOVALAIMedian.*Bareflag;
ANOVALAI5 = ANOVALAI5.*Bareflag;
ANOVAAntVOD = ANOVAAntVOD.*Bareflag;
ANOVAAntSM = ANOVAAntSM.*Bareflag;
ANOVASMPulse = ANOVASMPulse.*Bareflag;

ANOVALAI3_MODIS = ANOVALAI3_MODIS.*Bareflag;

ANOVALAICorr = ANOVALAICorr.*Bareflag;
ANOVALAICorr0 = ANOVALAICorr0.*Bareflag;
RZInfTimeMean = RZInfTimeMean.*Bareflag;
RZInfTimeMedian = RZInfTimeMedian.*Bareflag;
RZInfTimePer = RZInfTimePer.*Bareflag;

SMPulseSeasMean = SMPulseSeasMean.*Bareflag;
SMPulseSeasNShort = SMPulseSeasNShort.*Bareflag;
SMPulseSeasNLong = SMPulseSeasNLong.*Bareflag;
SMPulseSeasNRapid = SMPulseSeasNRapid.*Bareflag;

SMPulseSeasMeanLAI = SMPulseSeasMeanLAI.*Bareflag;
SMPulseSeasNShortLAI = SMPulseSeasNShortLAI.*Bareflag;
SMPulseSeasNLongLAI = SMPulseSeasNLongLAI.*Bareflag;
SMPulseSeasNRapidLAI = SMPulseSeasNRapidLAI.*Bareflag;

SMPulseSeasMeanLAIMODIS = SMPulseSeasMeanLAIMODIS.*Bareflag;
SMPulseSeasNShortLAIMODIS = SMPulseSeasNShortLAIMODIS.*Bareflag;
SMPulseSeasNLongLAIMODIS = SMPulseSeasNLongLAIMODIS.*Bareflag;
SMPulseSeasNRapidLAIMODIS = SMPulseSeasNRapidLAIMODIS.*Bareflag;

TTPFlag = TTPZeros(:,:,5);
TTPFlag(TTPFlag>0.5)=1;
TTPFlag(TTPFlag<=0.5)=nan;

GrowthFlag = ANOVALAI(:,:,9);
GrowthFlag(GrowthFlag<0)=nan;
GrowthFlag(isfinite(GrowthFlag))=1;

% Calculate percentages for bins 2 and 3
AntRZSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
LAICount = nan(size(ANOVAAntRZSM(:,:,7)));
AntVODCount = nan(size(ANOVAAntRZSM(:,:,7)));
AntSMCount = nan(size(ANOVAAntRZSM(:,:,7)));
SMPulseCount = nan(size(ANOVAAntRZSM(:,:,7)));

AntRZSMCount(TTPZeros(:,:,5)>0.5 & ANOVAAntRZSM(:,:,7)<0.05)=1;
LAICount(TTPZeros(:,:,5)>0.5 & ANOVALAI(:,:,7)<0.05 & ANOVALAI(:,:,9)>0)=1;
% LAICount(TTPZeros(:,:,5)>0.5 & ANOVALAI(:,:,7)<0.05)=1;
AntVODCount(TTPZeros(:,:,5)>0.5 & ANOVAAntVOD(:,:,7)<0.05)=1;
AntSMCount(TTPZeros(:,:,5)>0.5 & ANOVAAntSM(:,:,7)<0.05)=1;
SMPulseCount(TTPZeros(:,:,5)>0.5 & ANOVASMPulse(:,:,7)<0.05)=1;

PerRZSM = nansum(AntRZSMCount(:))/nansum(TTPFlag(:));
PerLAI = nansum(LAICount(:))/nansum(TTPFlag(:));
PerAntVOD = nansum(AntVODCount(:))/nansum(TTPFlag(:));
PerAntSM = nansum(AntSMCount(:))/nansum(TTPFlag(:));
PerSMPulse = nansum(SMPulseCount(:))/nansum(TTPFlag(:));

% Calculate percentages for bins 1 and 2
AntRZSMCount12 = nan(size(ANOVAAntRZSM(:,:,5)));
LAICount12 = nan(size(ANOVAAntRZSM(:,:,5)));
AntVODCount12 = nan(size(ANOVAAntRZSM(:,:,5)));
AntSMCount12 = nan(size(ANOVAAntRZSM(:,:,5)));
SMPulseCount12 = nan(size(ANOVAAntRZSM(:,:,5)));

AntRZSMCount12(TTPZeros(:,:,5)>0.5 & ANOVAAntRZSM(:,:,5)<0.05)=1;
LAICount12(TTPZeros(:,:,5)>0.5 & ANOVALAI(:,:,5)<0.05)=1;
AntVODCount12(TTPZeros(:,:,5)>0.5 & ANOVAAntVOD(:,:,5)<0.05)=1;
AntSMCount12(TTPZeros(:,:,5)>0.5 & ANOVAAntSM(:,:,5)<0.05)=1;
SMPulseCount12(TTPZeros(:,:,5)>0.5 & ANOVASMPulse(:,:,5)<0.05)=1;

PerRZSM12 = nansum(AntRZSMCount12(:))/nansum(TTPFlag(:));
PerLAI12 = nansum(LAICount12(:))/nansum(TTPFlag(:));
PerAntVOD12 = nansum(AntVODCount12(:))/nansum(TTPFlag(:));
PerAntSM12 = nansum(AntSMCount12(:))/nansum(TTPFlag(:));
PerSMPulse12 = nansum(SMPulseCount12(:))/nansum(TTPFlag(:));

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

% ANOVA mat description
% 1) full ANOVA p value
% 2) Group 1 minus 2 mean difference
% 3) Group 1 minus 3 mean difference
% 4) Group 2 minus 3 mean difference
% 5) Group 1 minus 2 mean difference p value
% 6) Group 1 minus 3 mean difference p value
% 7) Group 2 minus 3 mean difference p value

perZeroFlag = TTPZeros(:,:,5);
perZeroFlag(perZeroFlag<0.5)=nan;
perZeroFlag(isfinite(perZeroFlag))=1;

SMPulseSeasMeanDelay = SMPulseSeasMean.*perZeroFlag;
SMPulseSeasNShortDelay = SMPulseSeasNShort.*perZeroFlag;
SMPulseSeasNLongDelay = SMPulseSeasNLong.*perZeroFlag;
SMPulseSeasNRapidDelay = SMPulseSeasNRapid.*perZeroFlag;
SMPulseN = SMPulseSeasNShortDelay+SMPulseSeasNLongDelay+SMPulseSeasNRapidDelay;

SMPulseSeasMeanDelay(SMPulseN<30)=nan;
SMPulseSeasNShortDelay(SMPulseN<30)=nan;
SMPulseSeasNLongDelay(SMPulseN<30)=nan;
SMPulseSeasNRapidDelay(SMPulseN<30)=nan;
SMPulseN(SMPulseN<30)=nan;

SMPulseSeasNRapidDelayPer = SMPulseSeasNRapidDelay./SMPulseN;
SMPulseSeasNShortDelayPer = SMPulseSeasNShortDelay./SMPulseN;
SMPulseSeasNLongDelayPer = SMPulseSeasNLongDelay./SMPulseN;

ShortPerSeas = nan(3,6);
RapidPerSeas = nan(3,6);
LongPerSeas = nan(3,6);
NSeas = nan(1,6);
MeanDistSeas = nan(6,5000);

for i = 1:6
A = SMPulseSeasNRapidDelayPer(:,:,i); A=A(:); A(isnan(A))=[];
RapidPerSeas(1,i) = quantile(A,0.25);
RapidPerSeas(2,i) = quantile(A,0.5);
RapidPerSeas(3,i) = quantile(A,0.75);
A = SMPulseSeasNShortDelayPer(:,:,i); A=A(:); A(isnan(A))=[];
ShortPerSeas(1,i) = quantile(A,0.25);
ShortPerSeas(2,i) = quantile(A,0.5);
ShortPerSeas(3,i) = quantile(A,0.75);
A = SMPulseSeasNLongDelayPer(:,:,i); A=A(:); A(isnan(A))=[];
LongPerSeas(1,i) = quantile(A,0.25);
LongPerSeas(2,i) = quantile(A,0.5);
LongPerSeas(3,i) = quantile(A,0.75);
% A = SMPulseN(:,:,i); A=A(:); A(isnan(A))=[];
% NSeas(i) = nansum(A);
A = SMPulseN(:,:,i); A=A(:); A(isnan(A))=[];
NSeas(i) = quantile(A,0.5);
A = SMPulseSeasMeanDelay(:,:,i); A=A(:); A(isnan(A))=[];
MeanDistSeas(i,1:length(A))=A;
end
edge = -150:60:150;


ShortPerSeasN = nan(3,6);
RapidPerSeasN = nan(3,6);
LongPerSeasN = nan(3,6);

for i = 1:6
A = SMPulseSeasNRapidDelay(:,:,i); A=A(:); A(isnan(A))=[];
RapidPerSeasN(1,i) = quantile(A,0.25);
RapidPerSeasN(2,i) = quantile(A,0.5);
RapidPerSeasN(3,i) = quantile(A,0.75);
A = SMPulseSeasNShortDelay(:,:,i); A=A(:); A(isnan(A))=[];
ShortPerSeasN(1,i) = quantile(A,0.25);
ShortPerSeasN(2,i) = quantile(A,0.5);
ShortPerSeasN(3,i) = quantile(A,0.75);
A = SMPulseSeasNLongDelay(:,:,i); A=A(:); A(isnan(A))=[];
LongPerSeasN(1,i) = quantile(A,0.25);
LongPerSeasN(2,i) = quantile(A,0.5);
LongPerSeasN(3,i) = quantile(A,0.75);
end

SMPulseSeasMeanDelayLAI = SMPulseSeasMeanLAI.*perZeroFlag;
SMPulseSeasNShortDelayLAI = SMPulseSeasNShortLAI.*perZeroFlag;
SMPulseSeasNLongDelayLAI = SMPulseSeasNLongLAI.*perZeroFlag;
SMPulseSeasNRapidDelayLAI = SMPulseSeasNRapidLAI.*perZeroFlag;
SMPulseNLAI = SMPulseSeasNShortDelayLAI+SMPulseSeasNLongDelayLAI+SMPulseSeasNRapidDelayLAI;

SMPulseSeasMeanDelayLAI(SMPulseNLAI<30)=nan;
SMPulseSeasNShortDelayLAI(SMPulseNLAI<30)=nan;
SMPulseSeasNLongDelayLAI(SMPulseNLAI<30)=nan;
SMPulseSeasNRapidDelayLAI(SMPulseNLAI<30)=nan;
SMPulseNLAI(SMPulseNLAI<30)=nan;

SMPulseSeasNRapidDelayPerLAI = SMPulseSeasNRapidDelayLAI./SMPulseNLAI;
SMPulseSeasNShortDelayPerLAI = SMPulseSeasNShortDelayLAI./SMPulseNLAI;
SMPulseSeasNLongDelayPerLAI = SMPulseSeasNLongDelayLAI./SMPulseNLAI;

ShortPerSeasLAI = nan(3,6);
RapidPerSeasLAI = nan(3,6);
LongPerSeasLAI = nan(3,6);
NSeasLAI = nan(1,6);
MeanDistSeasLAI = nan(6,5000);

for i = 1:6
A = SMPulseSeasNRapidDelayPerLAI(:,:,i); A=A(:); A(isnan(A))=[];
RapidPerSeasLAI(1,i) = quantile(A,0.25);
RapidPerSeasLAI(2,i) = quantile(A,0.5);
RapidPerSeasLAI(3,i) = quantile(A,0.75);
A = SMPulseSeasNShortDelayPerLAI(:,:,i); A=A(:); A(isnan(A))=[];
ShortPerSeasLAI(1,i) = quantile(A,0.25);
ShortPerSeasLAI(2,i) = quantile(A,0.5);
ShortPerSeasLAI(3,i) = quantile(A,0.75);
A = SMPulseSeasNLongDelayPerLAI(:,:,i); A=A(:); A(isnan(A))=[];
LongPerSeasLAI(1,i) = quantile(A,0.25);
LongPerSeasLAI(2,i) = quantile(A,0.5);
LongPerSeasLAI(3,i) = quantile(A,0.75);
% A = SMPulseN(:,:,i); A=A(:); A(isnan(A))=[];
% NSeas(i) = nansum(A);
A = SMPulseNLAI(:,:,i); A=A(:); A(isnan(A))=[];
NSeasLAI(i) = quantile(A,0.5);
A = SMPulseSeasMeanDelayLAI(:,:,i); A=A(:); A(isnan(A))=[];
MeanDistSeasLAI(i,1:length(A))=A;
end
edge = -150:60:150;


 % MODIS LAI
SMPulseSeasMeanDelayLAIMODIS = SMPulseSeasMeanLAIMODIS.*perZeroFlag;
SMPulseSeasNShortDelayLAIMODIS = SMPulseSeasNShortLAIMODIS.*perZeroFlag;
SMPulseSeasNLongDelayLAIMODIS = SMPulseSeasNLongLAIMODIS.*perZeroFlag;
SMPulseSeasNRapidDelayLAIMODIS = SMPulseSeasNRapidLAIMODIS.*perZeroFlag;
SMPulseNLAIMODIS = SMPulseSeasNShortDelayLAIMODIS+SMPulseSeasNLongDelayLAIMODIS+SMPulseSeasNRapidDelayLAIMODIS;

SMPulseSeasMeanDelayLAIMODIS(SMPulseNLAIMODIS<30)=nan;
SMPulseSeasNShortDelayLAIMODIS(SMPulseNLAIMODIS<30)=nan;
SMPulseSeasNLongDelayLAIMODIS(SMPulseNLAIMODIS<30)=nan;
SMPulseSeasNRapidDelayLAIMODIS(SMPulseNLAIMODIS<30)=nan;
SMPulseNLAIMODIS(SMPulseNLAIMODIS<30)=nan;

SMPulseSeasNRapidDelayPerLAIMODIS = SMPulseSeasNRapidDelayLAIMODIS./SMPulseNLAIMODIS;
SMPulseSeasNShortDelayPerLAIMODIS = SMPulseSeasNShortDelayLAIMODIS./SMPulseNLAIMODIS;
SMPulseSeasNLongDelayPerLAIMODIS = SMPulseSeasNLongDelayLAIMODIS./SMPulseNLAIMODIS;

ShortPerSeasLAIMODIS = nan(3,6);
RapidPerSeasLAIMODIS = nan(3,6);
LongPerSeasLAIMODIS = nan(3,6);
NSeasLAIMODIS = nan(1,6);
MeanDistSeasLAIMODIS = nan(6,5000);

for i = 1:6
A = SMPulseSeasNRapidDelayPerLAIMODIS(:,:,i); A=A(:); A(isnan(A))=[];
RapidPerSeasLAIMODIS(1,i) = quantile(A,0.25);
RapidPerSeasLAIMODIS(2,i) = quantile(A,0.5);
RapidPerSeasLAIMODIS(3,i) = quantile(A,0.75);
A = SMPulseSeasNShortDelayPerLAIMODIS(:,:,i); A=A(:); A(isnan(A))=[];
ShortPerSeasLAIMODIS(1,i) = quantile(A,0.25);
ShortPerSeasLAIMODIS(2,i) = quantile(A,0.5);
ShortPerSeasLAIMODIS(3,i) = quantile(A,0.75);
A = SMPulseSeasNLongDelayPerLAIMODIS(:,:,i); A=A(:); A(isnan(A))=[];
LongPerSeasLAIMODIS(1,i) = quantile(A,0.25);
LongPerSeasLAIMODIS(2,i) = quantile(A,0.5);
LongPerSeasLAIMODIS(3,i) = quantile(A,0.75);
A = SMPulseNLAIMODIS(:,:,i); A=A(:); A(isnan(A))=[];
NSeasLAIMODIS(i) = quantile(A,0.5);
A = SMPulseSeasMeanDelayLAIMODIS(:,:,i); A=A(:); A(isnan(A))=[];
MeanDistSeasLAIMODIS(i,1:length(A))=A;
end
edge = -150:60:150;


%%
LAIChangeYesGrowth = TTPGrowth(:,:,3).*TTPFlag;
LAIChangeNoGrowth = TTPGrowth(:,:,4).*TTPFlag;
[p] = ranksum(LAIChangeYesGrowth(:),LAIChangeNoGrowth(:));

LAIChangeNoGrowth(isnan(LAIChangeNoGrowth))=[];
LAIChangeYesGrowth(isnan(LAIChangeYesGrowth))=[];

% %%
LAIBoxNoGrowth = quantile(LAIChangeNoGrowth(:),[0.05 0.25 0.5 0.75 0.95]);
LAIBoxYesGrowth = quantile(LAIChangeYesGrowth(:),[0.05 0.25 0.5 0.75 0.95]);

barwidth = 0.3;
barVec = 1+[-barwidth barwidth];
for i = 1:length(LAIBoxNoGrowth)
    if i == 1 | i==length(LAIBoxNoGrowth)
plot(barVec,ones(size(barVec))*LAIBoxNoGrowth(i),'-k','linewidth',1.5);
    else
plot(barVec,ones(size(barVec))*LAIBoxNoGrowth(i),'k','linewidth',2);        
    end
hold on
end

plot([barVec(1) barVec(1)],[LAIBoxNoGrowth(2) LAIBoxNoGrowth(4)],'k','linewidth',2)
hold on
plot([barVec(end) barVec(end)],[LAIBoxNoGrowth(2) LAIBoxNoGrowth(4)],'k','linewidth',2)
hold on
plot([barVec(1)+barwidth barVec(1)+barwidth],[LAIBoxNoGrowth(4) LAIBoxNoGrowth(5)],'--k','linewidth',2)
hold on
plot([barVec(1)+barwidth barVec(1)+barwidth],[LAIBoxNoGrowth(1) LAIBoxNoGrowth(2)],'--k','linewidth',2)
hold on
patch([barVec(1) barVec(1) barVec(2) barVec(2)],[LAIBoxNoGrowth(2) ...
    LAIBoxNoGrowth(4) LAIBoxNoGrowth(4) LAIBoxNoGrowth(2)],'r')
hold on
plot(barVec,ones(size(barVec))*LAIBoxNoGrowth(3),'-r','linewidth',3);
alpha(0.1)

hold on
barVec = 3+[-barwidth barwidth];
for i = 1:length(LAIBoxYesGrowth)
    if i == 1 | i==length(LAIBoxNoGrowth)
plot(barVec,ones(size(barVec))*LAIBoxYesGrowth(i),'-k','linewidth',1.5);        
    else
plot(barVec,ones(size(barVec))*LAIBoxYesGrowth(i),'k','linewidth',2); 
    end
hold on
end
plot([barVec(1) barVec(1)],[LAIBoxYesGrowth(2) LAIBoxYesGrowth(4)],'k','linewidth',2)
hold on
plot([barVec(end) barVec(end)],[LAIBoxYesGrowth(2) LAIBoxYesGrowth(4)],'k','linewidth',2)
hold on
plot([barVec(1)+barwidth barVec(1)+barwidth],[LAIBoxYesGrowth(4) LAIBoxYesGrowth(5)],'--k','linewidth',2)
hold on
plot([barVec(1)+barwidth barVec(1)+barwidth],[LAIBoxYesGrowth(1) LAIBoxYesGrowth(2)],'--k','linewidth',2)
hold on
patch([barVec(1) barVec(1) barVec(2) barVec(2)],[LAIBoxYesGrowth(2) ...
    LAIBoxYesGrowth(4) LAIBoxYesGrowth(4) LAIBoxYesGrowth(2)],'b')
alpha(0.1)
hold on
plot(barVec,ones(size(barVec))*LAIBoxYesGrowth(3),'-r','linewidth',3);

ylabel('Median {\itt_p} (Days)')
set(gca,'xtick',[1 3])
set(gca,'xticklabel',[{'\DeltaLAI/\Deltat<0'};{'\DeltaLAI/\Deltat>0'}])
xlim([0 4])
ylim([0 6])
set(gca,'fontsize',18)

%% LAI
Ncr = 200;
Ncb = 200;
% green1 = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;     
green1 = [ linspace(0,1,Ncb)'   linspace(0.7,1,Ncb)'  linspace(0,1,Ncb)'  ];     
brown1 = [ linspace(0.3,1,Ncb)'   linspace(0,1,Ncb)'  linspace(0,1,Ncb)'  ];     
% cmapRedBu = [red;flipud(blu)];
cmapGreBro = [brown1;flipud(green1)];

figure
% ax3 = subplot(2,4,[6 7]);
az5 = subplot(2,2,1);
% az5 = subplot(3,4,[9 10 11 12]);

LAIChangeTTPZero = ANOVALAI(:,:,8).*TTPFlag;
LAIChangeTTPShort = ANOVALAI(:,:,9).*TTPFlag;
LAIChangeTTPLong = ANOVALAI(:,:,10).*TTPFlag;

% LAILongPer = ANOVALAI3_Per(:,:,10).*TTPFlag;
% nanmean(LAILongPer(:))*7 

abox = boxplot([LAIChangeTTPZero(:) LAIChangeTTPShort(:) LAIChangeTTPLong(:)],...
'positions',[1:3],'Colors', [0 0 0; 1 0 0; 0 0 1],'Notch','on');
set(gca,'xticklabel','')
set(findobj(gca,'type','line'),'linew',1.5)
ylabel('\DeltaLAI/\Deltat (m^2 m^{-2} day^{-1})')

h = findobj(gca,'Tag','Box');
colors = flipud([0 0 0; 1 0 0; 0 0 1]);
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.1);
end
% grid on
set(gca,'fontsize',14)
h=findobj(gca,'tag','Outliers');
delete(h) 
ylim([-0.02 0.02])
hold on
plot([0 5],[0 0],'--k','linewidth',1)

xTick = [1 2 3];
set(gca,'xtick',xTick)
yTick = get(gca,'ytick');

text(1,-0.023,[{'{\itt_p}=0'}],'fontsize',14,'horizontalalignment','center')
text(2,-0.023,[{'1?{\itt_p}?3'}],'fontsize',14,'horizontalalignment','center')
text(3,-0.023,[{'{\itt_p}>3'}],'fontsize',14,'horizontalalignment','center')
text(2,-0.028,[{'Time to Peak Plant Water Content'}],'fontsize',14,'horizontalalignment','center')

% hLegend = legend(flipud(findall(gca,'Tag','Box')), ...
% {'Rapid Response ({\itt_p}=0)','Short VOD Increase (1?{\itt_p}?3)',...
% 'Long VOD Increase ({\itt_p}>3)'},'fontsize',14,'position',[0.228 0.07 0.1 0.1]);

% annotation('textarrow',[0.33 0.33],[0.59 0.75])
annotation('textarrow',[0.33 0.33],[0.77 0.85])
% text(2.5,0.015,[{'Growth'};{'Influences {\itt_p}'}],'fontsize',14,'horizontalalignment','center')
text(2.45,0.0155,[{'Growth'};{'Influences'}; {'VOD Increase'}],'fontsize',14,'horizontalalignment','center')
% annotation('textarrow',[0.33 0.33],[0.55 0.39])
annotation('textarrow',[0.33 0.33],[0.75 0.67])
% text(2.5,-0.015,[{'Only'};{'Rehydration'};{'Influences {\itt_p}'}],'fontsize',14,'horizontalalignment','center')
text(2.5,-0.0155,[{'Only'} ; {'Rehydration Influences'};{'VOD Increase'}],'fontsize',14,'horizontalalignment','center')

perZeroFlag = TTPZeros(:,:,5);
perZeroFlag(perZeroFlag<0.5)=nan;
perZeroFlag(isfinite(perZeroFlag))=1;
LAICorr = ANOVALAICorr(:,:,1).*perZeroFlag;
LAICorr0 = ANOVALAI0Corr(:,:,1).*perZeroFlag;

LAICorrPosSig = LAICorr0;
LAICorrPosSig(LAICorrPosSig<0 | ANOVALAICorr0(:,:,2)>0.05)=nan;
LAICorrPosSig(isfinite(LAICorrPosSig))=1;
LAICorrAntTot = LAICorr0;
LAICorrAntTot(isfinite(LAICorrAntTot))=1;
LAICorrPosSigPer0 = nansum(LAICorrPosSig(:))/nansum(LAICorrAntTot(:));

LAICorrPosSig = LAICorr;
LAICorrPosSig(LAICorrPosSig<0 | ANOVALAICorr(:,:,2)>0.05)=nan;
LAICorrPosSig(isfinite(LAICorrPosSig))=1;
LAICorrAntTot = LAICorr;
LAICorrAntTot(isfinite(LAICorrAntTot))=1;
LAICorrPosSigPer = nansum(LAICorrPosSig(:))/nansum(LAICorrAntTot(:));

az6 = subplot(2,2,2);
plot(edge,RapidPerSeas(2,:),'o-k','linewidth',2.5)
hold on
plot(edge,ShortPerSeas(2,:),'o-r','linewidth',2.5)
hold on
plot(edge,LongPerSeas(2,:),'o-','color',[0 0 1],'linewidth',2.5)
xlabel([{'Time After Seasonal'}; {'Soil Moisture Peak (Days)'}])
% ylim([0 0.5])
ylim([0 0.55])

% hLegend = legend(flipud(findall(gca,'Tag','Box')), ...
% {'Rapid Response','Short Delay',...
% 'Long Delay'},'fontsize',14,'position',[0.67 0.25 0.05 0.05]);
hLegend = legend(flipud(findall(gca,'Tag','Box')), ...
{'Rapid VOD Response ({\itt_p}=0)',...
'Short VOD Increase (1?{\itt_p}?3)',...
'Long VOD Increase ({\itt_p}>3)'},'fontsize',14,'position',[0.703 0.63 0.05 0.03]);
% hLegend = legend(flipud(findall(gca,'Tag','Box')), ...
% {'Short VOD Increase (1?{\itt_p}?3)',...
% 'Long VOD Increase ({\itt_p}>3)'},'fontsize',14,'position',[0.703 0.63 0.05 0.03]);

set(gca,'fontsize',14)
% set(gca,'xticklabel','')
xlim([-150 150])
ylabel('Fraction of Pulses')
% grid on


az1 = subplot(2,2,3);
pcolor(TopLeftCornerLon,TopLeftCornerLat,...
    ANOVALAI(:,:,9).*TTPFlag) ; shading flat
hold on
% title('\DeltaLAI per day (Short Delay)')
ylim([-35 38])
xlim([-20 55])
text(80,0,'\DeltaLAI/\Deltat (m^2 m^{-2} day^{-1})','rotation',270,...
    'fontsize',14,'horizontalalignment','center')
text(-3,-20,'1?{\itt_p}?3',...
    'fontsize',18,'horizontalalignment','center')
hold on
geoshow(clat,clon,'LineWidth',1,'Color','k')
colormap(cmapGreBro)
set(gca,'color',[0.8 0.8 0.8])
caxis([-0.02 0.02])
set(gca,'fontsize',14)
colorbar
set(gca,'xticklabel','')
set(gca,'yticklabel','')

az2 = subplot(2,2,4);
pcolor(TopLeftCornerLon,TopLeftCornerLat,...
    ANOVALAI(:,:,10).*TTPFlag) ; shading flat
hold on
% title('\DeltaLAI per day (Long Delay)')
text(80,0,'\DeltaLAI/\Deltat (m^2 m^{-2} day^{-1})','rotation',270,...
    'fontsize',14,'horizontalalignment','center')
text(-3,-20,'{\itt_p}>3',...
    'fontsize',18,'horizontalalignment','center')
ylim([-35 38])
xlim([-20 55])
hold on
geoshow(clat,clon,'LineWidth',1,'Color','k')
colormap(cmapGreBro)
set(gca,'color',[0.8 0.8 0.8])
caxis([-0.02 0.02])
set(gca,'fontsize',14)
colorbar
set(gca,'xticklabel','')
set(gca,'yticklabel','')

set(gcf,'position',[100 100 1000 550])


AMSA5 = get(az5,'Position'); % top right fig
set(az5,'Position',[AMSA5(1)-0.03 AMSA5(2)+0.01 AMSA5(3) AMSA5(4)]);
AMSA5 = get(az5,'Position'); % top right fig

AMSA6 = get(az6,'Position'); % top right fig
set(az6,'Position',[AMSA6(1) AMSA5(2) AMSA5(3) AMSA5(4)]);
AMSA6 = get(az6,'Position'); % top right fig

AMSA1 = get(az1,'Position'); % top right fig
set(az1,'Position',[AMSA5(1)-0.02 AMSA1(2) AMSA5(3)-0.04 AMSA5(4)]);
AMSA1 = get(az1,'Position'); % top right fig

AMSA2 = get(az2,'Position'); % top right fig
set(az2,'Position',[AMSA6(1)-0.02 AMSA1(2) AMSA1(3) AMSA1(4)]);
AMSA2 = get(az2,'Position'); % top right fig

set(gcf,'position',[10 10 800 650])

% AMSA6 = get(az6,'Position'); % top right fig
% set(az6,'Position',[AMSA6(1)-0.08 AMSA6(2)-0.0 AMSA6(3) AMSA6(4)]);
% AMSA6 = get(az6,'Position'); % top right fig

Atext = ['A'];
annotation('textbox',[0.02 0.86 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.48 0.86 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.02 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['D'];
annotation('textbox',[0.48 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
