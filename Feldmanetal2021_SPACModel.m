clc
clear

% Set directory
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy_TauEfolding/Figures_ToDistribute')

% Tuning variables
% 1) Root resistance
% 2) Plant capacitance
% 3) Stomatal resistance
% 4) Aerodynamic resistance
% 5) Initial psi_L
% 6) Initial psi_S
% 7) Soil moisture decay
% 8) VPD change

tdays = 8; %days
tdaysSpinUp = 3; %days
tdaysPreRain = 2; %days
tsteps = tdays*24; %hours
tstepsSpinUp = tdaysPreRain*24; %hours
% dt = 30; % minutes
% dt = 0.1; % minutes
dt = 0.01; % minutes
dtSave = dt;
% dt = 0.05; %time step (minutes/step)
% totaltime = 200; %total run time (hours)
% steps = totaltime*(60/dt); %number of steps
% tt = 0:dt:steps*dt; % minutes

minPerDay = 24*60;

% tt = 1:0.1:20; % days
tt = 0:dt:((tdays*minPerDay)-dt); % minutes
tt = tt./minPerDay; % days
ttSpinUp = 0:dt:((tdaysSpinUp*minPerDay)-dt); % minutes
ttSpinUp = ttSpinUp./minPerDay; % days
ttPreRain = 0:dt:((tdaysPreRain*minPerDay)-dt); % minutes
ttPreRain = ttPreRain./minPerDay; % days
ttFull = 0:dt:(((tdays+tdaysPreRain+tdaysSpinUp)*minPerDay)-dt); % minutes
ttFull = ttFull./minPerDay; % days
dt = dt./minPerDay; % days

% SMInit = 0.25; % Pulsed to this value
% SMInitPre = 0.18; % Start surface SM at this value
% SMInitDeep = 0.153; % Start rootzone SM at this value

psiS_Pulse = -0.01; % Pulsed to this value
psiS_Init = -0.5; % Full rootzone

clayAf = 0.2;
sandAf = 0.3;

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
SateffMat = SatMatSand.*sandAf+SatMatClay.*clayAf+SatMatLoam.*loamAf;

SMInit = neffMat*((psiS_Pulse/SateffMat)*((10^8)/9810)).^(-1/beffMat);
SMInitPre = neffMat*((psiS_Init/SateffMat)*((10^8)/9810)).^(-1/beffMat);
SMInitDeep = SMInitPre;

k = 10; % exponential decay (days)
kfact = 4; % Reduced decay of rootzone

smt1 = [SMInit*exp(-tt./k)]; %m3/m3
smtPreRain = [SMInitPre*exp(-ttPreRain./k)]; %m3/m3
smtRootPreRain = [SMInitDeep*exp(-ttPreRain./(k*kfact))]; %m3/m3
smtSpinUp = SMInitPre*ones(1,length(ttSpinUp),1); %m3/m3
smtRootSpinUp = SMInitDeep*ones(1,length(ttSpinUp)); %m3/m3

 [ttdays smZVec ttSave] = RichardsEq_InfFront_V1_FullRZ(k,kfact,...
    SMInit,smtRootPreRain(end),clayAf,sandAf,tsteps,dtSave);

%%%%% Random Draws %%%%%
% Set bounds of uniform distribution
dist_zrootzone = [0.1 0.9];
dist_VPD = [1 5];
dist_rX = [6 9];
dist_CW = [-6 -4];
dist_Wind = [1 8];
dist_psiW = [-4 -0.1];
dist_rXFactor = [-10 -1];
dist_rSFactor = [-10 -1];

% Draw parameter from distribution
draw_zrootzone = dist_zrootzone(1)+rand(1)*(dist_zrootzone(2)-dist_zrootzone(1));
draw_VPD = dist_VPD(1)+rand(1)*(dist_VPD(2)-dist_VPD(1));
draw_Wind = dist_Wind(1)+rand(1)*(dist_Wind(2)-dist_Wind(1));
draw_rX = 10^(dist_rX(1)+rand(1)*(dist_rX(2)-dist_rX(1)));
draw_CW = 10^(dist_CW(1)+rand(1)*(dist_CW(2)-dist_CW(1)));
draw_psiW = dist_psiW(1)+rand(1)*(dist_psiW(2)-dist_psiW(1));
draw_rXFactor = dist_rXFactor(1)+rand(1)*(dist_rXFactor(2)-dist_rXFactor(1));
draw_rSFactor = dist_rSFactor(1)+rand(1)*(dist_rSFactor(2)-dist_rSFactor(1));

% Fix Parameters
draw_rXFactor = -1;
draw_rSFactor = -3;
draw_CW = 10^(-6);
draw_rX = 1*10^(7);

psiLi = -3; %MPa ***
psiLCrit = -30;

draw_zrootzone = 0.4;

RCinithours = 1.5*((draw_rX*draw_CW)/3600)*3

% Based on rootzone, determine soil moisture series
zrootzone = [draw_zrootzone draw_zrootzone]; % (m) Get average rootzone moisture between bounds
kSave = k*24*60; % exponential decay (minutes)
kRoot = kSave*kfact;
% kRoot = k*kfact;
smtSave = [smt1 smt1(end)];
LayerInt = 0.05; % m
z = flipud((-1:LayerInt:0)');
[aRZ1] = find(abs(abs(z)-zrootzone(1)) == min(abs(abs(z)-zrootzone(1))));
[aRZ2] = find(abs(abs(z)-zrootzone(2)) == min(abs(abs(z)-zrootzone(2))));
smRZ = nanmean(smZVec(aRZ1:aRZ2,:),1);
tCross = find(abs(smRZ-smtSave) == min(abs(smRZ-smtSave)));
tCross = tCross(1);
ttCross = ttSave((tCross:length(ttSave))-tCross+1); % minutes
smtRoot1 = smRZ(1:tCross-1);
smtRoot2 = [smRZ(tCross)*exp(-ttCross./(kRoot))]; %m3/m3
smtRootZoneA = [smtRoot1 smtRoot2];

ttdays = ttdays(1:end-1);
smtRootZoneB = smtRootZoneA(1:end-1);

% Generate full soil moisture time series before and after spin up
smtRootZone = [smtRootSpinUp smtRootPreRain smtRootZoneB];
smt = [smtSpinUp smtPreRain smt1];

%%%% Initial Psi %%%%%
% Initial psiL

% Initial psiW
psiWi = draw_psiW; %MPa ***

CW = draw_CW; % m/MPa 

%%%%% Resistances %%%%%
rrinit = draw_rX;
rrfinal = draw_rX/abs(draw_rXFactor);
if rrfinal<(10^dist_rX(1))
    rrfinal = 10^dist_rX(1);
elseif rrfinal>(10^dist_rX(2))
    rrfinal = 10^dist_rX(2);
end

rr = [(rrinit)*ones(1,length(ttSpinUp)+length(ttPreRain)) ...
    linspace((rrinit),(rrfinal),length(ttSpinUp)) ...
    (rrfinal)*ones(1,length(tt)-length(ttSpinUp))]; % m/s
rX = rr; % Set rr = rX

TMultInit = 1; % transpiration multiplier
TMult = [(TMultInit)*ones(1,length(ttSpinUp)) (TMultInit)*ones(1,length(ttPreRain))...
    linspace(TMultInit,1,length(tt))]; % m/s
rW0 = rX; % MPa/(m/s) 

% Soil water potential and hydraulic conductivity (MPa)
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
SateffMat = SatMatSand.*sandAf+SatMatClay.*clayAf+SatMatLoam.*loamAf;
psiSSurface = ((SateffMat.*(smt./neffMat).^(-beffMat))/100).*(9810/10^6);
psiSRoot = ((SateffMat.*(smtRootZone./neffMat).^(-beffMat))/100).*(9810/10^6);

% Initial plant water content
RWC(1) = ((psiWi-psiLCrit)./(0-psiLCrit));

% Compute rW = f(RWC)
rW(1) = rW0(1)*(RWC(1).^(-1/4)); %Carlson and Lynn 1991

% Initial psiX
psiXinit = (psiSRoot(1)*rW(1)+psiWi*rX(1)+psiLi*rW(1))./(rX(1)+2*rW(1));

% Conditions After Rain
wind1 = linspace(draw_Wind,draw_Wind,length(tt)); % m/s
VPDlinDraw = linspace(draw_VPD,draw_VPD,length(tt)); % kPa

% Conditions Before Rain
windPreRain = linspace(draw_Wind,draw_Wind,length(ttPreRain)); % m/s
windSpinUp = linspace(draw_Wind,draw_Wind,length(ttSpinUp)); % m/s
VPDlinPreRain = linspace(draw_VPD,draw_VPD,length(ttPreRain)); % kPa
VPDlinSpinUp = linspace(draw_VPD,draw_VPD,length(ttSpinUp)); % kPa

wind = [windSpinUp windPreRain wind1];
VPDlin = [VPDlinSpinUp VPDlinPreRain VPDlinDraw];

A1 = sin((((1:(1./dt))./((1./dt)))*2*pi)-pi/2);%6AM
A = (1+repmat(A1,[1 tdays+tdaysPreRain+tdaysSpinUp]))/2;
VPD = VPDlin.*A;

%Psychrometric
cp = 1004; %J/(K kg)
MWeps1 = 0.622;
sp = 100; %kPa
rhoWater = 1000; %kg/m3
rhoAir = 1.2; %kg/m3

%Aerodynamic Resistance
h = 0.5; %m 
d = 0.75*h; %m
z0 = 0.1*h; %m
kVon = 0.41;
z = 2; %height at which wind is calculated
rA = (1./((kVon^2)*wind))*(log((z-d)./z0).^2); %s/m

%Stomatal Resistance
% Medlyn Model
LAI = 1; % m2/m2
g1 = 6; %kPa^1/2
g0 = -0.007; % mol/(m2*s)
Ca = 400; %umol/mol
PPreRain = linspace(4,4,length(ttPreRain)); % umolCo2/(m2*s)
PSpinUp = linspace(4,4,length(ttSpinUp)); % umolCo2/(m2*s)
P = linspace(6,4,length(tt)); % umolCo2/(m2*s)

P = [PSpinUp PPreRain P];
gs = g0+1.6.*(1+g1./sqrt(VPD)).*(P./Ca); % mol/(m2*s)

rs = 1./gs; % (m2*s)/mol
rs = rs.*(1/18).*(rhoAir*1000); % s/m ***check conversion factor
rS = rs./LAI; % s/m scale up to canopy scale

rSinitMult = abs(draw_rSFactor);
rSfinalMult = 1;
rSMult = [(rSinitMult)*ones(1,length(ttSpinUp)+length(ttPreRain)) ...
    linspace((rSinitMult),(rSfinalMult),length(ttSpinUp)) ...
    (rSfinalMult)*ones(1,length(tt)-length(ttSpinUp))]; % m/s
rS = rS.*rSMult;

psiL = nan(1,length(ttFull));
psiX = nan(1,length(ttFull));
psiW = nan(1,length(ttFull));
RWC = nan(1,length(ttFull));
rW = nan(1,length(ttFull));
dpsiW = nan(1,length(ttFull));

TSave = nan(1,length(ttFull));
RSave = nan(1,length(ttFull));
psiL(1) = psiLi;
psiX(1) = psiXinit;
psiW(1) = psiWi;

% Solve Numerically (Time Explicit)
for it = 1:length(ttFull)-1
% Convert to RWC (assuming linear relation between RWC and psiL)
RWC(it) = ((psiW(it)-psiLCrit)./(0-psiLCrit));

% Compute rW = f(RWC)
% Accounting for resistance if storage is dry
rW(it) = rW0(it)*(RWC(it).^(-1/4)); %Carlson and Lynn 1991
if isinf(rW(it))
    rW(it) = 1*10^20;
end

%1) Storage change
% dpsiL(it+1) = (1/CL)*(((psiX(it)-psiL(it))./rP)-T) * dt;
dpsiW(it+1) = ((psiX(it)-psiW(it))./(rW(it)*CW)) * (dt*86400);
psiW(it+1) = psiW(it)+dpsiW(it+1);
% if psiW(it+1)>psiX(it)
%     psiW(it+1) = psiX(it);
% end
if psiW(it+1)<psiLCrit
    psiW(it+1) = psiLCrit;
    dpsiW(it+1) = psiW(it+1)-psiW(it);
elseif psiW(it+1)>0
    psiW(it+1) = 0;
    dpsiW(it+1) = psiW(it+1)-psiW(it);        
end

% 2) Transpiration
% Setting a minimum leaf threshold
if psiL(it)>psiLCrit
T = (MWeps1./sp)*(rhoAir/rhoWater)*(VPD(it)./(rA(it)+rS(it)))*86400; %m/day See Margulis 
                    % example 8.3.1 -> 1000 is density of water
else
T = 0;
% rS(it+1:it+100)= 10^11;
end

T = T*TMult(it);
 
% 3) Leaf water change
psiL(it+1) = psiX(it)-(T/86400)*rX(it);
% Maintain critical leaf water potential
if psiL(it+1)<psiLCrit
    psiL(it+1) = psiLCrit;
end

% 4) Xylem water change
% psiX(it+1) = ((psiSRoot(it+1)*rW(it))+(psiW(it+1)*rX)+(psiL(it+1)*rW(it)))./...
%               (rX+2*rW(it));
psiX(it+1) = ((psiW(it+1)./rW(it))+(psiSRoot(it+1)./(rr(it+1)+rX(it+1)))+(psiL(it+1)./rX(it+1)))./...
       ((1./rW(it))+(1./(rr(it+1)+rX(it+1)))+(1./rX(it+1)));

      
if psiX(it+1)<psiLCrit
    psiX(it+1) = psiLCrit;
elseif psiX(it+1)>0
    psiX(it+1) = 0;
end

% if psiX(it+1)<psiL(it+1)
%     psiX(it+1) = psiL(it+1);
% end

TSave(it) = T;
if (psiSRoot(it)-psiX(it))>0
% RSave(it) = ((psiSRoot(it)-psiX(it))./rX(it))* 86400; %m/s
RSave(it) = ((psiSRoot(it)-psiX(it))./(rr(it+1)+rX(it+1)))* 86400; %m/s
else
RSave(it) = 0;    
end
end

TSavems = TSave;
TSave = (TSave.*1000); %m/day to mm/day
RSave = (RSave.*1000); %m/day to mm/day
dW = (RSave-TSave)*dt;

% dt = dtSave; % minutes
ttPlot1 = (-1*(tdaysPreRain+tdaysSpinUp)*minPerDay):dtSave:(((tdays)*minPerDay)-dtSave); % minutes
ttPlot = ttPlot1./minPerDay; % days

dayVec = [-(tdaysPreRain+tdaysSpinUp):1:tdays];
psiWmaxVec = nan(length(dayVec)-1,2);
for iday = dayVec(1):1:dayVec(end)-1
    [astart] = find(abs(ttPlot-(iday)) == min(abs(ttPlot-(iday))));
    [aend] = find(abs(ttPlot-(iday+1)) == min(abs(ttPlot-(iday+1))));
    da = round((aend-astart)/2);
    
    if iday == dayVec(1)
        astart1 = astart;
        aend1 = round(aend/2);  
    else
        astart1 = astart-da;
        aend1 = astart+da;
    end
    
    tfindVec = astart1:aend1;
[aSplit] = find(psiW(tfindVec)==max(psiW(tfindVec))); aSplit = aSplit(1);
    psiWmaxVec(iday+tdaysPreRain+tdaysSpinUp+1,1) = psiW(tfindVec(aSplit));
    psiWmaxVec(iday+tdaysPreRain+tdaysSpinUp+1,2) = tfindVec(aSplit);    
end

% equilibration psiW initial
psiWInit = psiWmaxVec(4,1);
psiWmaxVec = psiWmaxVec(5:end,:); % Don't include spin up period

% Version 3 tp
[psiWmax] = find(psiWmaxVec(:,1) == max(psiWmaxVec(:,1)));
if psiWmax>2
    while (psiWmaxVec(psiWmax,1)-psiWmaxVec(psiWmax-1,1))./abs(psiWmaxVec(2,1))<0.1 & psiWmax>2
        psiWmax = psiWmax-1;
    end   
end

StoreUpTime = ttPlot(psiWmaxVec(psiWmax,2));
StoreUpTime(StoreUpTime>=0&StoreUpTime<1)=1; % First day effects do not count as delay
StoreUpTime = StoreUpTime-1;
StoreUpTime(StoreUpTime<-1)=-1;

A = (diff(psiWmaxVec(:,1),1)./(abs(psiWmaxVec(2,1))))*100;

figure
plot(ttPlot,psiSSurface,'--k','linewidth',2)
hold on
plot(ttPlot,psiSRoot,'-k','linewidth',6)
hold on
plot(ttPlot,psiW,'-b','linewidth',6)
hold on
plot(ttPlot,psiX,'-','linewidth',2,'color',[0 0.7 0])
hold on
plot(ttPlot,psiL,'-r','linewidth',1)
legend('Surface Soil','Root Soil','Storage','Xylem','Leaf','location','Southeast')

ylabel('\psi (MPa)')
set(gca,'fontsize',18)
ylim([-6 0])
xlabel('Time After Rain Pulse (Days)')
xlim([-(tdaysPreRain) ttPlot(end)])
set(gcf,'position',[100 100 1000 400])
grid on

dim1 = [0.15 0.35 0.1 0.1];
str1 = ['VPD = ' ...
    sprintf('%0.2f',draw_VPD) ' kPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none') 
dim1 = [0.15 0.30 0.1 0.1];
str1 = ['Root Depth = ' ...
    sprintf('%0.2f',draw_zrootzone) ' m'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none') 
dim1 = [0.15 0.25 0.1 0.1];
str1 = ['R_{Plant} Change Factor = ' ...
    sprintf('%0.2f',draw_rXFactor) ];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')   
dim1 = [0.15 0.2 0.1 0.1];
str1 = ['R_{Plant} = ' ...
    sprintf('%0.2e',draw_rX) ' MPa/(m/s)'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')    
dim1 = [0.15 0.15 0.1 0.1];
str1 = ['C_{Plant} = ' ...
    sprintf('%0.2e',draw_CW) ' m/MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')    
dim1 = [0.15 0.1 0.1 0.1];
str1 = ['\psi_{W} Init = ' ...
    sprintf('%0.2f',psiWInit) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')   
dim1 = [0.15 0.40 0.1 0.1];
str1 = ['RS Change Factor = ' ...
    sprintf('%0.2f',draw_rSFactor)];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')  
dim1 = [0.15 0.45 0.1 0.1];
str1 = ['Time to Peak \psi_{W} = ' ...
    sprintf('%0.2f',StoreUpTime(1)) ' Days'];
    annotation('textbox',dim1,'String',str1,'fontsize',16,...
        'Fontname','arial','EdgeColor','none')    
