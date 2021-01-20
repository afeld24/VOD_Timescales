function [ttdays smZVec tt] = ...
    RichardsEq_InfFront_V1_FullRZ(k,kfact,...
    SMInit,SMInitDeep,clayAf,sandAf,totaltime,dt)

% dt = dtSave;
% totaltime = tsteps;
% SMInitDeep = smtRootPreRain(end);

% zrootzone = [0.2 0.6];

% dt = 0.05; %time step (minutes/step)
% totaltime = 200; %total run time (hours)
steps = totaltime*(60/dt); %number of steps
tt = 0:dt:steps*dt; % minutes
k = k*24*60; % exponential decay (minutes)
kRoot = k*kfact;
ttdays = tt/(60*24);
ddtLength = 1100000;

% hrs*dt/(24*60)
% tt = 0:0.01:20
% NumEvents = 1;
% SMInit = 0.3;
% SMInitDeep = 0.15;
% k = 15;
smt = [SMInit*exp(-tt./k)]; %m3/m3
% plot(tt,smt)

LayerInt = 0.05; % m
% LayerInt = 0.01; % m
% Initial Vertical profile
% z = flipud((-0.4:0.005:0)');
z = flipud((-1:LayerInt:0)');
kInf = 3;
% smZ = SMInitDeep+(SMInit-SMInitDeep)*exp(z./0.01);
smZVec = nan(length(z),length(tt));
smZVec(:,1) = SMInitDeep;
[z05] = find(abs(z+0.05)==min(abs(z+0.05)));
z05 = z05(1);
smZVec(1:z05,:) = repmat(smt,[z05 1]);
% plot(smZ,z)

% clayAf = 0.2;
% sandAf = 0.4;
SatMatSand = -12.1; SatMatClay = -40.5; SatMatLoam = -47.8;
bSand = 4.05;       bClay = 11.4;       bLoam = 5.39;
nSand = 0.395;      nClay = 0.482;      nLoam = 0.451;
KSSand = 63.36;     KSClay = 0.46;      KSLoam = 2.5;
% smVec1 = min(ssmv(:)):0.01:max(ssmv(:));
loamAf = 1-clayAf-sandAf;
beffMat = bSand.*sandAf+bClay.*clayAf+bLoam.*loamAf;
neffMat = nSand.*sandAf+nClay.*clayAf+nLoam.*loamAf;
SateffMat = SatMatSand.*sandAf+SatMatClay.*clayAf+SatMatLoam.*loamAf;
SatKeffMat = KSSand.*sandAf+KSClay.*clayAf+KSLoam.*loamAf;
% psiS = ((SateffMat.*(smt./neffMat).^(-beffMat))/100).*(9810/10^6); %MPa
ceffMat = 2*beffMat+3;
% KS = ((SatKeffMat.*(smt./neffMat).^(ceffMat))/100)/3600; %m/s
psiSat = SateffMat/100; %m
KSat = (SatKeffMat/100)/3600; %m/s

% % Check
% m = -10; % 1/m
% z1 = -0.4; % m
% AA1 = ((beffMat*m*psiSat/neffMat)*((SMInit+m*z1)/neffMat).^(-beffMat-1)+1);
% AA2 = (1*KSat*((SMInit+m*z1)/neffMat).^(ceffMat))*100*3600; 
% qz = -AA1*AA2; %cm/hr

% c= 0;
% for ddt = 1:hrs-1
for ddt = 1:ddtLength
% for ddt = 1:400
%         c = c+1;
AVec = nan(length(z),1);
for iz = z05+1:length(z)
% for iz = length(z)-1:length(z)-1
    itheta = smZVec(iz,ddt);
%     psiS = (SateffMat.*(itheta./neffMat).^(-beffMat))/100; %m
%     psiS = ((SateffMat.*(itheta./neffMat).^(-beffMat))/100).*(9810/10^6); %m    
    K = KSat*(itheta./neffMat)^ceffMat; % m/s
    dpsidtheta = ((-beffMat)./neffMat)*psiSat*...
        (itheta./neffMat)^(-beffMat-1); %m/(m3/m3)
    dthetadz = (smZVec(iz-1,ddt)-smZVec(iz,ddt))./((z(iz-1)-z(iz))); %(m3/m3)/m
    AVec(iz) = K*(dpsidtheta*dthetadz+1); %m/s
end
% plot(AVec,z)
% set(gca, 'XScale', 'log')
dAdz = [diff(AVec)./diff(z); nan];
% Assuming exfiltration is zero
dAdz(isnan(dAdz))=0;
dAdz(dAdz<0)=0; %NEED TO CHANGE THIS
dtheta = dAdz*dt*60; 
smNew = smZVec(:,ddt)+dtheta;
% Don't let richards equation not conserve mass of infiltration front
smNew(smNew>smZVec(1,ddt+1))= smZVec(1,ddt+1); %NEED TO CHANGE THIS

smZVec(:,ddt+1) = smNew;

%%%% Adds Decay In Rootzone %%%%
if ddt ~= ddtLength  
dsm = smZVec(:,ddt+1)-smZVec(:,ddt);
[dsmZeroFind] = find(dsm==0);
[dsmNoZeroFind] = find(dsm~=0);
dsmZeroFind(dsmZeroFind <= dsmNoZeroFind(end)) =[];
smZPrev = smZVec(:,ddt)*exp(-dt./(kRoot));
smZVec(dsmZeroFind,ddt+1) = smZPrev(dsmZeroFind,:);
end
% smZVec(smZVec>smZVec(1,ddt+1))=smZVec(1,ddt+1); %NEED TO CHANGE THIS

% smZVec(end-5:end,ddt+1) = smt(ddt); % boundary condition
% smZVec(end,ddt+1) = SMInit; % boundary condition
end
