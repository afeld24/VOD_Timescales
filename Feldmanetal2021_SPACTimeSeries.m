%%
clc
clear

cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy_TauEfolding/Figures_ToDistribute')
load('SPACTimeSeries_plantonly')
tpredawn = -5:1:8;
tpredawn1 = tpredawn-0.5;
dpsiWPrct = prctile(dpsiWPredawnVec,[25 50 75],1);

figure
iplotA = subplot(2,2,1);

patch([-1 0 0 -1],[-2 -2 2 2],[0.5 0.5 0.5],...
    'facealpha',0.5,'linestyle','none')
hold on
binwidth = 0.125;
for ip1 = 1:length(tpredawn1)
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(3,ip1) dpsiWPrct(3,ip1)],'-k')   
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(1,ip1)],'-k')  
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k')     
    hold on
    plot([tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k') 
    hold on
    patch([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],...
    [dpsiWPrct(1,ip1) dpsiWPrct(3,ip1) dpsiWPrct(3,ip1) dpsiWPrct(1,ip1)],...
    [0 0 0.8],'facealpha',0.1,'linestyle','none')
end
plot(tpredawn1,dpsiWPrct([2],:),'.-k','linewidth',2,'markersize',30)
hold on
plot([-1 15],[0 0],'--k','linewidth',1)
set(gca,'xtick',[-1:1:8])
xlim([-1 8])
ylim([-2 2])
ylabel('d\psi_{W}/dt (MPa/Day)')
xlabel('Time After Pulse (Day)')
grid on
set(gca,'fontsize',18)

dim1 = [0.14 0.56 0.1 0.1];
str1 = ['R_{P} Change Factor = ' ...
    sprintf('%0.1f',dist_rXFactor(1)) ' to ' ...
    sprintf('%0.1f',dist_rXFactor(2)) ''];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none') 
dim1 = [0.14 0.53 0.1 0.1];
str1 = ['R_{P} = ' ...
    sprintf('%0.1e',10.^dist_rX(1)) ' to ' ...
    sprintf('%0.1e',10.^dist_rX(2)) ' MPa/(m/s)'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')        
dim1 = [0.14 0.59 0.1 0.1];
str1 = ['Initial Rootzone \psi_{S} = ' ...
    sprintf('%0.2f',psiS_Init) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
dim1 = [0.14 0.62 0.1 0.1];
str1 = ['Surface \psi_{S} After Pulse = ' ...
    sprintf('%0.2f',psiS_Pulse) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')

text(2.5,1.5,['Plant Resistance Only'],'fontsize',18,'Fontname','arial')

load('SPACTimeSeries_soilonly')
tpredawn = -5:1:8;
tpredawn1 = tpredawn-0.5;
dpsiWPrct = prctile(dpsiWPredawnVec,[25 50 75],1);

iplotB = subplot(2,2,2);
patch([-1 0 0 -1],[-2 -2 2 2],[0.5 0.5 0.5],...
    'facealpha',0.5,'linestyle','none')
hold on
binwidth = 0.125;
for ip1 = 1:length(tpredawn1)
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(3,ip1) dpsiWPrct(3,ip1)],'-k')   
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(1,ip1)],'-k')  
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k')     
    hold on
    plot([tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k') 
    hold on
    patch([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],...
    [dpsiWPrct(1,ip1) dpsiWPrct(3,ip1) dpsiWPrct(3,ip1) dpsiWPrct(1,ip1)],...
    [0 0 0.8],'facealpha',0.1,'linestyle','none')
end
plot(tpredawn1,dpsiWPrct([2],:),'.-k','linewidth',2,'markersize',30)
hold on
plot([-1 15],[0 0],'--k','linewidth',1)
set(gca,'xtick',[-1:1:8])

xlim([-1 8])
ylim([-2 2])
ylabel('d\psi_{W}/dt (MPa/Day)')
xlabel('Time After Pulse (Day)')
grid on
set(gca,'fontsize',18)

dim1 = [0.6 0.53 0.1 0.1];
str1 = ['R_{P} = ' ...
    sprintf('%0.1e',10.^dist_rX(1)) ' to ' ...
    sprintf('%0.1e',10.^dist_rX(2)) ' MPa/(m/s)'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
dim1 = [0.6 0.56 0.1 0.1];
str1 = ['R_{P} Change Factor = ' ...
    sprintf('%0.1f',dist_rXFactor(1)) ' to ' ...
    sprintf('%0.1f',dist_rXFactor(2)) ''];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')      
dim1 = [0.6 0.59 0.1 0.1];
str1 = ['Initial Rootzone \psi_{S} = ' ...
    sprintf('%0.1f',dist_psiS_Init(1)) ' to ' ...
    sprintf('%0.1f',dist_psiS_Init(2)) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
dim1 = [0.6 0.62 0.1 0.1];
str1 = ['Surface \psi_{S} After Pulse = ' ...
    sprintf('%0.1f',dist_psiS_Pulse(1)) ' to ' ...
    sprintf('%0.1f',dist_psiS_Pulse(2)) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
    
text(2.5,1.5,['Soil Resistance Only'],'fontsize',18,'Fontname','arial')
    
load('SPACTimeSeries_soilplant')
tpredawn = -5:1:8;
tpredawn1 = tpredawn-0.5;
dpsiWPrct = prctile(dpsiWPredawnVec,[25 50 75],1);

iplotC = subplot(2,2,3);
patch([-1 0 0 -1],[-2 -2 2 2],[0.5 0.5 0.5],...
    'facealpha',0.5,'linestyle','none')
hold on
binwidth = 0.125;
for ip1 = 1:length(tpredawn1)
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(3,ip1) dpsiWPrct(3,ip1)],'-k')   
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(1,ip1)],'-k')  
    hold on
    plot([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k')     
    hold on
    plot([tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],[dpsiWPrct(1,ip1) dpsiWPrct(3,ip1)],'-k') 
    hold on
    patch([tpredawn1(ip1)-binwidth tpredawn1(ip1)-binwidth tpredawn1(ip1)+binwidth tpredawn1(ip1)+binwidth],...
    [dpsiWPrct(1,ip1) dpsiWPrct(3,ip1) dpsiWPrct(3,ip1) dpsiWPrct(1,ip1)],...
    [0 0 0.8],'facealpha',0.1,'linestyle','none')
end
plot(tpredawn1,dpsiWPrct([2],:),'.-k','linewidth',2,'markersize',30)
hold on
plot([-1 15],[0 0],'--k','linewidth',1)

xlim([-1 8])
ylim([-2 2])
ylabel('d\psi_{W}/dt (MPa/Day)')
xlabel('Time After Pulse (Day)')
grid on
set(gca,'fontsize',18)
set(gca,'xtick',[-1:1:8])
dim1 = [0.36 0.02 0.1 0.1];
str1 = ['R_{P} = ' ...
    sprintf('%0.1e',10.^dist_rX(1)) ' to ' ...
    sprintf('%0.1e',10.^dist_rX(2)) ' MPa/(m/s)'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
dim1 = [0.36 0.05 0.1 0.1];
str1 = ['R_{P} Change Factor = ' ...
    sprintf('%0.1f',dist_rXFactor(1)) ' to ' ...
    sprintf('%0.1f',dist_rXFactor(2)) ''];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')      
dim1 = [0.36 0.08 0.1 0.1];
str1 = ['Initial Rootzone \psi_{S} = ' ...
    sprintf('%0.1f',dist_psiS_Init(1)) ' to ' ...
    sprintf('%0.1f',dist_psiS_Init(2)) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')     
dim1 = [0.36 0.11 0.1 0.1];
str1 = ['Surface \psi_{S} After Pulse = ' ...
    sprintf('%0.1f',dist_psiS_Pulse(1)) ' to ' ...
    sprintf('%0.1f',dist_psiS_Pulse(2)) ' MPa'];
    annotation('textbox',dim1,'String',str1,'fontsize',14,...
        'Fontname','arial','EdgeColor','none')  
  
text(2.5,1.5,['Plant+Soil Resistance'],'fontsize',18,'Fontname','arial')
    
set(gcf,'position',[10 10 800 650])

AMSA = get(iplotA,'Position'); % top right fig
set(iplotA,'Position',[AMSA(1)-0.03 AMSA(2)-0. AMSA(3)+0.03 AMSA(4)+0.03]);
AMSA = get(iplotA,'Position'); % top right fig

AMSC = get(iplotB,'Position'); % top right fig
set(iplotB,'Position',[AMSC(1)+0.02 AMSA(2) AMSA(3) AMSA(4)]);
AMSC = get(iplotB,'Position'); % top right fig

AMSD = get(iplotC,'Position'); % top right fig
set(iplotC,'Position',[AMSD(1)+0.22 AMSD(2)-0.03 AMSA(3) AMSA(4)]);
AMSD = get(iplotC,'Position'); % top right fig

Atext = ['A'];
annotation('textbox',[0.02 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.51 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.26 0.39 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',28,'EdgeColor','none')
