%% Plot Figure 1 from Main_Figures_Data.mat
clear; close all;

script_dir = fileparts(mfilename('fullpath'));
source_dir = '/Users/jianzhidong/Documents/Research/C4transition/Remote_sensing_SIF';
addpath(source_dir);
figure_dir = fullfile(script_dir,'Generated_Figures');
if ~exist(figure_dir,'dir'); mkdir(figure_dir); end

S = load(fullfile(script_dir,'Main_Figures_Data.mat'),'Figure1');
F = S.Figure1;
lon1 = F.longitude;
lat1 = F.latitude;
kland_05 = F.land_index;
kclean = F.map_exclusion_index;
dsd_all = F.storage_change;
polygons = F.vegetation_boundaries;
dmeans = F.group_mean;
dstds = F.group_error;
dsums = F.group_fraction;
dsum_stds = F.group_fraction_error;

rd = viridis;
fig = figure;
width = 10.5;
height = 4.5;
set(gcf,'Position',[100 80 width*100 height*100]);
load coastlines.mat;

rd = brewermap(256,'OrRd');
tt = 1.8;
kk = round(linspace(1,256^tt,220).^(1/tt));

% Subplot 1: Map
h1 = subplot(3,1,1);
axesm('robinson','MapLatLimit',[-60 90],'Grid','off','Frame','off');
set(gca,'FontName','Arial','FontSize',12);
land = shaperead('landareas','UseGeoCoords',true);
geoshow(land,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
hold on;
A = nan(360,720);
A(kland_05) = dsd_all(:,1);
A(kclean) = nan;
A = double(A);
h = pcolorm(lat1,lon1,A);
caxis([-5 0]);
set(h,'EdgeColor','none');
colormap(h1,flipud(rd(kk,:)));
plotm(coastlat,coastlon,'Color',[0.25 0.25 0.25],'LineWidth',1);

for i = 1:length(polygons)
    a = polygons{i}.Vertices(:,1);
    b = polygons{i}.Vertices(:,2);
    a1 = nanmin(a):2:nanmax(a);
    b1 = nanmin(b):2:nanmax(b);
    [lon,lat] = meshgrid(a1,b1);
    IN = inpolygon(lon,lat,a,b);
    scatterm(lat(IN),lon(IN),1.5,'k','filled');
end

tightmap;
set(gca,'XTick',[],'YTick',[],'Box','off');
set(h1,'XColor','none','YColor','none','Color','none');
c1 = colorbar('southoutside');
xlabel(c1,'cm','FontSize',12,'FontName','Arial');
title('(a) E(TWSd | Z_p<0)','FontWeight','bold','FontName','Arial','FontSize',14);
set(gca,'FontSize',14,'FontName','Arial');
c1.TickLabels{end} = '\geq 0';
set(h1,'XColor','none','YColor','none','Color','none');

% Subplot 2: Mean TWSd
h2 = subplot(3,1,2);
hold on;
kkuse = 1:4;
group_shift = 0.225;
x = 1:4;
bar_width = 0.8;
hh = bar(x,dmeans(kkuse,:),'grouped','BarWidth',bar_width);
hh(1).FaceColor = [0.34 0.63 0.85];
hh(2).FaceColor = [0.45 0.75 0.45];
hh(3).FaceColor = [0.85 0.47 0.31];
errorbar(x-group_shift,dmeans(kkuse,1),dstds(kkuse,1),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
errorbar(x,dmeans(kkuse,2),dstds(kkuse,2),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
errorbar(x+group_shift,dmeans(kkuse,3),dstds(kkuse,3),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
grid on;
set(gca,'GridLineStyle','--','GridAlpha',0.3,'GridColor',[0.5 0.5 0.5]);
set(gca,'FontName','Arial','FontSize',12,'LineWidth',1);
ylim([-3 0]); xlim([0.5 4.5]);
xticks(x);
xticklabels({'Z_P<0','Z_P<-1','Z_P<0, VPD>1kpa','Z_P<0, VPD>1.5kpa'});
set(gca,'XTickLabelRotation',0);
ylabel('cm','FontSize',12,'FontName','Arial');
title('(b) Mean TWSd during droughts','FontWeight','bold','FontSize',14,'FontName','Arial');
set(gca,'YTick',-3:1:0);
box on;
legend({'Grass/Crop','Forest','Others'},'Location','SouthWest','FontSize',11,'Box','off','FontName','Arial');

% Subplot 3: Contribution
h3 = subplot(3,1,3);
hold on;
hh = bar(x,dsums(kkuse,:),'grouped','BarWidth',bar_width);
hh(1).FaceColor = [0.34 0.63 0.85];
hh(2).FaceColor = [0.45 0.75 0.45];
hh(3).FaceColor = [0.85 0.47 0.31];
errorbar(x-group_shift,dsums(kkuse,1),dsum_stds(kkuse,1),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
errorbar(x,dsums(kkuse,2),dsum_stds(kkuse,2),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
errorbar(x+group_shift,dsums(kkuse,3),dsum_stds(kkuse,3),'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
grid on;
set(gca,'GridLineStyle','--','GridAlpha',0.3,'GridColor',[0.5 0.5 0.5]);
set(gca,'FontName','Arial','FontSize',12,'LineWidth',1);
xlim([0.5 4.5]);
xticks(x);
xticklabels({'Z_P<0','Z_P<-1','Z_P<0, VPD>1kpa','Z_P<0, VPD>1.5kpa'});
set(gca,'XTickLabelRotation',0);
ylabel('Fraction','FontSize',12,'FontName','Arial');
title('(c) Contribution to global TWSd depletion','FontWeight','bold','FontSize',14,'FontName','Arial');
box on;

set(h1,'Position',[-0.22 0.15 0.9 0.62]);
set(c1,'Position',[0.03 0.15 0.45 0.015]);
set(h2,'Position',[0.57 0.55 0.42 0.26]);
set(h3,'Position',[0.57 0.15 0.42 0.26]);

filename = fullfile(figure_dir,'Figure1_drought_impacts_TWS_upd.png');
exportgraphics(fig,filename,'Resolution',300);
fprintf('Created: %s\n',filename);

