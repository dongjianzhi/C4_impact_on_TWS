%% Plot Figure 3 from Main_Figures_Data.mat
clear; close all;

script_dir = fileparts(mfilename('fullpath'));
source_dir = '/Users/jianzhidong/Documents/Research/C4transition/Remote_sensing_SIF';
addpath(source_dir);
figure_dir = fullfile(script_dir,'Generated_Figures');
if ~exist(figure_dir,'dir'); mkdir(figure_dir); end

S = load(fullfile(script_dir,'Main_Figures_Data.mat'),'Figure3');
F = S.Figure3;
lon1 = F.longitude;
lat1 = F.latitude;
kland_05 = F.land_index;
c4_diff_us = F.vegetation_fraction_change;
S_change = F.storage_change;
means = F.group_mean;
errors_lower = F.group_error_lower;
errors_upper = F.group_error_upper;
co2 = F.concentration;
ad = F.relative_response;

load coastlines.mat;
close all;
fig = figure;
width = 9.5;
height = 6.5;
set(gcf,'Position',[100 80 width*100 height*100],'Color','w', ...
    'PaperUnits','inches','PaperSize',[width height]);

color_use = cmocean('tarn',256);
color_use2 = cmocean('balance',256); %#ok<NASGU>
color_use_with_grey = [[0.7 0.7 0.7];color_use];
color_use_with_grey2 = [[0.7 0.7 0.7];flipud(color_use)];

% Subplot 1: C4 change
h1 = subplot(2,2,1);
axesm('MapProjection','robinson','Frame','off','Grid','off', ...
    'MapLatLimit',[-60 90],'MapLonLimit',[-180 180],'FontName','Helvetica');
A = nan(360,720);
k = isnan(c4_diff_us + S_change);
data = c4_diff_us;
data(k) = nan;
kknan_rec = find(isnan(data));
A(kland_05) = data;
kknans = isnan(A);
A_filt = A;
A_filt(~kknans) = imgaussfilt(A(~kknans),1);
A_smooth = A;
A_smooth(~kknans) = A_filt(~kknans);
A_smooth(kknans) = nan;
A_smooth(kland_05(kknan_rec)) = -10;
h = pcolorm(lat1,lon1,A_smooth);
set(h,'EdgeColor','none');
colormap(h1,color_use_with_grey);
caxis([-10 10]*0.6);
hold on;
plotm(coastlat,coastlon,'k-','LineWidth',1.5,'Color',[0 0 0]);
title('(a) \Delta C4','FontSize',16,'FontWeight','bold','FontName','Helvetica');
c1 = colorbar('southoutside');
c1.Label.String = 'Percentage (%)';
c1.FontSize = 12;
c1.FontName = 'Helvetica';

% Subplot 2: storage change
h2 = subplot(2,2,2);
axesm('MapProjection','robinson','Frame','off','Grid','off', ...
    'MapLatLimit',[-60 90],'MapLonLimit',[-180 180],'FontName','Helvetica');
A = nan(360,720);
k = isnan(c4_diff_us + S_change);
data = S_change;
data(k) = nan;
kknan_rec = find(isnan(data));
k = find(data < -0.8); data(k) = -0.8;
k = find(data > 0.8); data(k) = 0.8;
A(kland_05) = data;
kknans = isnan(A);
A_filt = A;
A_filt(~kknans) = imgaussfilt(A(~kknans),1);
A_smooth = A;
A_smooth(~kknans) = A_filt(~kknans);
A_smooth(kknans) = nan;
A_smooth(kland_05(kknan_rec)) = -1;
h = pcolorm(lat1,lon1,A_smooth);
set(h,'EdgeColor','none');
colormap(h2,color_use_with_grey2);
caxis([-1 1.001]*0.8);
hold on;
plotm(coastlat,coastlon,'k-','LineWidth',1.5,'Color',[0 0 0]);
title('(b) \Delta Z_{ P - \DeltaTWS}','FontSize',16,'FontWeight','bold','FontName','Helvetica');
c2 = colorbar('southoutside');
c2.Label.String = '[-]';
c2.FontSize = 12;
c2.FontName = 'Helvetica';
set(c2,'xtick',-1:0.2:1);

% Subplot 3: group means
h3 = subplot(2,2,3);
hh = bar(1:2,means,0.6,'FaceColor',[0.85 0.85 0.85], ...
    'EdgeColor','k','LineWidth',1.5);
hold on;
hh2 = errorbar(1:2,means,errors_lower,errors_upper,'k', ...
    'LineWidth',1.5,'CapSize',8,'Color',[0 0 0]);
set(gca,'xtick',[1 2],'xticklabel',{'\Delta C4<0','\Delta C4>0'},'FontName','Helvetica');
ylim([-0.06 0.02]); xlim([0.5 2.5]);
hh2.LineStyle = 'none';
hh.BarWidth = 0.5;
title('(c)','FontSize',16,'FontWeight','bold','FontName','Helvetica');
ylabel('\Delta Z_{ P - \DeltaTWS}','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',12);

% Subplot 4: relative response
h4 = subplot(2,2,4);
plot(co2,ad,'k-','LineWidth',2,'Color',[0 0 0]);
hold on;
yy = ones(size(co2));
plot(co2,yy,'k--','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
xx = [375 410 600 1130];
ads2 = spline(co2,ad,xx);
plot(xx,ads2,'ko','MarkerFaceColor',[0.9 0 0],'MarkerSize',8);
text(xx(1)+40,ads2(1),'2003','FontSize',10,'Color',[0 0 0.8],'FontName','Helvetica');
text(xx(2)+40,ads2(2),'2019','FontSize',10,'Color',[0 0 0.8],'FontName','Helvetica');
text(xx(3)-10,ads2(3)+0.085,{'2100','ssp245'},'FontSize',10,'Color',[0 0 0.8],'FontName','Helvetica');
text(xx(4)-150,ads2(4)+0.085,{'2100','ssp585'},'FontSize',10,'Color',[0 0 0.8],'FontName','Helvetica');
xlim([280 1200]); ylim([0.6 1.4]);
set(gca,'YTick',0.6:0.2:1.4,'FontName','Helvetica');
xlabel('Atmospheric CO_2 (ppm)','FontSize',12,'FontName','Helvetica');
ylabel('C4 Relative Advantage','FontSize',12,'FontName','Helvetica');
title('(d)','FontSize',16,'FontWeight','bold','FontName','Helvetica');
set(gca,'FontSize',12);

w1 = 0; hei = 0.475; wei = 0.575; hig = 0.52;
set(h1,'Position',[w1 hei wei hig]);
set(c1,'Position',[0.124 hei+0.035 0.35 0.02]);
set(h2,'Position',[w1+wei-0.125 hei wei hig]);
set(c2,'Position',[w1+wei-0.01 hei+0.035 0.35 0.02]);
set(h3,'Position',[w1+0.17 0.1 0.25 0.3]);
set(h4,'Position',[w1+wei+0.015 0.1 0.3 0.3]);
set(h1,'XColor','none','YColor','none','Color','none');
set(h2,'XColor','none','YColor','none','Color','none');

filename = fullfile(figure_dir,'Figure_3_C4_changes.png');
exportgraphics(fig,filename,'Resolution',300,'BackgroundColor','w');
fprintf('Created: %s\n',filename);

