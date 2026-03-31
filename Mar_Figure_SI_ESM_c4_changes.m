clear
load ../CPC_precp_daily/kland_05.mat
gd = 0.5;
a = -180+gd/2 :gd : 180-gd/2;
b = -90+gd/2 :gd : 90-gd/2;
[lon1,lat1] = meshgrid( a,b);
load coastlines.mat
gd = 1;
a = -180+gd/2 :gd : 180-gd/2;
b = -90+gd/2 :gd : 90-gd/2;
[lon2,lat2] = meshgrid( a,b);
% get a mask of significant C4 changes
load GIA_based_IA.mat
irr = nan(360,720); irr(kland_05) = IA;
load('../Modis_land_cover/LC_fraction.mat')
grass_fraction = LC_fraction;
name = 'C4_distribution_NUS_v2.2.nc'
a = flipud( ncread(name,'C4_area') );
% c4_diff = nanmean(a(:,:,11:end),3) - nanmean(a(:,:,1:10),3);
c4_diff = a(:,:,end) - a(:,:,1);
% c4_diff = a(:,:,15) - a(:,:,1);
k = find( grass_fraction< 50 | isnan(grass_fraction) | irr > 0.8); c4_diff(k) = nan;
for i = 1:19
    data = a(:,:,i);
    data(k) = nan;
    c4_time(i) = nanmean(nanmean(data,1),2);
end
c4_diff_us = c4_diff(kland_05);
a = prctile(c4_diff_us,99); b = prctile(c4_diff_us,1);
k = find(c4_diff_us>= a | c4_diff_us<=b); c4_diff_us(k) = nan;
co2_sens = readtable('../C4sensitivity/CO2_200-1400.csv');
co2 = co2_sens.CO2;
ad = co2_sens.Photosynthesis_C4./co2_sens.Photosynthesis_C3;
c4_his = matfile('../RAW_data/co2_C4/C4_fraction_historical.mat');
c4_245 = matfile('../RAW_data/co2_C4_245_585/C4_fraction_245.mat');
c4_585 = matfile('../RAW_data/co2_C4_245_585/C4_fraction_585.mat');
%%


models = {'ACCESS-CM2','ACCESS-ESM1-5',...
    'AWI-ESM','CanESM5-1','CESM2-WACCM',...
    'CMCC-CM2','IPSL-CM6A-LR','MPI-ESM1-2-HR','MPI-ESM1-2-LR'}
color_use = brewermap(256, 'BrBG');
color_use_with_grey =[[1 1 1]*0.75; color_use];
load coastlines.mat

%%
close all
fig = figure;
% Set figure dimensions and position
fig_width = 9.5; fig_height = 7.;
set(gcf, 'Position', [100 80 fig_width * 100 fig_height * 100])
diff_c4 = c4_245.C4_2020 - c4_his.C4_2000;
a = squeeze(nanmean(nanmean(diff_c4,1),2));
kkuse = find(~isnan(a));
ttest = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'}
for i = 1:9
    eval(['h',num2str(i),'=subplot(3,3,i)']);
    axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
        'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180], 'FontName', 'Helvetica');
    a = diff_c4(:,:,kkuse(i))';
    b = [a(:, 181:end), a(:, 1:180)];
    c = interp2(lon2, lat2, b, lon1, lat1);
    c4_median = c(kland_05);
    a = c4_median + c4_diff_us;
    k = find(isnan(a)); c4_median(k) = nan; c4_diff_us(k) = nan;
    % Plot the smoothed map
    %     plot_mask_025(0.75); hold on;
    A = nan(360, 720);
    A(kland_05) = c4_median; kloc = find(isnan(c4_median));
    kk = find(isnan(A)); A(kk) = nanmean(c4_median);
    A_smooth = imgaussfilt(A, 0.75);
    A_smooth(kk) = nan; A_smooth(kland_05(kloc)) = -100;
    h = pcolorm(lat1, lon1, A_smooth);
    set(h, 'EdgeColor', 'none');
    grid on;
    set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 12, 'LineWidth', 1.5);
    set(gca, 'Box', 'on', 'Layer', 'top');
    caxis([-1 1] * 10);
    colormap(color_use_with_grey);
    hold on;
    plotm(coastlat,coastlon, 'k-', 'LineWidth', 1.5);  % Thicker coastline
    title([ttest{i},' ',models{kkuse(i)}])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
    %     eval(['c',num2str(i),'=colorbar;'])
end
% Adjust subplot positions for consistency and aesthetics
vsp = 0.28;
w1 = -0.03; hei = 0.72; wei = 0.38; hight = 0.22; sp = -0.05;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h3, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h4, 'Position', [w1, hei, wei, hight]);
set(h5, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h6, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h7, 'Position', [w1, hei, wei, hight]);
set(h8, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h9, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
c1 = colorbar('southoutside')
set(c1, 'Position', [w1+(wei+sp)*0.75, hei-0.05, wei*1.5, 0.02]);
ylabel(c1,'C4 change (%)');

filename='.\Figures\Figure_SI_C4_reduction_historical.png'
exportgraphics(fig, filename, 'Resolution', 300);

%%
close all
fig = figure;
% Set figure dimensions and position
fig_width = 9.5; fig_height = 7.;
set(gcf, 'Position', [100 80 fig_width * 100 fig_height * 100])
diff_c4 = c4_245.C4_2100 - c4_245.C4_2020;
a = squeeze(nanmean(nanmean(diff_c4,1),2));
kkuse = find(~isnan(a));
ttest = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'}
for i = 1:9
    eval(['h',num2str(i),'=subplot(3,3,i)']);
    axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
        'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180], 'FontName', 'Helvetica');
    a = diff_c4(:,:,kkuse(i))';
    b = [a(:, 181:end), a(:, 1:180)];
    c = interp2(lon2, lat2, b, lon1, lat1);
    c4_median = c(kland_05);
    a = c4_median + c4_diff_us;
    k = find(isnan(a)); c4_median(k) = nan; c4_diff_us(k) = nan;
    % Plot the smoothed map
    %     plot_mask_025(0.75); hold on;
    A = nan(360, 720);
    A(kland_05) = c4_median; kloc = find(isnan(c4_median));
    kk = find(isnan(A)); A(kk) = nanmean(c4_median);
    A_smooth = imgaussfilt(A, 0.75);
    A_smooth(kk) = nan; A_smooth(kland_05(kloc)) = -100;
    h = pcolorm(lat1, lon1, A_smooth);
    set(h, 'EdgeColor', 'none');
    grid on;
    set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 12, 'LineWidth', 1.5);
    set(gca, 'Box', 'on', 'Layer', 'top');
    caxis([-1 1] * 10);
    colormap(color_use_with_grey);
    hold on;
    plotm(coastlat,coastlon, 'k-', 'LineWidth', 1.5);  % Thicker coastline
    title([ttest{i},' ',models{kkuse(i)}])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
    %     eval(['c',num2str(i),'=colorbar;'])
end
% Adjust subplot positions for consistency and aesthetics
vsp = 0.28;
w1 = -0.03; hei = 0.72; wei = 0.38; hight = 0.22; sp = -0.05;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h3, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h4, 'Position', [w1, hei, wei, hight]);
set(h5, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h6, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h7, 'Position', [w1, hei, wei, hight]);
set(h8, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h9, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
c1 = colorbar('southoutside')
set(c1, 'Position', [w1+(wei+sp)*0.75, hei-0.05, wei*1.5, 0.02]);
ylabel(c1,'C4 change (%)');

filename='.\Figures\Figure_SI_C4_reduction_245.png'
exportgraphics(fig, filename, 'Resolution', 300);

%%
close all
fig = figure;
% Set figure dimensions and position
fig_width = 9.5; fig_height = 7.;
set(gcf, 'Position', [100 80 fig_width * 100 fig_height * 100])
diff_c4 = c4_585.C4_2100 - c4_585.C4_2020;
a = squeeze(nanmean(nanmean(diff_c4,1),2));
kkuse = find(~isnan(a));
ttest = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'}
for i = 1:9
    eval(['h',num2str(i),'=subplot(3,3,i)']);
    axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
        'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180], 'FontName', 'Helvetica');
    a = diff_c4(:,:,kkuse(i))';
    b = [a(:, 181:end), a(:, 1:180)];
    c = interp2(lon2, lat2, b, lon1, lat1);
    c4_median = c(kland_05);
    a = c4_median + c4_diff_us;
    k = find(isnan(a)); c4_median(k) = nan; c4_diff_us(k) = nan;
    % Plot the smoothed map
    %     plot_mask_025(0.75); hold on;
    A = nan(360, 720);
    A(kland_05) = c4_median; kloc = find(isnan(c4_median));
    kk = find(isnan(A)); A(kk) = nanmean(c4_median);
    A_smooth = imgaussfilt(A, 0.75);
    A_smooth(kk) = nan; A_smooth(kland_05(kloc)) = -100;
    h = pcolorm(lat1, lon1, A_smooth);
    set(h, 'EdgeColor', 'none');
    grid on;
    set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 12, 'LineWidth', 1.5);
    set(gca, 'Box', 'on', 'Layer', 'top');
    caxis([-1 1] * 10);
    colormap(color_use_with_grey);
    hold on;
    plotm(coastlat,coastlon, 'k-', 'LineWidth', 1.5);  % Thicker coastline
    title([ttest{i},' ',models{kkuse(i)}])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
    %     eval(['c',num2str(i),'=colorbar;'])
end
% Adjust subplot positions for consistency and aesthetics
vsp = 0.28;
w1 = -0.03; hei = 0.72; wei = 0.38; hight = 0.22; sp = -0.05;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h3, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h4, 'Position', [w1, hei, wei, hight]);
set(h5, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h6, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
hei = hei - vsp;
set(h7, 'Position', [w1, hei, wei, hight]);
set(h8, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h9, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);
c1 = colorbar('southoutside')
set(c1, 'Position', [w1+(wei+sp)*0.75, hei-0.05, wei*1.5, 0.02]);
ylabel(c1,'C4 change (%)');

filename='.\Figures\Figure_SI_C4_reduction_585.png'
exportgraphics(fig, filename, 'Resolution', 300);

