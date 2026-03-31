clear
A = matfile('monthly_data.mat'); time = A.um;
load('GRACE_CSR.mat')
load ../CPC_precp_daily/kland_05.mat
Precp = A.Precp;
vpd_mon = A.vpd_mon;

gd = 0.5;
a = -180+gd/2 :gd : 180-gd/2;
b = -90+gd/2 :gd : 90-gd/2;
[lon1,lat1] = meshgrid( a,b);

latVec = b; lonVec = a;
load coastlines.mat


%% get a mask of significant C4 changes
load GIA_based_IA.mat
irr = nan(360,720); irr(kland_05) = IA;

load('../Modis_land_cover/LC_fraction.mat')
grass_fraction = LC_fraction;

A = matfile("../Modis_land_cover/LC_fraction_forest.mat")
forest_fraction = A.LC_fraction;

A = matfile("../Modis_land_cover/LC_fraction_ice_only.mat")
ice_fraction = A.LC_fraction;



name = 'C4_distribution_NUS_v2.2.nc'
a = flipud( ncread(name,'C4_area') );
c4_diff = nanmean(a(:,:,11:end),3) - nanmean(a(:,:,1:10),3);
c4_org = nanmean(a,3);

k = find( grass_fraction< 50 | isnan(grass_fraction) | irr > 0.01); c4_org(k) = nan;


c4_diff_us = c4_diff(kland_05);  c4_org_us = c4_org(kland_05);
a = prctile(c4_diff_us,95); b = prctile(c4_diff_us,5);
k = find(c4_diff_us>= a | c4_diff_us<=b); c4_diff_us(k) = nan;

% Step 4: Map Pixel Coordinates to Geographic Coordinates
binaryImage = imbinarize(c4_org);

% Step 3: Trace Boundaries
boundaries = bwboundaries(binaryImage);

% Display boundaries on the binary image
clear bd_corr polygons
N = 1;
for k = 1:length(boundaries)
    boundary = boundaries{k};
    bd_lon = lonVec(boundary(:,2)); bd_lat = latVec(boundary(:,1) );
    if length(bd_lat)>10
        
        a = polyshape(bd_lon,bd_lat);
        polygons{N} = a;
        N = N + 1;
    end
    
end

%% CRU data

B = matfile( 'Data_8day_glb_0421.mat');
pre_8d_cru = B.pre_8d_cru;
vpd_8d_cru = B.vpd_8d_cru;
time_8d = B.time_8d;
H_flux = B.H_flux;
LE_flux = B.LE_flux;
A = H_flux + LE_flux;

% into monthly
time_all = [time; time(end) + 31];
pre_mon = nan(size(csr_tws) );
vpd_mon = pre_mon;
A_mon  = pre_mon;  LE_mon = A_mon;
for i  = 1:length(time_all) - 1
    k = find( time_8d>= time_all(i) & time_8d<= time_all(i+1) );
    pre_mon(:,i) = nansum( pre_8d_cru(:,k),2)/10; % into cm/month
    vpd_mon(:,i) = nanmean( vpd_8d_cru(:,k),2);
    A_mon(:,i) = nanmean( A(:,k),2);
    LE_mon(:,i) = nanmean( LE_flux(:,k),2);
end

%% Calculate TWS
tws = nanmean( cat(3,csr_tws ,JPL_tws , gsfc_tws),3) ;
dS = [tws(:,2:end)-tws(:,1:end-1),nan(length(kland_05),1)];
vm = [vpd_mon(:,2:end)+vpd_mon(:,1:end-1),nan(length(kland_05),1)]/2;
[~, vpdz] = cal_mon_zcore_linear_by_month(time, vm);

[~, pz] = cal_mon_zcore_linear_by_month(time, Precp);
SPI = pz;

%% Only keep the period of 2003 - 2019
k = find( time>= datenum(2003,1,1) & time<= datenum(2019,12,31));
vpd_mon = vpd_mon(:,k); SPI = SPI(:,k); tws = tws(:,k);
time = time(k);
%% Calculate Z score
% kick out the edges
load kland_025.mat
A = nan(720,1440);  A(kland_025) = 1;
B = A(1:2:end,:) + A(2:2:end,:); C = B(:,1:2:end) + B(:,2:2:end);
kclean = find(isnan(C) );

dsd_all = nan( length(kland_05), 5);

lats = lat1(kland_05);
% tws(:,265:end) = nan;
tws_detrend = remove_linear_trend(time,tws);
for J = 1:5
    
    if J == 1; spi_thres = 0; vpd_thres = 0; end
    if J == 2; spi_thres = -1; vpd_thres = 0; end
    if J == 3; spi_thres = 0; vpd_thres = 1; end
    if J == 4; spi_thres = -1; vpd_thres = 1; end
    
%     if J == 2; spi_thres = 0; vpd_thres = 1; end
% %     if J == 3; spi_thres = 0; vpd_thres = 1.5; end
%     if J == 4; spi_thres = -1; vpd_thres = 0; end
%     if J == 4; spi_thres = -1; vpd_thres = 1; end


    [et_change] = cal_water_storage_change_warm(tws_detrend, SPI, time, lats, vpd_mon, vpdz,spi_thres, vpd_thres,20);
     dsd_all(:,J) = et_change;
end
% 
% kirr = irr(kland_05); k = find(kirr>0.1);
% dsd_all(k,:) = nan;
%%

A = forest_fraction(kland_05);
k_forests = A>50;

A = ice_fraction(kland_05);
k_ice = A>0; dsd_all(k_ice,:) = nan;


lats = lat1(kland_05);
W = cos( lats/180*pi);

k = find(isnan(c4_diff_us)  & ~ k_forests );
dsd_no_cp_fre = (dsd_all(k,:).*W(k) );
k = find(dsd_no_cp_fre < prctile(dsd_no_cp_fre,5) ); dsd_no_cp_fre(k) = nan;
k = find(dsd_no_cp_fre > prctile(dsd_no_cp_fre,95) ); dsd_no_cp_fre(k) = nan;

dsd_fre = (dsd_all(k_forests,:).*W(k_forests) );


k = find(~isnan(c4_diff_us) );
dsd_cp = (dsd_all(k,:).*W(k) );
k = find(dsd_cp < prctile(dsd_cp,5) ); dsd_cp(k) = nan;
k = find(dsd_cp > prctile(dsd_cp,95) ); dsd_cp(k) = nan;



%%
% land area
a = sum(lats.*W);
k = find(~isnan(c4_diff_us) );
b = sum(lats(k).*W(k));
c4_area = b./a

for i = 1:500
    k1 = randsample(1:length(dsd_cp),length(dsd_cp),true);
    k2 = randsample(1:length(dsd_no_cp_fre),length(dsd_no_cp_fre),true);
    k3 = randsample(1:length(dsd_fre),length(dsd_fre),true);

    dmeans_cp(i,:) = [nanmean(dsd_cp(k1,:))];
    dmeans_no_cp_fre(i,:) =  [nanmean(dsd_no_cp_fre(k2,:))];
    dmeans_fre(i,:) =  [nanmean(dsd_fre(k3,:))];

    dsums_cp(i,:) = [nansum(dsd_cp(k1,:))];
    dsums_no_cp_fre(i,:) =  [nansum(dsd_no_cp_fre(k2,:))];
    dsums_fre(i,:) =  [nansum(dsd_fre(k3,:))];

end

dmeans = [nanmean(dmeans_cp)' nanmean(dmeans_fre)' nanmean(dmeans_no_cp_fre)'];
dstds = [range(dmeans_cp)' range(dmeans_fre)' range(dmeans_no_cp_fre)']/2;

dsums_cp1 = dsums_cp./( dsums_cp + dsums_no_cp_fre +dsums_fre);
dsums_no_cp_fre1 = dsums_no_cp_fre./( dsums_cp + dsums_no_cp_fre +dsums_fre);
dsums_fre1 = dsums_fre./( dsums_cp + dsums_no_cp_fre +dsums_fre);

dsums = -[nanmean(dsums_cp1)' nanmean(dsums_fre1)' nanmean(dsums_no_cp_fre1)'];
dsum_stds = [range(dsums_cp1)' range(dsums_fre1)' range(dsums_no_cp_fre1)']/2;

% save('total_tws_change_2003_2019','dsums','dsum_stds')

%%
% Create the bar plot
rd = viridis;

close all;
fig = figure;
width = 10.5;
height = 4.5;
set(gcf, 'Position', [100 80 width*100 height*100]);
load coastlines.mat;

rd = brewermap(256, 'OrRd'); 
tt = 2;
tt = 1.8;
kk = round(linspace(1, 256^tt, 220).^(1/tt));

% Subplot 1: Map
h1 = subplot(3, 1, 1);
axesm('robinson', 'MapLatLimit', [-60 90], 'Grid', 'off', 'Frame', 'off');
set(gca, 'FontName', 'Arial', 'FontSize', 12);

% 1. Paint ALL land grey first
% This uses the built-in 'landareas' shapefile
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); 

hold on;

% 2. Prepare your data
A = nan(360, 720); 
A(kland_05) = dsd_all(:, 1); 
A(kclean) = nan;
A = double(A);

% 3. Plot the data on top of the grey land
h = pcolorm(lat1, lon1, A); 
% Ensure values higher than 0 don't get colored by the drought colormap
% or adjust caxis to ensure 0 is the limit
caxis([-5 0]); 
set(h, 'EdgeColor', 'none');

% 4. Colormap and Coastlines
colormap(h1, flipud(rd(kk, :)));
plotm(coastlat, coastlon, 'Color', [0.25 0.25 0.25], 'LineWidth', 1);

% Add hatched polygons manually
for i = 1:length(polygons)
    a = polygons{i}.Vertices(:, 1); % Longitude
    b = polygons{i}.Vertices(:, 2); % Latitude
    a1 = nanmin(a):2:nanmax(a);
    b1 = nanmin(b):2:nanmax(b);
    [lon, lat] = meshgrid(a1, b1);
    IN = inpolygon(lon, lat, a, b);
    lon_use = lon(IN);
    lat_use = lat(IN);
    scatterm(lat_use, lon_use, 1.5, 'k', 'filled');
end

tightmap;
set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
set(h1, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');

c1 = colorbar('southoutside');
xlabel(c1, 'cm', 'FontSize', 12, 'FontName', 'Arial');

title('(a) E(TWSd | Z_p<0)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
c1.TickLabels{end} = '\geq 0';
set(h1, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');


% Subplot 2: Bar Plot
h2 = subplot(3, 1, 2);
hold on;
kkuse = [1:4];
group_shift = 0.225; % Offset for grouped bars
x = [1:4]; 
bar_width = 0.8; % Slightly narrower bars for elegance

% Bar plot with improved colors
hh = bar(x, dmeans(kkuse, :), 'grouped', 'BarWidth', bar_width);
hh(1).FaceColor = [0.34, 0.63, 0.85]; % Soft blue (grass/crop)
hh(2).FaceColor = [0.45, 0.75, 0.45]; % Warm coral (others)
hh(3).FaceColor = [0.85, 0.47, 0.31]; % Soft green (e.g., forest/other vegetation)

% Add error bars with refined style
errorbar(x - group_shift, dmeans(kkuse, 1), dstds(kkuse, 1), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
errorbar(x , dmeans(kkuse, 2), dstds(kkuse, 2), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
errorbar(x + group_shift, dmeans(kkuse, 3), dstds(kkuse, 3), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);

% Aesthetic improvements
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5]); % Subtle grid
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1); % Clean font and axes
ylim([-3 0]); % Adjusted y-axis limits
xlim([0.5 4.5]);
xticks(x);
xticklabels({'Z_P<0','Z_P<-1', 'Z_P<0, VPD>1kpa', 'Z_P<0, VPD>1.5kpa'});
set(gca, 'XTickLabelRotation', 0); % Rotate labels for readability
ylabel('cm', 'FontSize', 12, 'FontName', 'Arial');

%    if J == 1; spi_thres = 0; vpd_thres = 1; end
%     if J == 2; spi_thres = 0; vpd_thres = 0; end
%     if J == 3; spi_thres = 0; vpd_thres = 1.5; end
%     if J == 4; spi_thres = -1; vpd_thres = 0; end
%     if J == 4; spi_thres = -1; vpd_thres = 1.5; end


% Title and final touches
title('(b) Mean TWSd during droughts', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
set(gca, 'YTick', -3:1:0);
box on;

% Add legend with refined style
legend({'Grass/Crop', 'Forest','Others'}, 'Location', 'SouthWest', 'FontSize', 11, ...
    'Box', 'off', 'FontName', 'Arial');



% Subplot 3: Bar Plot
h3 = subplot(3, 1, 3);
hold on;
kkuse = [1:4];
group_shift = 0.225; % Offset for grouped bars
x = [1:4]; 
bar_width = 0.8; % Slightly narrower bars for elegance

% Bar plot with improved colors
% dsums = -dsums./sum(dsums,2);
hh = bar(x, dsums(kkuse, :), 'grouped', 'BarWidth', bar_width);
hh(1).FaceColor = [0.34, 0.63, 0.85]; % Soft blue (grass/crop)
hh(2).FaceColor = [0.45, 0.75, 0.45]; % Warm coral (others)
hh(3).FaceColor = [0.85, 0.47, 0.31]; % Soft green (e.g., forest/other vegetation)

% Add error bars with refined style
errorbar(x - group_shift, dsums(kkuse, 1), dsum_stds(kkuse, 1), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
errorbar(x , dsums(kkuse, 2), dsum_stds(kkuse, 2), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
errorbar(x + group_shift, dsums(kkuse, 3), dsum_stds(kkuse, 3), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);

% Aesthetic improvements
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5]); % Subtle grid
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1); % Clean font and axes
% ylim([-3 0]); % Adjusted y-axis limits
xlim([0.5 4.5]);
xticks(x);
xticklabels({'Z_P<0','Z_P<-1', 'Z_P<0, VPD>1kpa', 'Z_P<0, VPD>1.5kpa'});
set(gca, 'XTickLabelRotation', 0); % Rotate labels for readability
ylabel('Fraction', 'FontSize', 12, 'FontName', 'Arial');


% Title and final touches
title('(c) Contribution to global TWSd depletion', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
% set(gca, 'YTick', -3:1:0);
box on;




% Adjust subplot positions
set(h1, 'Position', [-0.22 0.15 0.9 0.62]);
set(c1, 'Position', [0.03 0.15 0.45 0.015]);

set(h2, 'Position', [0.57 0.55 0.42 0.26]);

set(h3, 'Position', [0.57 0.15 0.42 0.26]);


filename='./Figures/Figure1_drought_impacts_TWS_upd.png'; 
exportgraphics(fig, filename, 'Resolution', 300);