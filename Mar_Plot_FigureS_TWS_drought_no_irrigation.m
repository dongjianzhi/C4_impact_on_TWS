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


name = 'C4_distribution_NUS_v2.2.nc'
a = flipud( ncread(name,'C4_area') );
c4_diff = nanmean(a(:,:,11:end),3) - nanmean(a(:,:,1:10),3);
c4_org = nanmean(a,3);

k = find( grass_fraction< 50 | isnan(grass_fraction) | irr > 0.); c4_org(k) = nan;


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
    
    
    if J == 1; spi_thres = 0; vpd_thres = 1; end
    if J == 2; spi_thres = 0; vpd_thres = 0; end
    if J == 3; spi_thres = 0; vpd_thres = 1.5; end
    if J == 4; spi_thres = -1; vpd_thres = 0; end
    if J == 4; spi_thres = -1; vpd_thres = 1.5; end


    [et_change] = cal_water_storage_change_warm(tws_detrend, SPI, time, lats, vpd_mon, vpdz,spi_thres, vpd_thres,20);
     dsd_all(:,J) = et_change;
end



lats = lat1(kland_05);
W = cos( lats/180*pi);

k = find(isnan(c4_diff_us) );
dsd_no_cp = (dsd_all(k,:).*W(k) );
k = find(dsd_no_cp < prctile(dsd_no_cp,5) ); dsd_no_cp(k) = nan;
k = find(dsd_no_cp > prctile(dsd_no_cp,95) ); dsd_no_cp(k) = nan;

k = find(~isnan(c4_diff_us) );
dsd_cp = (dsd_all(k,:).*W(k) );
k = find(dsd_cp < prctile(dsd_cp,5) ); dsd_cp(k) = nan;
k = find(dsd_cp > prctile(dsd_cp,95) ); dsd_cp(k) = nan;


for i = 1:500
    k1 = randsample(1:length(dsd_cp),length(dsd_cp),true);
    k2 = randsample(1:length(dsd_no_cp),length(dsd_no_cp),true);

    dmeans_cp(i,:) = [nanmean(dsd_cp(k1,:))];
    dmeans_no_cp(i,:) =  [nanmean(dsd_no_cp(k2,:))];

    dsums_cp(i,:) = [nansum(dsd_cp(k1,:))];
    dsums_no_cp(i,:) =  [nansum(dsd_no_cp(k2,:))];
end

dmeans = [nanmean(dmeans_cp)' nanmean(dmeans_no_cp)'];
dstds = [range(dmeans_cp)' range(dmeans_no_cp)']/2;

dsums = [nanmean(dsums_cp)' nanmean(dsums_no_cp)'];
dsum_stds = [range(dsums_cp)' range(dsums_no_cp)']/2;

% save('total_tws_change_2003_2019','dsums','dsum_stds')

%%
% Create the bar plot
rd = viridis;

close all;
fig = figure;
width = 7;
height = 3.;
set(gcf, 'Position', [100 80 width*100 height*100]);
load coastlines.mat;

hold on;
kkuse = [2 4 1 3 5];
group_shift = 0.15; % Offset for grouped bars
x = [1:5]; 
bar_width = 0.8; % Slightly narrower bars for elegance

% Bar plot with improved colors
hh = bar(x, dmeans(kkuse, :), 'grouped', 'BarWidth', bar_width);
hh(1).FaceColor = [0.34, 0.63, 0.85]; % Soft blue (grass/crop)
hh(2).FaceColor = [0.85, 0.47, 0.31]; % Warm coral (others)

% Add error bars with refined style
errorbar(x - group_shift, dmeans(kkuse, 1), dstds(kkuse, 1), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
errorbar(x + group_shift, dmeans(kkuse, 2), dstds(kkuse, 2), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);

% Aesthetic improvements
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5]); % Subtle grid
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1); % Clean font and axes
ylim([-1.6 0]); % Adjusted y-axis limits
xlim([0.5 5.5]);
xticks(x);
xticklabels({'Z_P<0', 'Z_P<-1', 'Z_P<0 & VPD>1', 'Z_P<0 & VPD>1.5', 'Z_P<-1 & VPD>1.5'});
set(gca, 'XTickLabelRotation', 0); % Rotate labels for readability
ylabel('cm', 'FontSize', 12, 'FontName', 'Arial');

% Add legend with refined style
legend({'Grass/Crop', 'Others'}, 'Location', 'SouthWest', 'FontSize', 11, ...
    'Box', 'off', 'FontName', 'Arial');


filename='./Figures/FigureS_drought_impacts_TWS_no_irrigation.png'; 
exportgraphics(fig, filename, 'Resolution', 300);



