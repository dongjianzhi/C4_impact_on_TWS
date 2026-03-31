clear; close all
A = matfile('monthly_data.mat'); time = A.um;
load('GRACE_CSR.mat')
load ../CPC_precp_daily/kland_05.mat
Precp = A.Precp;
vpd_mon = A.vpd_mon;
T_mon = A.T_mon;
Rs_mon = A.RS_mon;
gd = 0.5;
a = -180+gd/2 :gd : 180-gd/2;
b = -90+gd/2 :gd : 90-gd/2;
[lon1,lat1] = meshgrid( a,b);
load coastlines.mat
SIF = A.csfi_mon;
A = matfile('../MERRA2_VPD_data/MERRA2_monthly_pet.mat');
um = A.um;
data = A.pet_mon;
[~,ib,ic] = intersect(um,time);
pet = data(:,ib);
% CRU data
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

% Calculate TWS
tws = nanmean( cat(3,csr_tws ,JPL_tws , gsfc_tws),3) ;
% get a mask of significant C4 changes

%%
load GIA_based_IA.mat
irr = nan(360,720); irr(kland_05) = IA;
load('../Modis_land_cover/LC_fraction.mat')
grass_fraction = LC_fraction;
name = 'C4_distribution_NUS_v2.2.nc'
a = flipud( ncread(name,'C4_area') );
C4 = nanmean(a,3);
b = nanmean(a(:,:,10:end),3) - nanmean(a(:,:,1:9),3) ;
k = find( grass_fraction< 50 | isnan(grass_fraction) | irr>1  ); C4(k) = nan;
C4_use = C4(kland_05);
b(k) = nan;
C4_change = b(kland_05);
k = find(C4_change>= prctile(C4_change,95) | C4_change<= prctile(C4_change,5));
C4_change(k) = nan;  c4_diff_us = C4_change;

%%
% Calculate Z score
dS = [nan(length(kland_05),1) tws(:,2:end)-tws(:,1:end-1) ];
P = [nan(length(kland_05),1) Precp(:,2:end)+Precp(:,1:end-1) ]/20;
T = [nan(length(kland_05),1) T_mon(:,2:end)+T_mon(:,1:end-1) ]/2;
R = [nan(length(kland_05),1) Rs_mon(:,2:end)+Rs_mon(:,1:end-1) ]/2;
V = [nan(length(kland_05),1) vpd_mon(:,2:end)+vpd_mon(:,1:end-1) ]/2;

et = (P - dS);

% [~, tws_z] = cal_mon_zcore_linear_by_month(time, tws);
[~, dSz] = cal_mon_zcore_linear_by_month(time, et);
[~, Az] = cal_mon_zcore_linear_by_month(time, R);
[~, Tz] = cal_mon_zcore_linear_by_month(time, T);
[~, SPI] = cal_mon_zcore_linear_by_month(time, P);

% Calculate conditional TWS
% A = matfile('SPI3_data.mat'); SPI = A.SPI3;
lats = lat1(kland_05);
[vpdz, ~] = cal_mon_zcore_linear_by_month(time, V);
%  calculate the Z score of different factors
lats = lat1(kland_05);
vpd_thres = 1; spi_thres = 0; vpd_thres2 = 20;
% vpd_thres = 1.; spi_thres = 0; vpd_thres2 = 30;
% [ds_change,ds_spi,ds_vpd] = cal_water_storage_change(dSz, SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres);
% Tz = T_mon; Az = A_mon; vpdz = vpd_mon;
% [et_change,~,~] = cal_water_storage_change(etz, SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2003,1,1) & time<=datenum(2009,12,31));
[S_change1] = cal_water_storage_change_warm(dSz(:,kk), ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2010,1,1) & time<=datenum(2019,12,31));
[S_change2] = cal_water_storage_change_warm(dSz(:,kk), ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
S_change = S_change2 - S_change1;

kk = find( time>= datenum(2003,1,1) & time<=datenum(2009,12,31));
[spi_change1] = cal_water_storage_change_warm(SPI(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2010,1,1) & time<=datenum(2019,12,31));
[spi_change2] = cal_water_storage_change_warm(SPI(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
spi_change = spi_change2 - spi_change1;

kk = find( time>= datenum(2003,1,1) & time<=datenum(2009,12,31));
[A_change1] = cal_water_storage_change_warm(Az(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2010,1,1) & time<=datenum(2019,12,31));
[A_change2] = cal_water_storage_change_warm(Az(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
A_change = A_change2 - A_change1;

kk = find( time>= datenum(2003,1,1) & time<=datenum(2009,12,31));
[v_change1] = cal_water_storage_change_warm(vpdz(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2010,1,1) & time<=datenum(2019,12,31));
[v_change2] = cal_water_storage_change_warm(vpdz(:,kk),  ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
v_change = v_change2 - v_change1;

kk = find( time>= datenum(2003,1,1) & time<=datenum(2009,12,31));
[T_change1] = cal_water_storage_change_warm(Tz(:,kk),   ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
kk = find( time>= datenum(2010,1,1) & time<=datenum(2019,12,31));
[T_change2] = cal_water_storage_change_warm(Tz(:,kk),   ...
    SPI(:,kk), time(kk), lats, vpd_mon(:,kk), vpdz(:,kk),spi_thres,vpd_thres,vpd_thres2);
T_change = T_change2 - T_change1;

%% Bootstrapping the group means
k = ~isnan( c4_diff_us + spi_change + S_change + A_change + v_change  ); kkuse = kland_05(k); kk_ind = k;
% k = ~isnan( C4_change + spi_change + S_change + A_change + v_change  ); kkuse = kland_05(k); kk_ind = k;
y = S_change(k);  x2 = spi_change(k); x3 = A_change(k); x4 = v_change(k); x5 = T_change(k);
% train random forest
x6 = c4_diff_us(k);

X = [x6,x2,x3,x4,x5];

% Group y based on X(:,1) < 0 and X(:,1) > 0
group = X(:,1) < 0; % Logical index
y_neg = y(group);   % y values where X(:,1) < 0
y_pos = y(~group);  % y values where X(:,1) > 0

% Bootstrap settings
n_boot = 1000; % Number of bootstrap samples
ci_level = 0.95; % 95% confidence interval
alpha = (1 - ci_level) / 2; % For 2.5th and 97.5th percentiles

% Bootstrap for y_neg
boot_means_neg = zeros(n_boot, 1);
for i = 1:n_boot
    boot_sample = randsample(y_neg, length(y_neg), true); % Resample with replacement
    boot_means_neg(i) = mean(boot_sample);
end
ci_neg = prctile(boot_means_neg, [100*alpha, 100*(1-alpha)]); % 95% CI

% Bootstrap for y_pos
boot_means_pos = zeros(n_boot, 1);
for i = 1:n_boot
    boot_sample = randsample(y_pos, length(y_pos), true); % Resample with replacement
    boot_means_pos(i) = mean(boot_sample);
end
ci_pos = prctile(boot_means_pos, [100*alpha, 100*(1-alpha)]); % 95% CI

% Compute observed means and error bar lengths
means = [mean(y_neg), mean(y_pos)];
errors_lower = means - [ci_neg(1), ci_pos(1)]; % Distance to lower bound
errors_upper = [ci_neg(2), ci_pos(2)] - means; % Distance to upper bound

%% get a mask of significant C4 changes

co2_sens = readtable('../C4sensitivity/CO2_200-1400.csv');
co2 = co2_sens.CO2;
ad = co2_sens.Photosynthesis_C4./co2_sens.Photosynthesis_C3;
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


% color_use = brewermap(128, 'BrBG');
% load coastlines.mat
dd = slanCM('BrBG');
color_use = dd(16:1:240,:)
c4_245 = matfile('../RAW_data/co2_C4_245_585/C4_fraction_245.mat');
c4_585 = matfile('../RAW_data/co2_C4_245_585/C4_fraction_585.mat');
c4_his = matfile('../RAW_data/co2_C4/C4_fraction_historical.mat');

%%
close all;
fig = figure;
width = 9.5; % inches, suitable for Nature's column width (~90 mm or ~3.5 inches per column)
height = 6.5; % inches
set(gcf, 'Position', [100 80 width*100 height*100], 'Color', 'w', 'PaperUnits', 'inches', 'PaperSize', [width height]);

color_use = cmocean('tarn', 256); % High-resolution diverging colormap.\
color_use2 = cmocean('balance', 256); % High-resolution diverging colormap

color_use_with_grey = [[0.7 0.7 0.7]; (color_use)]; % Light gray for NaN regions
color_use_with_grey2 = [[0.7 0.7 0.7]; flipud(color_use)]; % Light gray for NaN regions


% Subplot 1: C4 Change
h1 = subplot(2, 2, 1);
axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
      'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180], 'FontName', 'Helvetica');
A = nan(360, 720);
k = isnan(c4_diff_us + S_change); 
data = c4_diff_us; data(k) = nan; kknan_rec = find(isnan(data));
% k = find(data<-6); data(k) = -6;
A(kland_05) = data; kknans = isnan(A); 
% Replace NaN with mean for filtering (manual NaN handling)
A_filt = A;
A_filt(~kknans) = imgaussfilt(A(~kknans), 1); % Filter non-NaN values
A_smooth = A; % Restore NaN structure
A_smooth(~kknans) = A_filt(~kknans);
A_smooth(kknans) = nan; A_smooth(kland_05(kknan_rec)) = -10;
h = pcolorm(lat1, lon1, A_smooth);
set(h, 'EdgeColor', 'none');
colormap(h1, color_use_with_grey);
caxis([-10 10]*0.6);
hold on;
plotm(coastlat, coastlon, 'k-', 'LineWidth', 1.5, 'Color', [0 0 0]);
title('(a) \Delta C4', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
c1 = colorbar('southoutside'); 
c1.Label.String = 'Percentage (%)'; 
c1.FontSize = 12; c1.FontName = 'Helvetica';
% set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');

% Subplot 2: Δ Z_L
h2 = subplot(2, 2, 2);
axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
      'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180], 'FontName', 'Helvetica');
A = nan(360, 720);
k = isnan(c4_diff_us + S_change); 
data = S_change; data(k) = nan; kknan_rec = find(isnan(data)); 
k = find(data<-0.8); data(k) = -0.8; k = find(data>0.8); data(k) = 0.8;

A(kland_05) = data; kknans = isnan(A); 
% Replace NaN with mean for filtering (manual NaN handling)
A_filt = A;
A_filt(~kknans) = imgaussfilt(A(~kknans), 1); % Filter non-NaN values
A_smooth = A; % Restore NaN structure
A_smooth(~kknans) = A_filt(~kknans);
A_smooth(kknans) = nan; A_smooth(kland_05(kknan_rec)) = -1;
h = pcolorm(lat1, lon1, A_smooth);
set(h, 'EdgeColor', 'none');
colormap(h2, color_use_with_grey2);
caxis([-1 1.001]*0.8);
hold on;
plotm(coastlat, coastlon, 'k-', 'LineWidth', 1.5, 'Color', [0 0 0]);
title('(b) \Delta Z_{ P - \DeltaTWS}', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
c2 = colorbar('southoutside'); 
c2.Label.String = '[-]'; 
c2.FontSize = 12; c2.FontName = 'Helvetica';
set(c2,'xtick',[-1:0.2:1])
% set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');


% Subplot 3: Bar Plot
h3 = subplot(2, 2, 3);
hh = bar(1:2, means, 0.6, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on;
hh2 = errorbar(1:2, means, errors_lower, errors_upper, 'k', 'LineWidth', 1.5, 'CapSize', 8, 'Color', [0 0 0]);
set(gca, 'xtick', [1 2], 'xticklabel', {'\Delta C4<0', '\Delta C4>0'}, 'FontName', 'Helvetica');
ylim([-0.06 0.02]); xlim([0.5 2.5]);
hh2.LineStyle = 'none'; hh.BarWidth = 0.5;
title('(c)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
ylabel('\Delta Z_{ P - \DeltaTWS}', 'FontSize', 12, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12);

% Subplot 4: Line Plot
h4 = subplot(2, 2, 4);
plot(co2, ad, 'k-', 'LineWidth', 2, 'Color', [0 0 0]);
hold on;
yy = ones(size(co2));
plot(co2, yy, 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xx = [375 410 600 1130];
ads2 = spline(co2, ad, xx);
plot(xx, ads2, 'ko', 'MarkerFaceColor', [0.9 0 0], 'MarkerSize', 8);
text(xx(1)+40, ads2(1), '2003', 'FontSize', 10, 'Color', [0 0 0.8], 'FontName', 'Helvetica');
text(xx(2)+40, ads2(2), '2019', 'FontSize', 10, 'Color', [0 0 0.8], 'FontName', 'Helvetica');
text(xx(3)-10, ads2(3)+0.085, {'2100', 'ssp245'}, 'FontSize', 10, 'Color', [0 0 0.8], 'FontName', 'Helvetica');
text(xx(4)-150, ads2(4)+0.085, {'2100', 'ssp585'}, 'FontSize', 10, 'Color', [0 0 0.8], 'FontName', 'Helvetica');
xlim([280 1200]); ylim([0.6 1.4]);
set(gca, 'YTick', 0.6:0.2:1.4, 'FontName', 'Helvetica');
xlabel('Atmospheric CO_2 (ppm)', 'FontSize', 12, 'FontName', 'Helvetica');
ylabel('C4 Relative Advantage', 'FontSize', 12, 'FontName', 'Helvetica');
title('(d)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(gca, 'FontSize', 12);

% Adjust subplot positions for better alignment
w1 = 0; hei = 0.475; wei = 0.575; hig = 0.52;
set(h1, 'Position', [w1 hei wei hig]);
set(c1, 'Position', [0.125-0.001 hei+0.035 0.35 0.02]);
set(h2, 'Position', [w1+wei-0.125 hei wei hig]);
set(c2, 'Position', [w1+wei-0.01 hei+0.035 0.35 0.02]);

set(h3, 'Position', [w1+0.17 0.1 0.25 0.3]);
set(h4, 'Position', [w1+wei+0.015 0.1 0.3 0.3]);

set(h1, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
set(h2, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background


% Export figure
filename = './Figures/Figure_3_C4_changes.png';
exportgraphics(fig, filename, 'Resolution', 300, 'BackgroundColor', 'w');

%%


close all;
fig = figure;
width = 9.5;
height = 5.5;
set(gcf,'position',[100 80 width*100 height*100])

A = nan(360,720); A(kland_05) = 1; kocean = isnan(A);
for i = 1:4
    eval(['h',num2str(i),'=subplot(2,2,i)'])
    
    % Select data based on subplot position
%     if i == 1; data = S_change; end
    if i == 1; data = T_change; end
    if i == 2; data = A_change; end
    if i == 3; data = v_change; end
    if i == 4; data = spi_change; end

    % Process data
    a = data + c4_diff_us;
    k = find(isnan(a)); 
    data(k) = nan; thres = ceil( max( [abs(prctile(data,5)), abs(prctile(data,5))]) );

    k = find(isnan(data)); data(k) = -thres;
    % Set up Robinson projection
    axesm('MapProjection', 'robinson', 'Frame', 'off', 'Grid', 'off', ...
          'MapLatLimit', [-60 90], 'MapLonLimit', [-180 180]);
    hold on;
    
    % Smooth and plot the data
    A = nan(360, 720);
    A(kland_05) = data;
    kk = find(isnan(A)); 
    A(kk) = nanmean(data);
    A_smooth = imgaussfilt(A, 0.75);
    A_smooth(kk) = nan;
    
    A(kocean) = nan;
    h = pcolorm(lat1, lon1, A);  % Note: pcolorm for mapping toolbox
    set(h, 'EdgeColor', 'none');
    
    % Styling
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    caxis([-thres thres]);
    colormap(color_use_with_grey);
    
    % Add coastline
    plotm(coastlat, coastlon, 'k-', 'LineWidth', 1.5);  % plotm for mapping toolbox
    
    % Titles
%     if i == 1; title('(a) \Delta Z_{TWS}', 'FontSize', 16, 'FontWeight', 'bold'); end
    if i == 1; title('(a) \Delta Z_{TA}', 'FontSize', 16, 'FontWeight', 'bold'); end
    if i == 2; title('(b) \Delta Z_{RS}', 'FontSize', 16, 'FontWeight', 'bold'); end
    if i == 3; title('(c) \Delta Z_{VPD}', 'FontSize', 16, 'FontWeight', 'bold'); end
    if i == 4; title('(d) \Delta Z_{PRE}', 'FontSize', 16, 'FontWeight', 'bold'); end
    
    eval(['c',num2str(i),'=colorbar(''southoutside'')'])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background

end

% Keep position settings
w1 = 0; hei = 0.55; wei = 0.5; hight = 0.4; sp = 0;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1 + wei + sp, hei, wei, hight]);

spacs = 0.16;
set(c1, 'Position', [w1+spacs, hei, wei*0.5, 0.02]);
set(c2, 'Position', [w1 + wei + sp+spacs, hei, wei*0.5, 0.02]);

 hei = 0.05; 
set(h3, 'Position', [w1, hei, wei, hight]);
set(h4, 'Position', [w1 + wei + sp, hei, wei, hight]);
set(c3, 'Position', [w1+spacs, hei, wei*0.5, 0.02]);
set(c4, 'Position', [w1 + wei + sp+spacs, hei, wei*0.5, 0.02]);

filename = './Figures/Figure_SI_all_factor_change.png';
exportgraphics(fig, filename, 'Resolution', 300);

%%
%% Bootstrapping the group means
k = ~isnan( c4_diff_us + spi_change + S_change + A_change + v_change  ); kkuse = kland_05(k); kk_ind = k;
% k = ~isnan( C4_change + spi_change + S_change + A_change + v_change  ); kkuse = kland_05(k); kk_ind = k;
y = S_change(k);  x2 = spi_change(k); x3 = A_change(k); x4 = v_change(k); x5 = T_change(k);
% train random forest
x6 = c4_diff_us(k);

X = [x6,x2,x3,x4,x5];

% Assume X is an n x 5 matrix (n observations, 5 predictors) and y is n x 1 vector
% Replace X and y with your actual data

% Example data (uncomment to test)
% n = 100;
% X = randn(n, 5); % 5 correlated columns
% X(:,2) = X(:,1) + 0.1*randn(n,1); % Induce correlation
% X(:,3) = X(:,2) + 0.2*randn(n,1);
% y = 2*X(:,1) + 1.5*X(:,2) + 0.5*X(:,3) + randn(n,1);

% Parameters
n = size(X, 1);         % Number of observations
p = size(X, 2);         % Number of predictors (5 in your case)
n_boot = 1000;          % Number of bootstrap samples
lambda = 0.1;           % Ridge regression penalty (tune this as needed)

% Standardize X and y (optional but recommended for ridge regression)
X_mean = mean(X);
X_std = std(X);
X_norm = (X - X_mean) ./ X_std;
y_mean = mean(y);
y_norm = y - y_mean;

% Initial ridge regression on full data
beta_full = ridge(y_norm, X_norm, lambda, 0); % 0 means no scaling of penalty

% Bootstrap sampling
boot_betas = zeros(n_boot, p+1); % Store bootstrap coefficients

for i = 1:n_boot
    % Generate bootstrap sample indices with replacement
    boot_idx = randsample(n, n, true);
    
    % Bootstrap samples
    X_boot = X_norm(boot_idx, :);
    y_boot = y_norm(boot_idx);
    
    % Ridge regression on bootstrap sample
    boot_betas(i, :) = ridge(y_boot, X_boot, lambda, 0);
end

% Calculate slopes and uncertainties
slopes = beta_full;                  % Slopes from full data
slope_std = std(boot_betas);         % Standard error from bootstrap
slope_ci = prctile(boot_betas, [2.5, 97.5]); % 95% confidence intervals


slopes_boot = boot_betas;
perc_95 = prctile(slopes_boot(:, 2), 95);


fig = figure;


histogram(slopes_boot(:, 2), 100, 'Normalization', 'probability', ...
    'FaceColor', '#4CAF50', 'EdgeColor', 'k', 'LineWidth', 0.1);
hold on;
xline(perc_95, 'b-', 'LineWidth', 2, 'Alpha', 0.8);
text(perc_95+0.0002, 0.9 * max(get(gca, 'YLim')), ...
    sprintf('95th percentile = %.3f', perc_95), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Color', 'b', 'FontSize', 12);
xlabel('Bootstrap Slope for \Delta C4');
ylabel('Probability');
% title('Distribution of Bootstrap Slopes for \Delta C4');
grid on;
set(gca, 'FontSize', 12);
xline(0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5);
% if perc_95 < 0
%     text(0.05, 0.95, '95th percentile significantly < 0', ...
%         'Units', 'normalized', 'Color', 'blue', 'FontSize', 10, ...
%         'FontWeight', 'bold');
% end
hold off;

xlim([-16 4]*10^-3)

% Export figure
filename = './Figures/FigureSI_slope_C4_changes.png';
exportgraphics(fig, filename, 'Resolution', 300, 'BackgroundColor', 'w');


