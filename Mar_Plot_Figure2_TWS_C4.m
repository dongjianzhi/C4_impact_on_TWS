clear
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
 
%% get a mask of significant C4 changes
load GIA_based_IA.mat
irr = nan(360,720); irr(kland_05) = IA;

load('../Modis_land_cover/LC_fraction.mat')
grass_fraction = LC_fraction;

name = 'C4_distribution_NUS_v2.2.nc'
a = flipud( ncread(name,'C4_area') );

C4 = nanmean(a,3);

k = find( grass_fraction< 50 | isnan(grass_fraction) | irr>1); C4(k) = nan;
C4_use = C4(kland_05);
%% Calculate Z score
% dS = [nan(length(kland_05),1) tws(:,2:end)-tws(:,1:end-1) ];

% [~, tws_z] = cal_mon_zcore_linear_by_month(time, tws);
[dSa, dSz] = cal_mon_zcore_linear_by_month(time, tws);
[Aa, Az] = cal_mon_zcore_linear_by_month(time, Rs_mon);
[Ta, Tz] = cal_mon_zcore_linear_by_month(time, T_mon);
[~, SPI] = cal_mon_zcore_linear_by_month(time, Precp);
[~, sif] = cal_mon_zcore_linear_by_month(time, SIF);


%% Calculate conditional TWS
% A = matfile('SPI3_data.mat'); SPI = A.SPI3;
lats = lat1(kland_05);
[vpda, vpdz] = cal_mon_zcore_linear_by_month(time, vpd_mon);

%%  calculate the Z score of different factors
lats = lat1(kland_05);
vpd_thres = 1.; spi_thres = 0; vpd_thres2 = 30;


%% only use a section of data
k = find( time>= datenum(2003,1,1) & time<= datenum(2019,12,31));
dSz = dSz(:,k); SPI = SPI(:,k); Az = Az(:,k); 
vpdz = vpdz(:,k); Tz = Tz(:,k);  time = time(k); vpd_mon = vpd_mon(:,k);
%%
[S_change] = cal_water_storage_change_warm(dSz, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);
[spi_change] = cal_water_storage_change_warm(SPI, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);
[A_change] = cal_water_storage_change_warm(Az, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);
[v_change] = cal_water_storage_change_warm(vpdz, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);
[T_change] = cal_water_storage_change_warm(Tz, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);



%%
k = ~isnan( C4_use + spi_change + S_change + A_change + v_change  ); kkuse = kland_05(k); kk_ind = k;
y = S_change(k); x1 = C4_use(k); x2 = spi_change(k); x3 = A_change(k); x4 = v_change(k); x5 = T_change(k);

% A = nan(360,720); A(kkuse) = v_sen; h = pcolor(A); set(h,'edgecolor','none')
%% train random forest 
X = [x1,x2,x3,x4,x5];
imp_all = nan(50,5); R2 = nan(50,1);
y_org = y; X_org = X;


P = regress(y,X)
imp_all = nan(50,5);

for i = 1:50
    i
    nTrees = 100; % Number of trees in the forest

    k = sort( randsample(1:length(X),length(X),'true') );

    model = TreeBagger(nTrees, X(k,:), y(k), 'Method', 'regression', ...
        'OOBPredictorImportance', 'on');
    imp_all(i,:) = model.OOBPermutedPredictorDeltaError;
end

%%

low_c4 = find(x1<median(x1));
high_c4 = find(x1>median(x1));

thres = 0.1;

x2_low = x2(low_c4); x3_low = x3(low_c4); x4_low = x4(low_c4); x5_low = x5(low_c4);
x2_high = x2(high_c4); x3_high = x3(high_c4); x4_high = x4(high_c4); x5_high = x5(high_c4);
y_low = y(low_c4); y_high = y(high_c4);


y_high_corr = nan(size(y_low));
y_low_corr = nan(size(y_low));

for i = 1:length(x2_low)
    %     try
    k =  abs(x2_high - x2_low(i))<thres & abs(x3_high - x3_low(i))<thres &...
        abs(x4_high - x4_low(i))<thres & abs(x5_high - x5_low(i))<thres;

    ys = y_high(k); xs = [x2_high(k) x3_high(k) x4_high(k) x5_high(k)];
    y_high_corr(i) = nanmean(y_high(k) );


    k =  abs(x2_low - x2_low(i))<thres & abs(x3_low - x3_low(i))<thres &...
        abs(x4_low - x4_low(i))<thres & abs(x5_low - x5_low(i))<thres;
    y_low_corr(i) = nanmean(y_low(k) );
end

%%
positions = [0, -0.2, -0.5, -0.8, -1];  % Uneven spacing, more focus on dry range
colors = [
    255/255, 248/255, 231/255;   % Pale Cream (#FFF8E7) - Neutral (0)
    250/255, 235/255, 215/255;   % Light Tan (#FAEBD7) - Slightly Dry (-0.2)
    198/255, 142/255,  23/255;   % Medium Brown (#C68E17) - Moderately Dry (-0.5)
    139/255,  90/255,  43/255;   % Rich Brown (#8B5A2B) - Dry (-0.8)
     60/255,  47/255,  47/255    % Dark Brown (#3C2F2F) - Very Dry (-1)
];
% Number of interpolation points
n = 256;  % Smooth gradient with 256 levels

% Create a custom colormap by interpolating with uneven spacing
cmap = interp1(positions, colors, linspace(0, -1, n), 'linear');

%%
%%
% Plot figure
load coastlines.mat
close all


% Set figure size in inches (double-column width for Nature)
figure('Units', 'inches', 'Position', [1 1 7.2 5]*1.5, 'Color', 'w');

% Define a refined nature-inspired colormap
nature_cmap = [0.2, 0.4, 0.2; 0.35, 0.6, 0.4; 0.6, 0.8, 0.5; ...
               0.9, 0.85, 0.7; 0.8, 0.6, 0.4; 0.5, 0.3, 0.2];
color_use = interp1(1:6, nature_cmap, linspace(1, 6, 256), 'pchip');
gray_color = [0.7 0.7 0.7]  ; % Light gray for NaN regions
color_use_with_gray = [gray_color; flipud(color_use)];

color_use_with_gray2 = [gray_color; flipud(cmap)];


% Subplot (a): Z_{TWS} with Robinson projection
h1 = subplot(2, 5, [1 2]);
axesm('robinson', 'MapLatLimit', [-60 90], 'Frame', 'off', 'Grid', 'off');
A = nan(360, 720);
k = find(isnan(C4_use + S_change)); S_change(k) = nan;
A(kland_05) = S_change;
k_nan = isnan(S_change);

ths = -0.9; 
A(kland_05(k_nan)) = ths; 
k_fill = isnan(A); A = max(ths,A);
A(k_fill) = nan; % Non-land areas remain NaN (transparent)
Iblur = imgaussfilt(A, 1);
h = pcolorm(lat1, lon1, A); set(h, 'EdgeColor', 'none');
colormap(h1, color_use_with_gray2);
plotm(coastlat, coastlon, 'k-', 'LineWidth', 1);
c1 = colorbar('southoutside', 'FontSize', 11, 'FontName', 'Helvetica', 'LineWidth', 0.5);
title('(a) Z_{TWSd}', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(c1, 'Ticks', [-0.9:0.3:-0.], 'TickLength', 0.015);
tightmap;
set(h1, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
caxis([-0.9 -0.0]); c1.TickLabels{end} ='\geq 0'



A = matfile("../Modis_land_cover/LC_fraction_ice_desert.mat")
data = A.LC_fraction; ice_des = data(kland_05);
lats = lat1(kland_05);
W = cos( lats/180*pi);
k1 = ~isnan(S_change);
k2 = isnan(ice_des);
sum(W(k1))./sum(W(k2))



% Subplot (b): C4 fraction with Robinson projection
h2 = subplot(2, 5, [3 4]);
axesm('robinson', 'MapLatLimit', [-60 90], 'Frame', 'off', 'Grid', 'off');
k = find(isnan(C4_use + S_change)); C4_use(k) = nan;
A = nan(360, 720); k = find(isnan(C4_use)); c4 = C4_use; c4(k) = nan;
A(kland_05) = c4;
k_nan = isnan(c4); A(kland_05(k_nan)) = 0;
k_fill = isnan(A);
A(k_fill) = nan; % Non-land areas remain NaN (transparent)
h = pcolorm(lat1, lon1, A); set(h, 'EdgeColor', 'none');
colormap(h2, color_use_with_gray);
caxis([10 80]);
plotm(coastlat, coastlon, 'k-', 'LineWidth', 1);
c2 = colorbar('southoutside', 'FontSize', 11, 'FontName', 'Helvetica', 'LineWidth', 0.5);
xlabel(c2, 'Percentage (%)', 'FontSize', 11, 'FontName', 'Helvetica');
title('(b) C4 fraction', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(c2, 'Ticks', [10:20:80], 'TickLength', 0.015);
tightmap;
set(h2, 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Transparent background
caxis([0 80])
set(c2,'xtick',[0:20:80])

% Subplot (c): Scatter plot with regression
h3 = subplot(2, 5, 5);
% Define x and y (assuming x1 and y are already defined)
x = x1; 
x_range = linspace(min(x), max(x), 100); % For plotting smooth lines

% Prepare design matrix for quantile regression
X = [ones(length(x), 1), x]; % Design matrix with intercept and x
X_range = [ones(length(x_range), 1), x_range']; % Design matrix for x_range

% Define quantiles
quantiles = [0.01 0.10]; % 10th, 25th, 50th percentiles
y_quantiles = zeros(length(x_range), length(quantiles));
slopes = zeros(length(quantiles), 1); % To store slopes
t_values = zeros(length(quantiles), 1); % To store t-values

% Perform quantile regression
n = length(x);
opts = optimoptions('linprog', 'Display', 'off');
texts = {'1th Slope = ','10th Slope = ','50th Slope = '}

for i = 1:length(quantiles)
    tau = quantiles(i); % Quantile level
    
    % Quantile regression
    f = [zeros(2, 1); tau * ones(n, 1); (1 - tau) * ones(n, 1)]; % Objective
    Aeq = [X, eye(n), -eye(n)]; % Equality: X * beta + u - v = y
    beq = y; % Right-hand side: y
    lb = [-Inf; -Inf; zeros(2*n, 1)]; % Lower bounds
    
    beta_uv = linprog(f, [], [], Aeq, beq, lb, [], opts);
    beta = beta_uv(1:2); % Intercept and slope
    slopes(i) = beta(2); % Record slope
    
    % Compute fitted values and residuals
    y_fit = X * beta;
    residuals = y - y_fit;
    
    % Estimate density at the quantile (sparsity)
    h = 1.06 * std(residuals) * n^(-1/5); % Bandwidth (Silverman's rule)
    f_hat = mean(normpdf(residuals / h)) / h; % Kernel density estimate
    
    % Asymptotic standard error of slope
    V = inv(X' * X); % Variance-covariance matrix (unscaled)
    se_slope = sqrt(tau * (1 - tau) / (n * f_hat^2) * V(2, 2));
    t_values(i) = slopes(i) / se_slope;
    
    % Compute quantile line over x_range
    y_quantiles(:, i) = X_range * beta;
end

% Mean regression (for reference) with t-test
p = polyfit(x, y, 1);
y_pred = polyval(p, x_range);
slope_mean = p(1);
y_fit_mean = p(1) * x + p(2);
residuals_mean = y - y_fit_mean;
MSE = sum(residuals_mean.^2) / (n - 2); % Mean squared error
se_slope_mean = sqrt(MSE / sum((x - mean(x)).^2)); % Standard error of slope
t_value_mean = slope_mean / se_slope_mean; % t-value for mean slope


hold on
% Scatter plot
scatter(x, y, 30, 'MarkerFaceColor', [0.4 0.6 0.4]*0+0.8, 'MarkerEdgeColor', [0.2 0.3 0.2]*0+0.6, ...
    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Plot regression lines
line_colors = [[0.45 0.25 0.1]; [0.5 0.3 0.2]; [0.6 0.4 0.3]; [0.7 0.5 0.4]]*0;
plot(x_range, y_pred, 'Color', line_colors(1, :), 'LineWidth', 3); % Mean
plot(x_range, y_quantiles(:, 1), '--', 'Color', line_colors(2, :), 'LineWidth', 2); % 10th
plot(x_range, y_quantiles(:, 2), '--', 'Color', line_colors(3, :), 'LineWidth', 2); % 25th

% Customize axes
set(gca, 'Color', 'w', 'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1, ...
    'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.9 0.9 0.9], 'GridAlpha', 0.7);
set(gca, 'XTick', 0:20:100, 'YTick', -1.5:0.5:0.5);
xlim([0 100]); ylim([-1.5 0.5]);

% Labels and title
xlabel('C4 fraction (%)', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('Z_{TWSd}', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
title('(c)', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

% Add text annotations near regression lines
x_text = 20; % Position text near x = 70
text_offset = 0.05; % Vertical offset
t_crit = tinv(0.975, n-2); % Critical t-value for α = 0.05, two-tailed
for i = 1:length(quantiles)
    slope_str = sprintf('%.3f', slopes(i));
    if abs(t_values(i)) > t_crit % Significance test
        slope_str = [slope_str, '*'];
    end
    y_pos = y_quantiles(find(x_range >= x_text, 1), i) + (i-3.5) * text_offset;
    text(x_text, y_pos, [texts{i}, slope_str], ...
        'FontSize', 12, 'FontName', 'Arial', 'Color', line_colors(i+1, :));
end

% Mean slope annotation with star if significant
slope_mean_str = sprintf('%.3f', slope_mean);
if abs(t_value_mean) > t_crit
    slope_mean_str = [slope_mean_str, '*'];
end
y_pos_mean = y_pred(find(x_range >= x_text, 1));
text(x_text, y_pos_mean + (i+2) * text_offset, ['Mean Slope = ', slope_mean_str], ...
    'FontSize', 12, 'FontName', 'Arial', 'Color', line_colors(1, :));
box on


% Subplot (d): Bar plot
h5 = subplot(2, 5, [7 8]);
a = [imp_all(:,1) imp_all(:,2) imp_all(:,4) imp_all(:,5) imp_all(:,3)];
a = a ./ sum(a, 2);
imp = nanmean(a); imp_s = nanstd(a);
h = barh(imp, 'FaceColor', [0.35 0.6 0.4], 'EdgeColor', 'k', 'LineWidth', 0.5, 'BarWidth', 0.7);
hold on;
h = errorbar(imp, 1:5, imp_s, 'horizontal', 'k', 'LineWidth', 1.5, 'CapSize', 3);
h.LineStyle = 'none';
% plot([mean(a(:,1)) mean(a(:,1))], [0 6], '--', 'Color', [0.8 0.6 0.4], 'LineWidth', 1.5);
set(gca, 'YTick', 1:5, 'YTickLabel', {'C4', 'Z_P', 'VPD', 'Rs', 'T'}, ...
    'FontSize', 11, 'FontName', 'Helvetica', 'LineWidth', 0.5);
xlim([0 0.4]); ylim([0.5 5.5]);
title('(d)', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');
xlabel('Relative Importance', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
box on;

h6 = subplot(2, 5, 9);

data = matfile("ZL_slopes.mat");
slopes_boot = data.slopes_boot;
perc_95 = data.perc_95;

% Histogram with a bolder, earthy brown color
histogram(slopes_boot(:, 2), 100, 'Normalization', 'probability', ...
    'FaceColor', [1 1 1]*0.8, ... % Rich brown for a strong, natural look
    'EdgeColor', [1 1 1]*0.7, ... % Slightly lighter brown for soft edges
    'LineWidth', 0.3, ... % Thin edges to keep it clean
    'FaceAlpha', 0.9); % Higher opacity for emphasis
hold on;

% Vertical lines for P95 and zero
xline(perc_95, '-', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.5, 'Alpha', 0.9, 'Label', sprintf('P95 = %.3f', perc_95), ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'FontSize', 11, 'FontName', 'Helvetica');
xline(0, '--', 'Color', [0.5 0.3 0.2], 'LineWidth', 1.5, 'Alpha', 0.7); % Matches histogram for consistency

% Customize axes to match other subplots
set(gca, 'FontSize', 11, 'FontName', 'Helvetica', 'LineWidth', 0.5, ...
    'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.9 0.9 0.9], 'GridAlpha', 0.7);
xlabel('Bootstrapped Slopes', 'FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
title('(e)', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');

% Set limits and ticks for clarity
xlim([-0.015 0]); % Consistent with original scaling
ylim([0 0.045]); % Slightly expanded for visibility
set(gca, 'XTick', -0.015:0.005:0, 'YTick', 0:0.01:0.04); % Granular ticks

% Add a subtle box
box on;
hold off;


% Adjust subplot positions

% Adjust subplot positions
w1 = -0.02; hei = 0.55; wei = 0.52; hight = 0.375; sp = -0.04;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1+wei+sp, hei, wei, hight]);
set(c1, 'Position', [w1 + wei*0.18, hei, wei*0.7, 0.02]);
set(c2, 'Position', [w1+wei*1.18+sp, hei, wei*0.7, 0.02]);

w1 = 0.075; hei = 0.1; wei = 0.22; hight = 0.325; sp = 0.02;
set(h3, 'Position', [w1, hei, wei*1.2, hight]);
set(h5, 'Position', [w1+1.15*(wei+0.05), hei, wei*1.2, hight]);
set(h6, 'Position', [w1+2.5*(wei+0.04), hei, wei*1.15, hight]);

% set(c1,'xtick',[-1:0.2:0.6])


% Final polish
set(gcf, 'Renderer', 'painters');


filename='./Figures/Mar_Figure_2_TWS_droughts.png'; 
exportgraphics(gcf, filename, 'Resolution', 300);


