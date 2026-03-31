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

k = find( grass_fraction< 50 | isnan(grass_fraction) | irr > 1); c4_org(k) = nan;


c4_diff_us = c4_diff(kland_05);  c4_org_us = c4_org(kland_05);
a = prctile(c4_diff_us,95); b = prctile(c4_diff_us,5);
k = find(c4_diff_us>= a | c4_diff_us<=b); c4_diff_us(k) = nan;

% Step 4: Map Pixel Coordinates to Geographic Coordinates
binaryImage = imbinarize(c4_org);

% Step 3: Trace Boundaries



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
k = find( time>= datenum(2003,1,1) & time<= datenum(2021,12,31));
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


    [et_change] = cal_water_storage_change_warm_count(tws_detrend, SPI, time, lats, vpd_mon, vpdz,spi_thres, vpd_thres,20);
     dsd_all(:,J) = et_change;
end


dd_use = dsd_all(:,[2 4 1 3 5]);
%%
% Create the bar plot
rd = viridis;

close all;
fig = figure;
width = 8;
height = 5.5;
set(gcf, 'Position', [100 80 width*100 height*100]);
load coastlines.mat;

rd = brewermap(256, 'OrRd');
tt = 2;
tt = 1.8;
kk = round(linspace(1, 256^tt, 220).^(1/tt));

% Subplot 1: Map
for i = 1:5
    eval(['h',num2str(i),'=subplot(2, 3, i)'])

    axesm('robinson', 'MapLatLimit', [-60 90], 'Grid', 'off', 'Frame', 'off');
    set(gca, 'FontName', 'Arial', 'FontSize', 12);

    % Use existing 2D lat1 and lon1 directly
    A = nan(360, 720);
    A(kland_05) = dd_use(:, i);


    % Plot the data with its colormap
    h = pcolorm(lat1, lon1, A);
    k = find(A > 5); A(k) = 5;
    set(h, 'EdgeColor', 'none');
    caxis([0 100]);
    % colormap(h1, [[1 1 1]*0.75; flipud(rd(kk, :))]);
    colormap(rd(kk, :));

    % Add grey shading only over land NaN regions
    hold on;
    plotm(coastlat, coastlon, 'Color', [0.25 0.25 0.25], 'LineWidth', 1);
    eval(['c',num2str(i), '=colorbar(''southoutside'')'])


    tightmap;
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');

    if i == 1; title('(a) Z_P<0 ', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14); end
    if i == 2; title('(b) Z_P<-1', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14); end
    if i == 3; title('(c) Z_P<0 & VPD>1kpa', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14); end
    if i == 4; title('(d) Z_P<0 & VPD>1.5kpa', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14); end
    if i == 5; title('(e) Z_P<-1 & VPD>1.5kpa', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14); end
    set(gca, 'FontSize', 14, 'FontName', 'Arial');
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');


end


% Adjust subplot positions
w1 = 0; hei = 0.55; wei = 0.37; hight = 0.375; sp = -0.05;
set(h1, 'Position', [w1, hei, wei, hight]);
set(h2, 'Position', [w1+wei+sp, hei, wei, hight]);
set(h3, 'Position', [w1+(wei+sp)*2, hei, wei, hight]);

ssp = 0.08;
set(c1, 'Position', [w1+ssp, hei+0.05, wei*0.6, 0.015]);
set(c2, 'Position', [w1+ssp+wei+sp, hei+0.05, wei*0.6, 0.015]);
set(c3, 'Position', [w1+ssp+(wei+sp)*2, hei+0.05, wei*0.6, 0.015]);

hei = 0.1; 
set(h4, 'Position', [w1, hei, wei, hight]);
set(h5, 'Position', [w1+wei+sp, hei, wei, hight]);
set(c4, 'Position', [w1+ssp, hei+0.05, wei*0.6, 0.015]);
set(c5, 'Position', [w1+ssp+wei+sp, hei+0.05, wei*0.6, 0.015]);


filename='./Figures/Figure1_counts.png'; 
exportgraphics(fig, filename, 'Resolution', 300);

%%
load /Users/jianzhidong/Documents/Research/C4transition/Modis_land_cover/LC_fraction_forest.mat

forests = LC_fraction(kland_05);

dd = dd_use( forests > 10,:);




