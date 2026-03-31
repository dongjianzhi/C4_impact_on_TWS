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

% tws(:,265:end) = nan;
tws_detrend = remove_linear_trend(time,tws);
lats = lat1(kland_05);
for J = 1:5
    
    
    if J == 1; spi_thres = 0; vpd_thres = 1; end
    if J == 2; spi_thres = 0; vpd_thres = 0; end
    if J == 3; spi_thres = 0; vpd_thres = 1.5; end
    if J == 4; spi_thres = -1; vpd_thres = 0; end
    if J == 5; spi_thres = -1; vpd_thres = 1.5; end


    [et_change] = cal_water_storage_change_warm(tws_detrend, SPI, time, lats, vpd_mon, vpdz,spi_thres, vpd_thres,20);
     dsd_all(:,J) = et_change;
end


A = matfile('../Modis_land_cover/LC_fraction_ice_only.mat')

LC_fraction = A.LC_fraction;
kice = LC_fraction(kland_05);

dsd_all(kice > 1,:) = nan;
dsd_all(:,2) = [];
%%
% Create the bar plot
rd = viridis;

close all;
fig = figure;
width = 9;
height = 7.5;
set(gcf, 'Position', [100 80 width*100 height*100]);
load coastlines.mat;

rd = brewermap(256, 'OrRd');
tt = 2;
tt = 1.8;
kk = round(linspace(1, 256^tt, 220).^(1/tt));

% Subplot 1: Map
for J = 1:4
    eval(['h',num2str(J),'=subplot(2,2,J);'])

    axesm('robinson', 'MapLatLimit', [-60 90], 'Grid', 'off', 'Frame', 'off');
    set(gca, 'FontName', 'Arial', 'FontSize', 12);

    % 1. Paint ALL land grey first
    % This uses the built-in 'landareas' shapefile
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(land, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');

    hold on;

    % 2. Prepare your data
    A = nan(360, 720);
    A(kland_05) = dsd_all(:, J);
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
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');


    if J ==1
        title('(a) E(TWSd | Z_p<0, VPD>1.0)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14);
    end

    if J ==2
        title('(b) E(TWSd | Z_p<0, VPD>1.5)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14);
    end

    if J ==3
        title('(c) E(TWSd | Z_p<-1, VPD>0)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14);
    end

    if J ==4
        title('(c) E(TWSd | Z_p<-1, VPD>1.5)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 14);
    end

    set(gca, 'FontSize', 14, 'FontName', 'Arial');
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');

end
colormap(h1, [[1 1 1]*0.75; flipud(rd(kk, :))]);
colormap(h2, [[1 1 1]*0.75; flipud(rd(kk, :))]);
colormap(h3, [[1 1 1]*0.75; flipud(rd(kk, :))]);
colormap(h4, [[1 1 1]*0.75; flipud(rd(kk, :))]);

c1 = colorbar('southoutside');
xlabel(c1, 'cm', 'FontSize', 12, 'FontName', 'Arial');
c1.TickLabels{end} = '\geq 0';

w1 = -0.05; hei = 0.475; wei = 0.57; hig = 0.52; spc = 0.1;
set(h1, 'Position', [w1 hei wei hig]);
set(h2, 'Position', [w1+wei-spc hei wei hig]);


hei = hei - 0.4;
set(h3, 'Position', [w1 hei wei hig]);
set(h4, 'Position', [w1+wei-spc hei wei hig]);

c1 = colorbar('southoutside')
set(c1,'position',[0.2 hei+0.05 wei 0.02])

xlabel(c1, 'cm', 'FontSize', 12, 'FontName', 'Arial');
c1.TickLabels{end} = '\geq 0';


filename='./Figures/FigureSI_drought_impacts_TWS.png'; 
exportgraphics(fig, filename, 'Resolution', 300);