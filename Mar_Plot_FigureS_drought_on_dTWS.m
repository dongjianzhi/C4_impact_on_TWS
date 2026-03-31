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

k = find( grass_fraction< 80 | isnan(grass_fraction) | irr > 1); c4_org(k) = nan;


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

pm = [Precp(:,2:end)+Precp(:,1:end-1),nan(length(kland_05),1)]/2;
[~, pz] = cal_mon_zcore_linear_by_month(time, pm);
SPI = pz;

%% Only keep the period of 2003 - 2019
k = find( time>= datenum(2003,1,1) & time<= datenum(2019,12,31));
vm = vm(:,k); SPI = SPI(:,k); tws = tws(:,k); dS = dS(:,k);
time = time(k);
%% Calculate Z score
% kick out the edges
load kland_025.mat
A = nan(720,1440);  A(kland_025) = 1;
B = A(1:2:end,:) + A(2:2:end,:); C = B(:,1:2:end) + B(:,2:2:end);
kclean = find(isnan(C) );

dsd_all = nan( length(kland_05), 5);

% tws(:,265:end) = nan;
tws_detrend = remove_linear_trend(time,dS);

for J = 1:5
    varb = tws_detrend;
    
    if J == 1;  k = find(SPI>0 |vm<1); end
    if J == 2;  k = find(SPI>0 |vm<0); end
    if J == 3;  k = find(SPI>0 |vm<1.5); end
    if J == 4;  k = find(SPI>-1 |vm<0); end
    if J == 5;  k = find(SPI>-1 |vm<1); end
    
    varb(k) = nan;
    
    dsd = nanmean(varb,2) - nanmean(tws_detrend,2);
%     k = find(dsd>0);   dsd(k) = nan;
    k = find( sum(~isnan(varb),2)<4); dsd(k) = nan;
    dsd_all(:,J) = dsd;
end



lats = lat1(kland_05);
W = cos( lats/180*pi);

k = find(isnan(c4_diff_us) );
dsd_no_cp = (dsd_all(k,:).*W(k) );
k = find(dsd_no_cp < prctile(dsd_no_cp,5) ); dsd_no_cp(k) = nan;
k = find(dsd_no_cp > prctile(dsd_no_cp,95) ); dsd_no_cp(k) = nan;

k = find(~isnan(c4_diff_us) & ~isnan(dsd) );
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
fig = figure
width=7.5
height=5
set(gcf,'position',[100 80 width*100 height*100])
load coastlines.mat

rd = brewermap(256,'OrRd'); tt = 2;
kk = round( linspace(1,220^tt,32).^(1/tt) );

h1 = subplot(2,1,1)

plot_mask_025(0.75)
hold on
A = nan(360,720); A(kland_05) = dsd_all(:,1); A(kclean) = nan;
h = pcolor(lon1,lat1,A); set(h,'edgecolor','none')
caxis([-5 0]); colormap(flipud(rd(kk,:) ))
c1 = colorbar;
ylim([-60 90]); grid on; set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
hold on
for i = 1:length(polygons)
    %     plot(polygons{i} ,'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
    a = polygons{i}.Vertices(:,1); b = polygons{i}.Vertices(:,2);
    a1 = nanmin(a):2:nanmax(a); b1 = nanmin(b):2:nanmax(b);
    [lon,lat] = meshgrid(a1,b1);
    IN = inpolygon(lon,lat,a,b); lon_use = lon(IN); lat_use = lat(IN);

    scatter(lon_use, lat_use, 3, 'k', 'filled');
end
load coastlines.mat
hold on; plot(coastlon,coastlat,'color',[1 1 1]*.25,'linewidth',1)
ylim([-60 90])

xlabel(c1,'cm','FontSize',12, 'FontName', 'Arial');

title('(a) E(dTWS/dt | Z_p<0, VPD>1.5) - E(dTWS/dt)', 'FontWeight', 'bold','FontName', 'Arial','FontSize',14);
set(gca,'FontSize', 14, 'FontName', 'Arial');
c1.TickLabels{end} = '\geq 0'
%-------------------------------------------

h2 = subplot(2,1,2)
hold on;
kkuse = [2 4 1 3 5];
group_shift = 0.15; % Offset for grouped bars
x = [1:5]; bar_width = 0.9;
hh = bar(x, dmeans(kkuse, :), 'grouped', 'BarWidth', bar_width);

% Use clean colors
hh(1).FaceColor = [0.2, 0.4, 0.8]; % Blue
hh(2).FaceColor = [0.8, 0.4, 0.2]; % Red

% Add error bars
errorbar(x - group_shift, dmeans(kkuse, 1), dstds(kkuse, 1), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(x + group_shift, dmeans(kkuse, 2), dstds(kkuse, 2), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

% Aesthetic improvements
grid on;
% set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6); % Light gridlines
set(gca, 'FontName', 'Arial', 'FontSize', 12); % Professional font
ylim([-0.8 0]); % Adjust y-axis limits
xlim([0.5, 5.5]); % Adjust x-axis limits
xticks(x); % Set x-axis ticks
xticklabels({'Z_P<0', 'Z_P<-1', 'Z_P<0 & VPD>1', ...
    'Z_P<0 & VPD>1.5', 'Z_P<-1 & VPD>1.5'}); % Label x-axis
ylabel('cm'); % Label y-axis
% xlabel('Group', 'FontWeight', 'Bold'); % Label x-axis
% xticklabels({'SPI<-1', 'SPI<0 & VPD>1', ...
%     'SPI<0 & VPD>1.5', 'SPI<-1 & VPD>1.5'}); % Label x-axis

% Add legend
legend({'Grass/crop', 'Others'}, 'Location', 'SouthWest', 'FontSize', 12, 'Box', 'off');
% Finalize the plot
box on; % Add border to the axes
% set(gcf, 'Color', 'w'); % Set white background
title('(b) dTWS/dt during droughts', 'FontWeight', 'Bold', 'FontSize', 14,'FontName', 'Arial'); % Add title
% ylabel('cm','FontWeight', 'Bold', 'FontSize', 12,'FontName', 'Arial');
ylim([-1.6 0])
set(gca,'YTick',[-1.6:0.4:0])
%-------------------------------------------

set(h1,'position',[0.15 0.51 0.7 0.425])
set(h2,'position',[0.15 0.1 0.725 0.35])

set(c1,'position',[0.855 0.51 0.02 0.425])
% 
% set(gcf,'InvertHardcopy','on');
% set(gcf,'paperUnits','inches');
% papersize=get(gcf,'paperSize');
% left=(papersize(1)-width)/2;
% bottom=(papersize(2)-height)/2;
% myfiguresize=[left,bottom,width,height]'
% set(gcf,'paperPosition',myfiguresize);
% 
% filename='./Figures/Figure_drought_impacts_TWS'; 
% print(filename,'-dpng','-r300')

filename='./Figures/FigureS2_drought_impacts_dTWS.png'; 
exportgraphics(fig, filename, 'Resolution', 300);



