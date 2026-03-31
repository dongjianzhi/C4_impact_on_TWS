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

% SIF = A.csfi_mon;

A = matfile('/Users/jianzhidong/Documents/Research/C4transition/LAI_025degree/data_lai_col.mat');
time_lai = A.time;

[~,ib,ic] = intersect(time,time_lai);
LAI = nan(length(kland_05),length(time));
LAI(:,ib) = A.data_lai_col;

SIF = LAI;

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

[SIF_change] = cal_water_storage_change_warm(sif, ...
    SPI, time, lats, vpd_mon, vpdz,spi_thres,vpd_thres,vpd_thres2);

%%
k = ~isnan( C4_use + spi_change + S_change + A_change + v_change + SIF_change ); kkuse = kland_05(k); kk_ind = k;
y = S_change(k); x1 = C4_use(k); x2 = spi_change(k); x3 = A_change(k); x4 = v_change(k); x5 = T_change(k);

x6 = SIF_change(k);
% A = nan(360,720); A(kkuse) = v_sen; h = pcolor(A); set(h,'edgecolor','none')
%% train random forest 
X = [x1,x2,x3,x4,x5, x6];
imp_all = nan(50,6); R2 = nan(50,1);
y_org = y; X_org = X;


P = regress(y,X)
imp_all = nan(50,6);

for i = 1:50
    i
    nTrees = 100; % Number of trees in the forest

    k = sort( randsample(1:length(X),length(X),'true') );

    model = TreeBagger(nTrees, X(k,:), y(k), 'Method', 'regression', ...
        'OOBPredictorImportance', 'on');
    imp_all(i,:) = model.OOBPermutedPredictorDeltaError;
end

%%
figure('Units', 'inches', 'Position', [1 1 7.2 5]*1., 'Color', 'w');

a = [imp_all(:,1) imp_all(:,6) imp_all(:,2) imp_all(:,4) imp_all(:,5) imp_all(:,3) ];
a = a ./ sum(a, 2);
imp = nanmean(a); imp_s = nanstd(a);
h = barh(imp, 'FaceColor', [0.35 0.6 0.4], 'EdgeColor', 'k', 'LineWidth', 0.5, 'BarWidth', 0.7);
hold on;
h = errorbar(imp, 1:6, imp_s, 'horizontal', 'k', 'LineWidth', 1.5, 'CapSize', 3);
h.LineStyle = 'none';
% plot([mean(a(:,1)) mean(a(:,1))], [0 6], '--', 'Color', [0.8 0.6 0.4], 'LineWidth', 1.5);
set(gca, 'YTick', 1:6, 'YTickLabel', {'C4', 'LAI' ,'Z_P', 'VPD', 'Rs', 'T'}, ...
    'FontSize', 11, 'FontName', 'Helvetica', 'LineWidth', 0.5);
xlim([0 0.4]); ylim([0.5 5.5]);
% title('(d)', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');
xlabel('Relative Importance', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
box on;

set(gcf, 'Renderer', 'painters');


filename='./Figures/Mar_FigureS_with_LAI.png'; 
exportgraphics(gcf, filename, 'Resolution', 300);



