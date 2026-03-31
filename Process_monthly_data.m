clear
load F:\CPC_precp_daily\kland_05.mat
A = nan(360,720);
A(kland_05) = kland_05; h = pcolor(A); set(h,'edgecolor','none')
kuse = find( kland_05 == A(240,580) );

A(kland_05) = kland_05; h = pcolor(A); set(h,'edgecolor','none')
% kuse = find( kland_05 == A(310,500) );
kuse = find( kland_05 == A(150,410) );

% load Resampled_8day_long_vpd.mat
load Resampled_8day.mat
A = matfile('Temp_rs_8day.mat')
rs = A.Rs_8d; Tair = A.Tmean_8d;


dd = datevec(Time);
mm = datenum(dd(:,1),dd(:,2),1);
um = unique(mm);


load F:\CPC_precp_daily\cpc_pre_mon.mat

[~,ib,ic] = intersect(cpc_m,um);
cpc_mom = double( cpc_pre_mon(:,ib) );



load F:\Research\MSWEP_monthly\mswep_month_05.mat

for i = 1:length(um)
    k = find( mm == um(i) );
    csfi_mon(:,i) = nanmean( CSFI(:,k),2);
    geo_mon(:,i) = nanmean(GEO_SIF(:,k),2);
    vpd_mon(:,i) = nanmean(vpd_8d(:,k),2);
    et_mon(:,i) = nanmean(ET_8d(:,k),2);
    
    T_mon(:,i) = nanmean(Tair(:,k),2);
    RS_mon(:,i) = nanmean(rs(:,k),2);
end

load F:\Fluxcom_extended\Fluxcom_month_RS_metero_05_long.mat

[~,ib,ic] = intersect(um,time_um);
LE_flux = nan(length(kland_05),length(um)); LE_flux(:,ib) = LE_05;
H_flux = nan(length(kland_05),length(um)); H_flux(:,ib) = H_05;


save('monthly_data','csfi_mon','geo_mon','vpd_mon','et_mon','T_mon','RS_mon','Precp','um','H_flux','LE_flux','-v7.3')


