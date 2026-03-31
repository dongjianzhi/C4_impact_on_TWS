function [et_change] = cal_water_storage_change_warm(tws_z, SPI, time, lats, vpd_mon, vpdz,spi_thres, vpd_thres,vpd_thres2)
 

et_change = nan(length(tws_z),1);

num_thres = 2;

%For tropical regions (23.25°N to 23.25°S),all 12 months were included. 
varb = tws_z; k = find(SPI>spi_thres| vpd_mon<vpd_thres | vpd_mon> vpd_thres2); varb(k) = nan;
dd = datevec(time); ksummer = find(dd(:,2)>=1 & dd(:,2)<=12);
k_space = find(lats<= 23.25 & lats>=-23.25 );

data = varb(k_space,ksummer);
k = sum(~isnan(data),2) <num_thres; data(k,:) = nan;

et_change(k_space) = nanmean(data,2);
% et_change(k_space) = prctile(data,50,2);


% In the subtropics (23.25°N to 35°N and 23.25°S to 35°S), 
%  May–October data (in the Northern Hemisphere) or 
% November–April data (in the Southern Hemisphere).
varb = tws_z; k = find(SPI>spi_thres| vpd_mon<vpd_thres | vpd_mon> vpd_thres2); varb(k) = nan;
dd = datevec(time); ksummer = find(dd(:,2)>=5 & dd(:,2)<=10);
k_space = find(lats<= 35 & lats>=23.25 );

data = varb(k_space,ksummer);
k = sum(~isnan(data),2) <num_thres; data(k,:) = nan;
et_change(k_space) = nanmean(data,2);
% et_change(k_space) = prctile(data,50,2);


dd = datevec(time); ksummer = find(dd(:,2)>=11 | dd(:,2)<=4);
k_space = find(lats<= -23.25 & lats>=-35 );

data = varb(k_space,ksummer);
k = sum(~isnan(data),2) <num_thres; data(k,:) = nan;
et_change(k_space) = nanmean(data,2);
% et_change(k_space) = prctile(data,50,2);

% In temperate zones (>35°N and >35°S), 
% results were based on June-July-August data (in the Northern Hemisphere) 
% or December-January-February data (in the Southern Hemisphere). 
varb = tws_z; k = find(SPI>spi_thres| vpd_mon<vpd_thres | vpd_mon> vpd_thres2); varb(k) = nan;
dd = datevec(time); ksummer = find(dd(:,2)>=6 & dd(:,2)<=8);
k_space = find(lats >=35 );

data = varb(k_space,ksummer);
k = sum(~isnan(data),2) <num_thres; data(k,:) = nan;
et_change(k_space) = nanmean(data,2);
% et_change(k_space) = prctile(data,50,2);

dd = datevec(time); ksummer = find(dd(:,2)>=12 | dd(:,2)<=2);
k_space = find(lats<= -35 );

data = varb(k_space,ksummer);
k = sum(~isnan(data),2) <num_thres; data(k,:) = nan;
et_change(k_space) = nanmean(data,2);
% et_change(k_space) = prctile(data,50,2);

