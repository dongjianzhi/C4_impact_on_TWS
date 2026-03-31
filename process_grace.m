clear
A = matfile('monthly_data.mat');
time_mon = A.um;

name = 'F:\GRACE\CSR_TWS_05.mat'
load(name)
time_tws = time;

name = 'F:\GRACE\gsfc_TWS_05.mat'
load(name)
time_gsfc = time;

name = 'F:\GRACE\JPL_TWS_05.mat'
load(name)


load F:\CPC_precp_daily\kland_05.mat
clear tws
for i = 1:length(time_tws)
    data = TWS_CSR(:,:,i);
    a = data(kland_05);
    tws(:,i) = a;
end

dd = datevec(time_tws); time = datenum(dd(:,1),dd(:,2),1);
csr_tws = nan(length(tws),length(time_mon));
[~,ib,ic] = intersect(time,time_mon);
csr_tws(:,ic) = tws(:,ib);


%%
clear tws
for i = 1:length(time_tws)
    data = TWS_gsfc(:,:,i);
    a = data(kland_05);
    tws(:,i) = a;
end
dd = datevec(time_tws); time = datenum(dd(:,1),dd(:,2),1);
gsfc_tws = nan(length(tws),length(time_mon));
[~,ib,ic] = intersect(time,time_mon);
gsfc_tws(:,ic) = tws(:,ib);

%%

name = 'F:\GRACE\JPL_TWS_05.mat'
load(name)

clear tws
for i = 1:length(time)
    data = TWS_JPL(:,:,i);
    a = data(kland_05);
    tws(:,i) = a;
end
dd = datevec(time); time = datenum(dd(:,1),dd(:,2),1);
JPL_tws = nan(length(tws),length(time_mon));
[~,ib,ic] = intersect(time,time_mon);
JPL_tws(:,ic) = tws(:,ib);


save('GRACE_CSR','csr_tws','gsfc_tws','JPL_tws','-v7.3')
