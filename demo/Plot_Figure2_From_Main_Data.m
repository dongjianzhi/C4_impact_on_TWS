%% Plot Figure 2 from Main_Figures_Data.mat
clear; close all;

script_dir = fileparts(mfilename('fullpath'));
source_dir = '/Users/jianzhidong/Documents/Research/C4transition/Remote_sensing_SIF';
addpath(source_dir);
figure_dir = fullfile(script_dir,'Generated_Figures');
if ~exist(figure_dir,'dir'); mkdir(figure_dir); end

S = load(fullfile(script_dir,'Main_Figures_Data.mat'),'Figure2');
F = S.Figure2;
lon1 = F.longitude;
lat1 = F.latitude;
kland_05 = F.land_index;
C4_use = F.vegetation_fraction;
S_change = F.storage_change;
x1 = F.valid_vegetation_fraction;
y = F.valid_storage_change;
imp_all = F.predictor_importance;
slopes_boot = F.bootstrap_slopes;
perc_95 = F.slope_percentile;

positions = [0 -0.2 -0.5 -0.8 -1];
colors = [
    255/255 248/255 231/255;
    250/255 235/255 215/255;
    198/255 142/255 23/255;
    139/255 90/255 43/255;
    60/255 47/255 47/255];
n = 256;
cmap = interp1(positions,colors,linspace(0,-1,n),'linear');

load coastlines.mat;
fig = figure('Units','inches','Position',[1 1 7.2 5]*1.5,'Color','w');
nature_cmap = [0.2 0.4 0.2; 0.35 0.6 0.4; 0.6 0.8 0.5; ...
               0.9 0.85 0.7; 0.8 0.6 0.4; 0.5 0.3 0.2];
color_use = interp1(1:6,nature_cmap,linspace(1,6,256),'pchip');
gray_color = [0.7 0.7 0.7];
color_use_with_gray = [gray_color; flipud(color_use)];
color_use_with_gray2 = [gray_color; flipud(cmap)];

% Subplot (a)
h1 = subplot(2,5,[1 2]);
axesm('robinson','MapLatLimit',[-60 90],'Frame','off','Grid','off');
A = nan(360,720);
k = find(isnan(C4_use + S_change));
S_change(k) = nan;
A(kland_05) = S_change;
k_nan = isnan(S_change);
ths = -0.9;
A(kland_05(k_nan)) = ths;
k_fill = isnan(A);
A = max(ths,A);
A(k_fill) = nan;
Iblur = imgaussfilt(A,1); %#ok<NASGU>
h = pcolorm(lat1,lon1,A);
set(h,'EdgeColor','none');
colormap(h1,color_use_with_gray2);
plotm(coastlat,coastlon,'k-','LineWidth',1);
c1 = colorbar('southoutside','FontSize',11,'FontName','Helvetica','LineWidth',0.5);
title('(a) Z_{TWSd}','FontSize',14,'FontName','Helvetica','FontWeight','bold');
set(c1,'Ticks',-0.9:0.3:0,'TickLength',0.015);
tightmap;
set(h1,'XColor','none','YColor','none','Color','none');
caxis([-0.9 0]);
c1.TickLabels{end} = '\geq 0';

% Subplot (b)
h2 = subplot(2,5,[3 4]);
axesm('robinson','MapLatLimit',[-60 90],'Frame','off','Grid','off');
k = find(isnan(C4_use + S_change));
C4_use(k) = nan;
A = nan(360,720);
c4 = C4_use;
A(kland_05) = c4;
k_nan = isnan(c4);
A(kland_05(k_nan)) = 0;
k_fill = isnan(A);
A(k_fill) = nan;
h = pcolorm(lat1,lon1,A);
set(h,'EdgeColor','none');
colormap(h2,color_use_with_gray);
caxis([10 80]);
plotm(coastlat,coastlon,'k-','LineWidth',1);
c2 = colorbar('southoutside','FontSize',11,'FontName','Helvetica','LineWidth',0.5);
xlabel(c2,'Percentage (%)','FontSize',11,'FontName','Helvetica');
title('(b) C4 fraction','FontSize',14,'FontName','Helvetica','FontWeight','bold');
set(c2,'Ticks',10:20:80,'TickLength',0.015);
tightmap;
set(h2,'XColor','none','YColor','none','Color','none');
caxis([0 80]);
set(c2,'xtick',0:20:80);

% Subplot (c): original quantile-regression code
h3 = subplot(2,5,5);
x = x1;
x_range = linspace(min(x),max(x),100);
X = [ones(length(x),1),x];
X_range = [ones(length(x_range),1),x_range'];
quantiles = [0.01 0.10];
y_quantiles = zeros(length(x_range),length(quantiles));
slopes = zeros(length(quantiles),1);
t_values = zeros(length(quantiles),1);
n = length(x);
opts = optimoptions('linprog','Display','off');
texts = {'1th Slope = ','10th Slope = ','50th Slope = '};
for i = 1:length(quantiles)
    tau = quantiles(i);
    f = [zeros(2,1); tau*ones(n,1); (1-tau)*ones(n,1)];
    Aeq = [X,eye(n),-eye(n)];
    beq = y;
    lb = [-Inf;-Inf;zeros(2*n,1)];
    beta_uv = linprog(f,[],[],Aeq,beq,lb,[],opts);
    beta = beta_uv(1:2);
    slopes(i) = beta(2);
    y_fit = X*beta;
    residuals = y-y_fit;
    h = 1.06*std(residuals)*n^(-1/5);
    f_hat = mean(normpdf(residuals/h))/h;
    V = inv(X'*X);
    se_slope = sqrt(tau*(1-tau)/(n*f_hat^2)*V(2,2));
    t_values(i) = slopes(i)/se_slope;
    y_quantiles(:,i) = X_range*beta;
end
p = polyfit(x,y,1);
y_pred = polyval(p,x_range);
slope_mean = p(1);
y_fit_mean = p(1)*x+p(2);
residuals_mean = y-y_fit_mean;
MSE = sum(residuals_mean.^2)/(n-2);
se_slope_mean = sqrt(MSE/sum((x-mean(x)).^2));
t_value_mean = slope_mean/se_slope_mean;
hold on;
scatter(x,y,30,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.6 0.6 0.6], ...
    'LineWidth',0.5,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
line_colors = zeros(4,3);
plot(x_range,y_pred,'Color',line_colors(1,:),'LineWidth',3);
plot(x_range,y_quantiles(:,1),'--','Color',line_colors(2,:),'LineWidth',2);
plot(x_range,y_quantiles(:,2),'--','Color',line_colors(3,:),'LineWidth',2);
set(gca,'Color','w','FontSize',12,'FontName','Arial','LineWidth',1, ...
    'XGrid','on','YGrid','on','GridColor',[0.9 0.9 0.9],'GridAlpha',0.7);
set(gca,'XTick',0:20:100,'YTick',-1.5:0.5:0.5);
xlim([0 100]); ylim([-1.5 0.5]);
xlabel('C4 fraction (%)','FontSize',14,'FontName','Arial','FontWeight','bold');
ylabel('Z_{TWSd}','FontSize',14,'FontName','Arial','FontWeight','bold');
title('(c)','FontSize',16,'FontName','Arial','FontWeight','bold');
x_text = 20;
text_offset = 0.05;
t_crit = tinv(0.975,n-2);
for i = 1:length(quantiles)
    slope_str = sprintf('%.3f',slopes(i));
    if abs(t_values(i)) > t_crit; slope_str = [slope_str,'*']; end
    y_pos = y_quantiles(find(x_range >= x_text,1),i)+(i-3.5)*text_offset;
    text(x_text,y_pos,[texts{i},slope_str],'FontSize',12,'FontName','Arial','Color',line_colors(i+1,:));
end
slope_mean_str = sprintf('%.3f',slope_mean);
if abs(t_value_mean) > t_crit; slope_mean_str = [slope_mean_str,'*']; end
y_pos_mean = y_pred(find(x_range >= x_text,1));
text(x_text,y_pos_mean+(i+2)*text_offset,['Mean Slope = ',slope_mean_str], ...
    'FontSize',12,'FontName','Arial','Color',line_colors(1,:));
box on;

% Subplot (d)
h5 = subplot(2,5,[7 8]);
a = [imp_all(:,1) imp_all(:,2) imp_all(:,4) imp_all(:,5) imp_all(:,3)];
a = a./sum(a,2);
imp = nanmean(a);
imp_s = nanstd(a);
h = barh(imp,'FaceColor',[0.35 0.6 0.4],'EdgeColor','k','LineWidth',0.5,'BarWidth',0.7);
hold on;
h = errorbar(imp,1:5,imp_s,'horizontal','k','LineWidth',1.5,'CapSize',3);
h.LineStyle = 'none';
set(gca,'YTick',1:5,'YTickLabel',{'C4','Z_P','VPD','Rs','T'}, ...
    'FontSize',11,'FontName','Helvetica','LineWidth',0.5);
xlim([0 0.4]); ylim([0.5 5.5]);
title('(d)','FontSize',14,'FontName','Helvetica','FontWeight','bold');
xlabel('Relative Importance','FontSize',14,'FontName','Arial','FontWeight','bold');
box on;

% Subplot (e)
h6 = subplot(2,5,9);
histogram(slopes_boot(:,2),100,'Normalization','probability', ...
    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.7 0.7 0.7], ...
    'LineWidth',0.3,'FaceAlpha',0.9);
hold on;
xline(perc_95,'-','Color',[0.2 0.4 0.6],'LineWidth',1.5,'Alpha',0.9, ...
    'Label',sprintf('P95 = %.3f',perc_95),'LabelVerticalAlignment','top', ...
    'LabelHorizontalAlignment','left','FontSize',11,'FontName','Helvetica');
xline(0,'--','Color',[0.5 0.3 0.2],'LineWidth',1.5,'Alpha',0.7);
set(gca,'FontSize',11,'FontName','Helvetica','LineWidth',0.5, ...
    'XGrid','on','YGrid','on','GridColor',[0.9 0.9 0.9],'GridAlpha',0.7);
xlabel('Bootstrapped Slopes','FontSize',12,'FontName','Helvetica','FontWeight','bold');
ylabel('Probability','FontSize',12,'FontName','Helvetica','FontWeight','bold');
title('(e)','FontSize',14,'FontName','Helvetica','FontWeight','bold');
xlim([-0.015 0]); ylim([0 0.045]);
set(gca,'XTick',-0.015:0.005:0,'YTick',0:0.01:0.04);
box on; hold off;

w1 = -0.02; hei = 0.55; wei = 0.52; hight = 0.375; sp = -0.04;
set(h1,'Position',[w1 hei wei hight]);
set(h2,'Position',[w1+wei+sp hei wei hight]);
set(c1,'Position',[w1+wei*0.18 hei wei*0.7 0.02]);
set(c2,'Position',[w1+wei*1.18+sp hei wei*0.7 0.02]);
w1 = 0.075; hei = 0.1; wei = 0.22; hight = 0.325;
set(h3,'Position',[w1 hei wei*1.2 hight]);
set(h5,'Position',[w1+1.15*(wei+0.05) hei wei*1.2 hight]);
set(h6,'Position',[w1+2.5*(wei+0.04) hei wei*1.15 hight]);
set(gcf,'Renderer','painters');

filename = fullfile(figure_dir,'Mar_Figure_2_TWS_droughts.png');
exportgraphics(fig,filename,'Resolution',300);
fprintf('Created: %s\n',filename);

