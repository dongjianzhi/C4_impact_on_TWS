% Clear workspace and close figures
clear; close all;

% Read the data from three CSV files
data_vpd = readtable('/Users/jianzhidong/Documents/Research/C4transition/C4sensitivity/VPD with interval.csv');
data_T = readtable('/Users/jianzhidong/Documents/Research/C4transition/C4sensitivity/Temperature with interval.csv');
data_S = readtable('/Users/jianzhidong/Documents/Research/C4transition/C4sensitivity/SoilWaterPotential with interval.csv');

data_vpd.VPD = data_vpd.VPD*100;
% Swap Resistance_C3 and Resistance_C3_rangeL for Temperature data
dd1 = data_T.Resistance_C3;
dd2 = data_T.Resistance_C3_rangeL;
data_T.Resistance_C3 = dd2;
data_T.Resistance_C3_rangeL = dd1;

% Define Okabe-Ito color palette
c3_color = [0.000, 0.447, 0.741]; % Blue for C3
c4_color = [0.850, 0.325, 0.098]; % Orange for C4

% Create figure with Nature-compliant settings
close all
fig = figure('Color', 'white', 'Renderer', 'painters', ...
       'Units', 'inches', 'Position', [1 2 13.2 5.5]);

% Subplot 1: Resistance vs VPD
h1 = subplot(1,3,1);
hold on;
fill([data_vpd.VPD; flipud(data_vpd.VPD)], ...
     [data_vpd.Resistance_C3_rangeH; flipud(data_vpd.Resistance_C3_rangeL)], ...
     c3_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(data_vpd.VPD, data_vpd.Resistance_C3, '-', 'Color', c3_color, 'LineWidth', 2);
fill([data_vpd.VPD; flipud(data_vpd.VPD)], ...
     [data_vpd.Resistance_C4_rangeH; flipud(data_vpd.Resistance_C4_rangeL)], ...
     c4_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(data_vpd.VPD, data_vpd.Resistance_C4, '--', 'Color', c4_color, 'LineWidth', 2);
hold off;
xlabel('VPD (kPa)', 'FontSize', 10, 'FontName', 'Helvetica');
ylabel('Resistance (s m^{-1})', 'FontSize', 10, 'FontName', 'Helvetica');
title('(a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
      'Units', 'normalized');
set(gca, 'FontSize', 14, 'FontName', 'Helvetica', 'LineWidth', 1, ...
         'Box', 'off', 'TickDir', 'out', 'GridLineStyle', '--');
grid on;box on; 

% Subplot 2: Resistance vs Temperature
h2 = subplot(1,3,2);
hold on;
fill([data_T.Temp; flipud(data_T.Temp)], ...
     [data_T.Resistance_C3_rangeH; flipud(data_T.Resistance_C3_rangeL)], ...
     c3_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ll1 = plot(data_T.Temp, data_T.Resistance_C3, '-', 'Color', c3_color, 'LineWidth', 2);
fill([data_T.Temp; flipud(data_T.Temp)], ...
     [data_T.Resistance_C4_rangeH; flipud(data_T.Resistance_C4_rangeL)], ...
     c4_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ll2 = plot(data_T.Temp, data_T.Resistance_C4, '--', 'Color', c4_color, 'LineWidth', 2);
hold off;
xlabel('Temperature (°C)', 'FontSize', 10, 'FontName', 'Helvetica');
title('(b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
      'Units', 'normalized');
ylabel('Resistance (s m^{-1})', 'FontSize', 10, 'FontName', 'Helvetica');
set(gca, 'FontSize', 14, 'FontName', 'Helvetica', 'LineWidth', 1, ...
         'Box', 'off', 'TickDir', 'out', 'GridLineStyle', '--');
grid on; box on

% Subplot 3: Resistance vs Soil Water Potential
h3 = subplot(1,3,3);
hold on;
fill([data_S.SoilWaterPotential; flipud(data_S.SoilWaterPotential)], ...
     [data_S.Resistance_C3_rangeH; flipud(data_S.Resistance_C3_rangeL)], ...
     c3_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(data_S.SoilWaterPotential, data_S.Resistance_C3, '-', 'Color', c3_color, 'LineWidth', 2);
fill([data_S.SoilWaterPotential; flipud(data_S.SoilWaterPotential)], ...
     [data_S.Resistance_C4_rangeH; flipud(data_S.Resistance_C4_rangeL)], ...
     c4_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(data_S.SoilWaterPotential, data_S.Resistance_C4, '--', 'Color', c4_color, 'LineWidth', 2);
hold off;
xlabel('Soil Water Potential (MPa)', 'FontSize', 10, 'FontName', 'Helvetica');
title('(c)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
      'Units', 'normalized');
ylabel('Resistance (s m^{-1})', 'FontSize', 10, 'FontName', 'Helvetica');
set(gca, 'FontSize', 14, 'FontName', 'Helvetica', 'LineWidth', 1, ...
         'Box', 'off', 'TickDir', 'out', 'GridLineStyle', '--');
grid on; box on;

% Add shared legend above subplots
lgd = legend(subplot(1,3,3), [ll1,ll2],'C3', 'C4', 'Location', 'northoutside', ...
             'Orientation', 'horizontal', 'FontSize', 14, 'FontName', 'Helvetica');

% Add main title
% sgtitle('C3 and C4 Resistance Responses to Environmental Factors', ...
%         'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');

sp = 0.075; wei = 0.25; high1 = 0.18;
set(h1,'position',[0.075 high1 wei 0.7])
set(h2,'position',[0.075+wei+sp high1 wei 0.7])
set(h3,'position',[0.075+(wei+sp)*2 high1 wei 0.7])
lgd.Position = [0.1 0.81 0.1 0.05]; % Centered above

% Export as high-quality PDF
filename='./Figures/Figure3_stomatal_conductance.png'; 
exportgraphics(fig, filename, 'Resolution', 300);








