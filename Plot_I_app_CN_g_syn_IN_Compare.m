close all;
clc;
clear;

FontSize = 14;

filename4 = 'Data_I_app_CN_g_syn_IN_Comparisson_patterns.txt';
path4 = '/common/home/makovkin_s/HH_Recogn_Rev_4/HH_Recogn_2023_03_15_4x5_50x50_r_P';

Data = importdata( [ path4 '/' filename4 ] );
I_app_CN = Data(:, 1);
g_syn_IN = Data(:, 2);
Compare_flag  = Data(:, 3);

start_I_app_SN = 0.00;
d_I_app_SN     = 0.024;
end_I_app_SN   = 1.20;

start_g_syn_IN = 0.00;
d_g_syn_IN     = 0.0064;
end_g_syn_IN   = 0.32;

last_g_syn_IN = g_syn_IN(end);
if last_g_syn_IN < end_g_syn_IN
    zero_array = zeros( round( ((end_g_syn_IN - start_g_syn_IN) / d_g_syn_IN + 1) * ( (end_I_app_SN - start_I_app_SN) / d_I_app_SN + 1) ) - round( size(g_syn_IN, 1) ), 1);
    I_app_CN = [I_app_CN' zero_array']';
    g_syn_IN = [g_syn_IN' zero_array']';
    Compare_flag = [Compare_flag' zero_array']';
end

[X, Y] = meshgrid(start_g_syn_IN : d_g_syn_IN : end_g_syn_IN, start_I_app_SN : d_I_app_SN : end_I_app_SN);
reshape_Compare_flag = reshape(Compare_flag, round( (end_g_syn_IN - start_g_syn_IN) / d_g_syn_IN) + 1, round( (end_I_app_SN - start_I_app_SN) / d_I_app_SN) + 1)';

fig_results_I_app_CN_g_syn_IN_Compare_flag = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'visible', 'on');
pcolor(X, Y, reshape_Compare_flag);
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
box on;
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(gray(2));
cbh = colorbar('v');
title(cbh, 'Compare flag');
caxis([0.0 1.0]);
set(cbh, 'YTick', [0: 1: 1]);
ylim([0 1.2]);
cbh.Ruler.TickLabelFormat = '%.0f';
%set(cbh, 'YTick', [0.0:0.2:2.0]);
set(gca, 'YTick', [0: 0.1: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.02: 0.32]);
xtickformat('%.2f');
ylabel('I_{app CN}, \muA', 'FontSize', FontSize);
xlabel('g_{syn IN}, \muA/cm^2', 'FontSize', FontSize);
title(['CN recognition regions'], 'FontSize', FontSize);
%title('I_{app RN} = 1.05', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);

return

saveas(gcf, 'fig_results_I_app_CN_g_syn_IN_Compare_flag.png');
%saveas(gcf, 'fig_results_I_app_CN_g_syn_IN_Compare_flag', 'epsc');
close(fig_results_I_app_CN_g_syn_IN_Compare_flag);
disp('Picture fig_results_I_app_CN_g_syn_IN_Compare_flag save');

return

fig_results_I_app_SN_g_syn_RN_delta_phi_SN_RN = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'visible', 'on');
pcolor(X, Y, reshape_delta_phi_SN_RN);
hold on;
%plot(0.12, 0.8, '.r', 'MarkerSize', 20);
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(hsv);
h1 = colorbar('v');
title(h1, '\Delta \phi_{SN-RN}');
caxis([0.0 2*pi]);
set(h1, 'YTick', [0: pi/4: 2*pi], 'YTickLabel', {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})
h1.Ruler.TickLabelFormat = '%.0f';
ylim([0 1.2]);
set(gca, 'YTick', [0: 0.1: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.02: 0.32]);
xtickformat('%.2f');
ylabel('I_{app SN}, \muA', 'FontSize', FontSize);
xlabel('g_{syn RN}, \muA/cm^2', 'FontSize', FontSize);
title(['SN excitation regions'], 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);

saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_delta_phi_SN_RN.png');
%saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_delta_phi_SN_RN', 'epsc');
close(fig_results_I_app_SN_g_syn_RN_delta_phi_SN_RN);
disp('Picture fig_results_I_app_SN_g_syn_RN_delta_phi_SN_RN save');

fig_results_I_app_SN_g_syn_RN_disp_phi_SN_RN = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'visible', 'on');
pcolor(X, Y, reshape_disp_phi_SN_RN);
hold on;
%plot(0.12, 0.8, '.r', 'MarkerSize', 20);
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(hot);
h2 = colorbar('v');
title(h2, '$$\sqrt{D}$$', 'interpreter', 'latex');
caxis([0.0 2.0]);
set(h2, 'YTick', [0: 0.5: 2.0]);
h2.Ruler.TickLabelFormat = '%.1f';
ylim([0 1.2]);
set(gca, 'YTick', [0: 0.1: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.02: 0.32]);
xtickformat('%.2f');
ylabel('I_{app SN}, \muA', 'FontSize', FontSize);
xlabel('g_{syn RN}, \muA/cm^2', 'FontSize', FontSize);
title(['SN excitation regions'], 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);

saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_disp_phi_SN_RN.png');
%saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_disp_phi_SN_RN', 'epsc');
close(fig_results_I_app_SN_g_syn_RN_disp_phi_SN_RN);
disp('Picture fig_results_I_app_SN_g_syn_RN_disp_phi_SN_RN save');

FontSize2 = 12;
%FontSize2 = 15;

fig_results_I_app_SN_g_syn_RN_freq_SN_all = figure('units', 'normalized', 'outerposition', [0 0 1 0.6], 'visible', 'on');
ax131 = subplot(1, 3, 1);
pcolor(X, Y, reshape_Compare_flag);
hold on;
%plot(0.30, 1.12, '.m', 'MarkerSize', 18);
%hold on;
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(ax131, hot);
cbh = colorbar('v');
title(cbh, '\nu_{SN}, s^{-1}');
caxis([0.0 40.0]);
set(cbh, 'YTick', [0: 10: 40]);
cbh.Ruler.TickLabelFormat = '%.0f';
ylim([0 1.2]);
set(gca, 'YTick', [0: 0.3: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.08: 0.32]);
xtickformat('%.2f');
ylabel('I_{app SN}, \muA', 'FontSize', FontSize2);
xlabel('g_{syn RN}, \muA/cm^2', 'FontSize', FontSize2);
title(['SN excitation regions'], 'FontSize', FontSize2);
set(gca, 'FontSize', FontSize2);

ax132 = subplot(1, 3, 2);
pcolor(X, Y, reshape_delta_phi_SN_RN);
hold on;
%plot(0.30, 1.12, '.m', 'MarkerSize', 18);
%hold on;
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(ax132, hsv);
h1 = colorbar('v');
title(h1, '\Delta \phi_{SN-RN}');
caxis([0.0 2*pi]);
set(h1, 'YTick', [0: pi/2: 2*pi], 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%set(h1, 'YTick', [0: pi/4: 2*pi], 'YTickLabel', {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})
h1.Ruler.TickLabelFormat = '%.0f';
ylim([0 1.2]);
set(gca, 'YTick', [0: 0.3: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.08: 0.32]);
xtickformat('%.2f');
%ylabel('I_{app SN}', 'FontSize', FontSize2);
xlabel('g_{syn RN}, \muA/cm^2', 'FontSize', FontSize2);
title(['SN excitation regions'], 'FontSize', FontSize2);
set(gca, 'FontSize', FontSize2);

ax133 = subplot(1, 3, 3);
pcolor(X, Y, reshape_disp_phi_SN_RN);
hold on;
%plot(0.30, 1.12, '.m', 'MarkerSize', 18);
%hold on;
%surf(X, Y, reshape_corr);
shading flat;
pbaspect([1 1 1]);
%clrmap = [1 0 0 0 1 0];
%colormap(clrmap);
colormap(ax133, hot);
h2 = colorbar('v');
title(h2, '$$\sqrt{D}$$', 'interpreter', 'latex');
caxis([0.0 2.0]);
set(h2, 'YTick', [0: 0.5: 2.0]);
h2.Ruler.TickLabelFormat = '%.1f';
ylim([0 1.2]);
set(gca, 'YTick', [0: 0.3: 1.2]);
ytickformat('%.1f');
set(gca, 'XTick', [0.00: 0.08: 0.32]);
xtickformat('%.2f');
%ylabel('I_{app SN}', 'FontSize', FontSize2);
xlabel('g_{syn RN}, \muA/cm^2', 'FontSize', FontSize2);
title(['SN excitation regions'], 'FontSize', FontSize2);
set(gca, 'FontSize', FontSize2);

saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_freq_SN_all.png');
%saveas(gcf, 'fig_results_I_app_SN_g_syn_RN_freq_SN_all', 'epsc');
close(fig_results_I_app_SN_g_syn_RN_freq_SN_all);
disp('Picture fig_results_I_app_SN_g_syn_RN_freq_SN_all save');