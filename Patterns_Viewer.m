close all;
clc;
clear;
warning('off', 'all');

Neuron_width = 4;
Neuron_height = 5;

path1 = '/home/makovkin_s/HH_Recogn_Rev_4/HH_Recogn_2023_03_15_4x5_50x50_r_P';

FontSize = 8;
FontSize2 = 24;

filename1 = 'N_A.txt';
Data_N_A = importdata( [ path1 '/' filename1 ] );
disp('Data_N_A load complete');

filename1 = 'tau.txt';
Data_tau = importdata( [ path1 '/' filename1 ] );
disp('Data_tau load complete');

filename1 = 'E_syn_matrix.txt';
Data_E_syn_matrix = importdata( [ path1 '/' filename1 ] );
disp('Data_tau load complete');

filename1 = 'g_syn.txt';
Data_g_syn = importdata( [ path1 '/' filename1 ] );
disp('g_syn load complete');

g_syn_IN = Data_g_syn(end-Neuron_width*Neuron_height, 1);
g_syn_IN_string = num2str(g_syn_IN, '%.4f');

filename1 = 'I_app.txt';
Data_I_app = importdata( [ path1 '/' filename1 ] );
disp('I_app load complete');

I_app_CN = Data_I_app(end, 1);
I_app_CN_string = num2str(I_app_CN, '%.4f');

fig_tau_N_A_E_syn = figure('units', 'normalized', 'outerposition', [0 0 1.0 0.66], 'visible', 'on');
ax131 = subplot(1, 3, 1);
%pcolor(Data_tau');
imshow(Data_tau);
%axis ij;
axis('on', 'image');
%axis equal;
daspect([1 1 1])
%%axis off;
colormap(ax131, hot);
caxis([0 30]);
cbh1 = colorbar('v');
title(cbh1, 'ms');
%set(cbh, 'YTick', [0:1:1]);
%set(gca, 'XTick', [1 size(Data_A_N, 1)/2 size(Data_A_N, 1)]);
%set(gca, 'YTick', [1 size(Data_A_N, 2)/2 size(Data_A_N, 2)]);
set(gca, 'XTick', [1, [10: 10: size(Data_N_A, 1)]]);
set(gca, 'YTick', [1, [10: 10: size(Data_N_A, 2)]]);
xtickangle(90);
ax = gca;
ax.XAxisLocation = 'top';
title('Neurons \tau_{lag delay} matrix');
set(gca, 'FontSize', FontSize);

ax132 = subplot(1, 3, 2);
imshow(Data_N_A);
%axis ij;
axis('on', 'image');
%axis equal;
daspect([1 1 1])
%%axis off;
colormap(ax132, gray(2));
caxis([0 1]);
cbh2 = colorbar('v');
title(cbh2, '\xi_i');
set(cbh2, 'YTick', [0:1:1]);
%set(gca, 'XTick', [1 size(Data_A_N, 1)/2 size(Data_A_N, 1)]);
%set(gca, 'YTick', [1 size(Data_A_N, 2)/2 size(Data_A_N, 2)]);
set(gca, 'XTick', [1, [10: 10: size(Data_N_A, 1)]]);
set(gca, 'YTick', [1, [10: 10: size(Data_N_A, 2)]]);
xtickangle(90);
ax = gca;
ax.XAxisLocation = 'top';
title('Neurons adjacency matrix');
set(gca, 'FontSize', FontSize);

ax133 = subplot(1, 3, 3);
imshow(Data_E_syn_matrix);
%axis ij;
axis('on', 'image');
%axis equal;
daspect([1 1 1])
%%axis off;
colormap(ax133, hot);
caxis([-90 0]);
cbh3 = colorbar('v');
title(cbh3, 'mV');
set(cbh3, 'YTick', [-90:10:0]);
%set(gca, 'XTick', [1 size(Data_A_N, 1)/2 size(Data_A_N, 1)]);
%set(gca, 'YTick', [1 size(Data_A_N, 2)/2 size(Data_A_N, 2)]);
set(gca, 'XTick', [1, [10: 10: size(Data_N_A, 1)]]);
set(gca, 'YTick', [1, [10: 10: size(Data_N_A, 2)]]);
xtickangle(90);
ax = gca;
ax.XAxisLocation = 'top';
title('E_{syn} value matrix');
set(gca, 'FontSize', FontSize);

saveas(gcf, 'fig_tau_N_A_E_syn.png');
%saveas(gcf, 'fig_tau_N_A_E_syn', 'epsc');
close(fig_tau_N_A_E_syn);
disp('Picture fig_tau_N_A_E_syn save');

filename1 = 'SN_input_array.txt';
Data_SN_input_array = importdata( [ path1 '/' filename1 ] );
disp('Data_SN_input_array load complete');

filename1 = 'Patterns.txt';
Data_patterns = importdata( [ path1 '/' filename1 ] );
disp('Data_patterns load complete');

Layer_size = Neuron_width * Neuron_height;
Number_of_patterns = size(Data_patterns, 1) / Layer_size;
Patterns = reshape(Data_patterns, [Layer_size, Number_of_patterns]);

Data_Pattern_1x5_12 = Patterns(:, 1);
Data_Pattern_1x5_3 = Patterns(:, 2);
Data_Pattern_1x5_5 = Patterns(:, 3);

xi_12_ij = Data_Pattern_1x5_12 * Data_Pattern_1x5_12';
xi_3_ij = Data_Pattern_1x5_3 * Data_Pattern_1x5_3';
xi_5_ij = Data_Pattern_1x5_5 * Data_Pattern_1x5_5';

xi_all = xi_12_ij + xi_3_ij + xi_5_ij;

xi_all_signum = xi_all;
xi_all_signum(xi_all_signum >= 0) = 1;
xi_all_signum(xi_all_signum < 0) = -1;

a2 = 1;
b2 = 0;
E_syn_IN_matlab_binary_2 = a2 * xi_all_signum + b2;

Correct_recognized_pattern = Data_SN_input_array' * E_syn_IN_matlab_binary_2;
Correct_binary_pattern = Correct_recognized_pattern;
Correct_binary_pattern(Correct_binary_pattern >= 0) = 1;
Correct_binary_pattern(Correct_binary_pattern < 0) = -1;
%{
fig_E_syn = figure('units', 'normalized', 'outerposition', [0 0 1.0 0.5], 'visible', 'on');
ax121 = subplot(1, 2, 1);
imshow(E_syn_IN_matlab_binary_2);
%axis ij;
axis('on', 'image');
%axis equal;
daspect([1 1 1])
%%axis off;
colormap(ax121, gray(2));
caxis([-1 1]);
cbh3 = colorbar('v');
title(cbh3, '\xi_i');
set(cbh3, 'YTick', [-1:1:1]);
set(gca, 'XTick', [1, [5: 5: size(E_syn_IN_matlab_binary_2, 1)]]);
set(gca, 'YTick', [1, [5: 5: size(E_syn_IN_matlab_binary_2, 2)]]);
xtickangle(90);
ax = gca;
ax.XAxisLocation = 'top';
title('H_{ij}');
set(gca, 'FontSize', FontSize);

saveas(gcf, 'fig_E_syn.png');
%saveas(gcf, 'fig_E_syn', 'epsc');
close(fig_E_syn);
disp('Picture fig_E_syn save');
%}
filename1 = 'results_V.txt';
Data = importdata( [ path1 '/' filename1 ] );
time_V = Data(round(0.25*end):end, 1);
V = Data(round(0.25*end):end, 2:end);
%dt = time_V(2) - time_V(1);
disp('results_V load complete');

time_Debug_begin = 0.0;
d_time_Debug = 0.5;
time_Debug_end = 10.0;
xtckfrmt = '%.1f';

Neuron_RN = 1;
%CN_num = 1 + 2 * Layer_size;
%CN_num = 1 + 3 * Layer_size;

%Begin Freqs IN
omega_cur = 0;
omega_total = [];
time_total = [];
freqs_full = [];
%yyaxis left;
regular_V_spiking_counter = 0;
for i = 1 + Layer_size + 1: 1: 1 + 2 * Layer_size
%for i = 5: 1: 7
%i = 85;
    V_p = V(:, i);
    [pks, locs] = findpeaks(V_p, time_V, 'MinPeakHeight', -10, 'MinPeakDistance', 0.001);

    for pks_counter = 1:size(pks, 1)
        pks(pks_counter) = i;
    end

    if isempty([pks, locs]) || size(pks, 1) == 1
        omega_cur = 0;
        omega_total = [omega_total omega_cur'];
        time_locs = locs;
        time_locs = time_locs(2:end);
        %time_locs = time_locs(find(omega_cur'<=60));
        time_total = [time_total time_locs'];
        %plot(omega_cur, i, 'd', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        %hold on;
    end
    if ~isempty([pks, locs]) && size(pks, 1) > 1
        regular_V_spiking_counter = regular_V_spiking_counter + 1;
        omega_cur = 1./diff(locs);
        %omega_cur = omega_cur(find(omega_cur'<=60));
        omega_total = [omega_total omega_cur'];
        time_locs = locs;
        time_locs = time_locs(2:end);
        %time_locs = time_locs(find(omega_cur'<=60));
        %time_total = [time_total time_locs'];
        %plot(omega_cur, pks(2:end), 'd', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        %hold on;
    end

	%{
	if ~isempty([pks, locs]) && size(pks, 1) > 1
        plot(time_locs, omega_cur,  '-', 'LineWidth', 1);
        hold on;
    end
	%}
    
	nu = omega_cur';
    nu_time = time_locs';
    time = time_V';
    
    if isempty(nu_time)
        Nu_new = zeros(size(time));    
    else
        if length(nu)<2
            nu = repmat(nu, 1, 2);
            nu_time = [nu_time nu_time + dt];
        end
        Nu_new = interp1(nu_time, nu, time, 'previous');
        Nu_new(1: find(Nu_new > 0, 1) - 1) = 0;
        Nu_new(find(Nu_new > 0, 1, 'last') + 1 : end) = nu(end);
    end

    %plot(time, Nu_new,  '-', 'LineWidth', 1);
	%hold on;
    
    freqs_full = [freqs_full Nu_new'];

end
% End Freqs IN

% Begin Phase IN-RN
%omega_cur = 0;
%omega_total = [];
%yyaxis left;
phi_new_IN_RN_all_data = [];
regular_V_spiking_counter = 0;
%for i = 1:size(V, 2)
i = Neuron_RN;
V_p = V(:, i);
[pks, locs] = findpeaks(V_p, 'MinPeakHeight', -10);
%end
time_locs_RN = time_V(locs);
period_RN = mean(diff(time_locs_RN), 1);

for i = 1 + Layer_size + 1: 1: 1 + 2 * Layer_size
    %i = CN_num;
    V_p = V(:, i);
    [pks, locs] = findpeaks(V_p, 'MinPeakHeight', -10);
    %end
    time_locs_IN = time_V(locs);
    
    if isempty(time_locs_IN)
        nu_IN = 0;
        delta_phi_IN_RN = zeros(size(time));
        phi_IN_RN_0_2pi_Mean = 0;
        phi_IN_RN_0_2pi_Disp = 0;
    else
        period_IN = mean(diff(time_locs_IN), 1);
        nu_IN = 1 / period_IN;
        
        if time_locs_IN(1) < time_locs_RN(1)
            time_locs_IN(1) = [];
        end
        
        min_time_size_IN_RN = min(size(time_locs_IN, 1), size(time_locs_RN, 1));
        time_locs_IN = time_locs_IN(1:min_time_size_IN_RN);
        time_locs_RN = time_locs_RN(1:min_time_size_IN_RN);
        
        phi_IN_RN = 2*pi*(time_locs_IN - time_locs_RN) / period_RN;
        phi_IN_RN_0_2pi = mod(phi_IN_RN, 2*pi);
        
        phi_IN_RN_0_2pi_Mean = mean(phi_IN_RN_0_2pi, 1);
        phi_IN_RN_0_2pi_Disp = sqrt(var(phi_IN_RN_0_2pi, 1));
        
        nu = phi_IN_RN_0_2pi';
        nu_time = time_locs_IN';
        time = time_V';
        
        if isempty(nu_time)
            delta_phi_IN_RN = zeros(size(time));
        else
            delta_phi_IN_RN = interp1(nu_time, nu, time, 'previous');
            delta_phi_IN_RN(1: find(delta_phi_IN_RN > 0, 1) - 1) = 0;
            delta_phi_IN_RN(find(delta_phi_IN_RN > 0, 1, 'last') + 1 : end) = nu(end);
        end
    end
    
    phi_new_IN_RN_all_data = [phi_new_IN_RN_all_data delta_phi_IN_RN'];
end

d_time_Debug2 = 1.0;
xtckfrmt2 = '%.0f';

fig_IN_Phases = figure('units', 'normalized', 'outerposition', [0 0 0.7 1.0], 'visible', 'on');

plot(time, phi_new_IN_RN_all_data,  '-', 'LineWidth', 2);
grid on;
box on;
hold on;
xlim([time_Debug_begin time_Debug_end]);
ylim([0 2*pi]);
ax = gca;
ax.XTick = [time_Debug_begin: d_time_Debug2: time_Debug_end];
xtickformat(xtckfrmt2);
ax.YTick = [0: pi/2: 2*pi];
%ytickformat('%.2f');
yticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
%legend('<\nu_i>');
xlabel('Time, s', 'FontSize', FontSize);
ylabel('\Delta \phi_{V_{IN}-V_{RN}}', 'FontSize', FontSize);
title('\Delta \phi between V_{IN} and V_{RN}', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize2);

saveas(gcf, ['fig_IN_Phases_' I_app_CN_string '_' g_syn_IN_string '.png']);
%saveas(gcf, 'fig_IN_Phases', 'epsc');
close(fig_IN_Phases);
disp('Picture fig_IN_Phases save');

% End Phase IN-RN

% Begin Phase CN-RN
%omega_cur = 0;
%omega_total = [];
%yyaxis left;
phi_new_CN_RN_all_data = [];
regular_V_spiking_counter = 0;
%for i = 1:size(V, 2)
i = Neuron_RN;
V_p = V(:, i);
[pks, locs] = findpeaks(V_p, 'MinPeakHeight', -10);
%end
time_locs_RN = time_V(locs);
period_RN = mean(diff(time_locs_RN), 1);

for i = 1 + 2 * Layer_size + 1: 1: 1 + 3 * Layer_size
    %i = CN_num;
    V_p = V(:, i);
    [pks, locs] = findpeaks(V_p, 'MinPeakHeight', -10);
    %end
    time_locs_CN = time_V(locs);
    
    if isempty(time_locs_CN)
        nu_CN = 0;
        delta_phi_CN_RN = zeros(size(time));
        phi_CN_RN_0_2pi_Mean = 0;
        phi_CN_RN_0_2pi_Disp = 0;
    else
        period_CN = mean(diff(time_locs_CN), 1);
        nu_CN = 1 / period_CN;
        
        if time_locs_CN(1) < time_locs_RN(1)
            time_locs_CN(1) = [];
        end
        
        min_time_size_CN_RN = min(size(time_locs_CN, 1), size(time_locs_RN, 1));
        time_locs_CN = time_locs_CN(1:min_time_size_CN_RN);
        time_locs_RN = time_locs_RN(1:min_time_size_CN_RN);
        
        phi_CN_RN = 2*pi*(time_locs_CN - time_locs_RN) / period_RN;
        phi_CN_RN_0_2pi = mod(phi_CN_RN, 2*pi);
        
        phi_CN_RN_0_2pi_Mean = mean(phi_CN_RN_0_2pi, 1);
        phi_CN_RN_0_2pi_Disp = sqrt(var(phi_CN_RN_0_2pi, 1));
        
        nu = phi_CN_RN_0_2pi';
        nu_time = time_locs_CN';
        time = time_V';
        
        if isempty(nu_time)
            delta_phi_CN_RN = zeros(size(time));
        else
            delta_phi_CN_RN = interp1(nu_time, nu, time, 'previous');
            delta_phi_CN_RN(1: find(delta_phi_CN_RN > 0, 1) - 1) = 0;
            delta_phi_CN_RN(find(delta_phi_CN_RN > 0, 1, 'last') + 1 : end) = nu(end);
        end
    end
    
    phi_new_CN_RN_all_data = [phi_new_CN_RN_all_data delta_phi_CN_RN'];
end

fig_CN_Phases = figure('units', 'normalized', 'outerposition', [0 0 0.7 1.0], 'visible', 'on');

plot(time, phi_new_CN_RN_all_data,  '-', 'LineWidth', 2);
grid on;
box on;
hold on;
xlim([time_Debug_begin time_Debug_end]);
ylim([0 2*pi]);
ax = gca;
ax.XTick = [time_Debug_begin: d_time_Debug2: time_Debug_end];
xtickformat(xtckfrmt2);
ax.YTick = [0: pi/2: 2*pi];
%ytickformat('%.2f');
yticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
%legend('<\nu_i>');
xlabel('Time, s', 'FontSize', FontSize);
ylabel('\Delta \phi_{V_{CN}-V_{RN}}', 'FontSize', FontSize);
title('\Delta \phi between V_{CN} and V_{RN}', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize2);

saveas(gcf, ['fig_CN_Phases_' I_app_CN_string '_' g_syn_IN_string '.png']);
%saveas(gcf, 'fig_CN_Phases', 'epsc');
close(fig_CN_Phases);
disp('Picture fig_CN_Phases save');

% End Phase CN-RN

Result_binary_pattern = phi_new_CN_RN_all_data(end, :);
Result_binary_pattern(Result_binary_pattern >= 0 & Result_binary_pattern < pi) = -1;
Result_binary_pattern(Result_binary_pattern >= pi & Result_binary_pattern < 2 * pi) = 1;

%Result_binary_pattern == Correct_binary_pattern

fig_SN_IN_CN_Patterns = figure('units', 'normalized', 'outerposition', [0 0 1.0 1.0], 'visible', 'on');
subplot(3, 5, 1);
image_V = reshape(Data_SN_input_array, [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
SNh1 = colorbar;
set(SNh1, 'YTick', [-1:1:1]);
title(SNh1, '\xi_i');
title(['SN, \xi_i. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 2);
image_VP = reshape(phi_new_IN_RN_all_data(end, :), [Neuron_width, Neuron_height]);
imagesc(image_VP');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, hsv);
caxis([0 2*pi]);
CNh1 = colorbar;
%set(h1, 'YTick', [0: pi/2: 2*pi]);
set(CNh1,'YTick',[0: pi/2: 2*pi], 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title(CNh1, 'rad');
title(['IN, \Delta\phi_{IN-RN}. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 3);
image_VP = reshape(phi_new_CN_RN_all_data(end, :), [Neuron_width, Neuron_height]);
imagesc(image_VP');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, hsv);
caxis([0 2*pi]);
CNh1 = colorbar;
%set(h1, 'YTick', [0: pi/2: 2*pi]);
set(CNh1,'YTick',[0: pi/2: 2*pi], 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title(CNh1, 'rad');
title(['CN, \Delta\phi_{CN-RN}. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 4);
image_V = reshape(Result_binary_pattern, [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
SNh4 = colorbar;
set(SNh4, 'YTick', [-1:1:1]);
title(SNh4, '\xi_i');
title(['Result Binary Pattern. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 5);
image_V = reshape(Correct_binary_pattern, [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
SNh4 = colorbar;
set(SNh4, 'YTick', [-1:1:1]);
title(SNh4, '\xi_i');
title(['Corrected recognized pattern. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 6);
image_V = reshape(Patterns(:, 1), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph1 = colorbar;
set(ph1, 'YTick', [-1:1:1]);
title(ph1, '\xi_i');
title(['Pattern #1. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 7);
image_V = reshape(Patterns(:, 2), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph2 = colorbar;
set(ph2, 'YTick', [-1:1:1]);
title(ph2, '\xi_i');
title(['Pattern #2. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 8);
image_V = reshape(Patterns(:, 3), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph3 = colorbar;
set(ph3, 'YTick', [-1:1:1]);
title(ph3, '\xi_i');
title(['Pattern #3. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);
%{
subplot(3, 5, 9);
image_V = reshape(Patterns(:, 4), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph3 = colorbar;
set(ph3, 'YTick', [-1:1:1]);
title(ph3, '\xi_i');
title(['Pattern #4. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 10);
image_V = reshape(Patterns(:, 5), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph3 = colorbar;
set(ph3, 'YTick', [-1:1:1]);
title(ph3, '\xi_i');
title(['Pattern #5. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 11);
image_V = reshape(Patterns(:, 6), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph3 = colorbar;
set(ph3, 'YTick', [-1:1:1]);
title(ph3, '\xi_i');
title(['Pattern #6. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);

subplot(3, 5, 12);
image_V = reshape(Patterns(:, 7), [Neuron_width, Neuron_height]);
imagesc(image_V');
axis equal tight;
set(gca,'xtick',[],'ytick',[])
box on;
colormap(gca, gray(2));
caxis([-1 1]);
ph3 = colorbar;
set(ph3, 'YTick', [-1:1:1]);
title(ph3, '\xi_i');
title(['Pattern #7. ' num2str(Neuron_width) 'x' num2str(Neuron_height)]);
set(gca, 'FontSize', FontSize);
%}
%saveas(gcf, 'fig_SN_IN_CN_Patterns.png');
saveas(gcf, ['fig_SN_IN_CN_Patterns_' I_app_CN_string '_' g_syn_IN_string '.png']);
%saveas(gcf, 'fig_SN_IN_CN_Patterns', 'epsc');
close(fig_SN_IN_CN_Patterns);
disp('Picture fig_SN_IN_CN_Patterns save');

Comparisson_flag = all(Result_binary_pattern == Correct_binary_pattern);

return

Data_I_app_CN_g_syn_IN_Comparisson_patterns = [I_app_CN, g_syn_IN, Comparisson_flag];
formatSpec = '%2.4f\t%2.4f\t%2.4f\n';
fileID = fopen('Data_I_app_CN_g_syn_IN_Comparisson_patterns.txt', 'A');
fprintf(fileID, formatSpec, Data_I_app_CN_g_syn_IN_Comparisson_patterns);
fclose(fileID);

return