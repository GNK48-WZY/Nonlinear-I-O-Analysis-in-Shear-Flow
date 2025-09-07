% Figure 8a
clear all
M = [5	10	20	30	40	50	60	70	80	90];
forcing_list = readmatrix("forcing_list.csv");
figure;
scatter(M, forcing_list);
xlabel('$M$', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('$f_\Pi$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'FontSize', 12);
grid on;



% Figure 8b
clear all

Re_list = logspace(log10(200), log10(2000), 16);  
delta_list = logspace(-6, 0, 40);                
forcing_bi = readmatrix("forcing_bi_M_40_delta_40.csv");


figure('Position', [100 100 800 600], 'Color', 'w');
contourf(delta_list, Re_list, log10(forcing_bi), 20, 'LineStyle', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(parula);
cbar = colorbar;
cbar.Label.String = '';
cbar.Label.FontSize = 12;
xlabel('$\delta$', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('Re', 'Interpreter', 'latex', 'FontSize', 28);
title('', 'Interpreter', 'latex', 'FontSize', 18);
y = [200, 500, 1000, 1500, 2000];
ylim([200, 2000])
set(gca, 'XTick', 10.^(-6:2:0), 'YTick', y);
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'ColorScale','log');


