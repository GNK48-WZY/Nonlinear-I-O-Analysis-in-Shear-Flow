clear all

Re_list = logspace(log10(200), log10(2000), 16);  
delta_list = logspace(-6, 0, 400);                
forcing_bi = readmatrix("forcing_bi.csv");


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


