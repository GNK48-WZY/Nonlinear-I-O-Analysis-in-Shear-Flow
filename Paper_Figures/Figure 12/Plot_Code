Re_list = readmatrix("Re_list.csv");
delta_list = readmatrix("delta_list.dat");
Gamma_theorem_e_6 = readmatrix("Gamma_theorem_e_6.csv");
u_upper_bound_e_6 = readmatrix("u_upper_bound_e_6.csv");

figure('Position', [100 100 800 600], 'Color', 'w');
contourf(delta_list, Re_list, log10(Gamma_theorem_e_6), 20, 'LineStyle', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(parula);
cbar = colorbar;
cbar.Label.String = '';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 12;
xlabel('$\delta$', 'Interpreter', 'latex', 'FontSize', 25);
ylabel('Re', 'Interpreter', 'latex', 'FontSize', 25);
title('', 'Interpreter', 'latex', 'FontSize', 16);
y = [200, 500, 1000, 1500, 2000];
ylim([200, max(Re_list)]); % Set lower y-axis limit to 200
set(gca, 'XTick', 10.^(-6:2:0), 'YTick', y);
grid on;
set(gca, 'FontSize', 18, 'LineWidth', 1.2,'ColorScale','log');
caxis([9, 16]);

disp('Figure 12a');


figure('Position', [100 100 800 600], 'Color', 'w');
contourf(delta_list, Re_list, log10(u_upper_bound_e_6), 20, 'LineStyle', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(parula);
cbar = colorbar;
cbar.Label.String = '';
cbar.Label.FontSize = 12;
xlabel('$\delta$', 'Interpreter', 'latex', 'FontSize', 25);
ylabel('Re', 'Interpreter', 'latex', 'FontSize', 25);
title('', 'Interpreter', 'latex', 'FontSize', 16);
y = [200, 500, 1000, 1500, 2000];
ylim([200, max(Re_list)]); % Set lower y-axis limit to 200
set(gca, 'XTick', 10.^(-6:2:0), 'YTick', y);
grid on;
set(gca, 'FontSize', 18, 'LineWidth', 1.2,'ColorScale','log');
caxis([-20, -11]);

disp('Figure 12b');
