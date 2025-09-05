% Figure 9a
clear all

line_width = 1.5;
colors = lines(16);
line_styles = {"-","--","-.",":"};
markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
line_width = 1.5;
figure;
max_forcing_bi_M_5_delta_400 = max(forcing_bi_M_5_delta_400, [], 2); 
max_forcing_bi_M_5_delta_40 = max(forcing_bi_M_5_delta_40, [], 2);
max_forcing_bi_M_40_delta_20 = max(forcing_bi_M_40_delta_20, [], 2);
% max_forcing_bi_M_40_delta_40 = max(forcing_bi_M_40_delta_40, [], 2);
loglog(Re_list, max_forcing_bi_M_5_delta_400(:,:), line_styles{3}, 'LineWidth', line_width, 'Color', colors(3, :),'Marker', markers{1});
% Retain the current plot and add the next plot
hold on;
loglog(Re_list, max_forcing_bi_M_5_delta_40(:, :), line_styles{1}, 'LineWidth', line_width, 'Color', colors(1, :),'Marker', markers{5});
loglog(Re_list, max_forcing_bi_M_40_delta_20(:, :), line_styles{2}, 'LineWidth', line_width, 'Color', colors(2, :),'Marker', markers{2});
% loglog(Re_list, max_forcing_bi_M_40_delta_40(:, :), line_styles{4}, 'LineWidth', line_width, 'Color', colors(4, :),'Marker', markers{3});

xlabel('Re', 'FontSize', 12,'Interpreter', 'latex');
xlim([1e2, 2000]);  
legend('$f_{c,\Pi}$, 400 $\delta$, M = 5','$f_{c,\Pi}$, 40 $\delta$, M = 5','$f_{c,\Pi}$, 20 $\delta$, M = 40','Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);

hold off;

legend("Position",[0.73899,0.49841,0.14821,0.14048])



% Figure 9b
clear all

Re_list = logspace(log10(200), log10(2000), 16);  
forcing_bi = readmatrix("forcing_bi.csv");
u_upper_bound_max = readmatrix("u_upper_bound_max.csv");
u_upper_bound_max_SOS = readmatrix("u_upper_bound_max_SOS.csv");

colors = lines(16);
line_styles = {"-","--","-.",":"};
markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
line_width = 1.5;
max_forcing_bi = max(forcing_bi, [],2); 

figure;

% Plot original data with DisplayNames
loglog(Re_list, u_upper_bound_max, line_styles{3}, 'LineWidth', line_width, ...
    'Color', colors(3, :), 'Marker', markers{1}, 'DisplayName', '$f_{c,LMI}$');
hold on;
loglog(Re_list, u_upper_bound_max_SOS, line_styles{2}, 'LineWidth', line_width, ...
    'Color', colors(2, :), 'Marker', markers{3}, 'DisplayName', '$f_{c,SOS}$');
loglog(Re_list, max_forcing_bi, line_styles{4}, 'LineWidth', line_width, ...
    'Color', colors(4, :), 'Marker', markers{2}, 'DisplayName', '$f_{c,\Pi}$');

% Compute log values for fitting
log_Re = log10(Re_list(:));

% Fit for u_upper_bound_max
log_data1 = log10(u_upper_bound_max(:));
coeffs1 = polyfit(log_Re, log_data1, 1);
fit_line1 = 10.^polyval(coeffs1, log_Re);
loglog(Re_list, fit_line1, '-', 'LineWidth', line_width, 'Color', colors(3, :), ...
    'DisplayName', sprintf('Fit: $10^{%.2f} Re^{%.2f}$', coeffs1(2), coeffs1(1)));

% Fit for u_upper_bound_max_SOS
log_data2 = log10(u_upper_bound_max_SOS(:));
coeffs2 = polyfit(log_Re, log_data2, 1);
fit_line2 = 10.^polyval(coeffs2, log_Re);
loglog(Re_list, fit_line2, '-', 'LineWidth', line_width, 'Color', colors(2, :), ...
    'DisplayName', sprintf('Fit: $10^{%.2f} Re^{%.2f}$', coeffs2(2), coeffs2(1)));

% Fit for max_forcing_bi
log_data3 = log10(max_forcing_bi(:));
coeffs3 = polyfit(log_Re, log_data3, 1);
fit_line3 = 10.^polyval(coeffs3, log_Re);
loglog(Re_list, fit_line3, '-', 'LineWidth', line_width, 'Color', colors(4, :), ...
    'DisplayName', sprintf('Fit: $10^{%.2f} Re^{%.2f}$', coeffs3(2), coeffs3(1)));

hold off;

% Labels and legend
xlabel('Re', 'FontSize', 10, 'Interpreter', 'latex');
xlim([1e2, 2000]);
% legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend('Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
