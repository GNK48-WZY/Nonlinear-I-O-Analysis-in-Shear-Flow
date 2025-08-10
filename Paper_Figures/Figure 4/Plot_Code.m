clear all

unzip("no_forcing_Re_1000_1e_6.zip")
unzip("no_forcing_Re_1900_1e_6.zip")
unzip("no_forcing_Re_200_1e_6.zip")
unzip("no_initial_Re_200_1e_6.zip")
unzip("no_initial_Re_1000_1e_6.zip")
unzip("no_initial_Re_1900_1e_6.zip")
unzip("t.zip")

Upper_bound_1e_6 = readmatrix("Upper_bound_1e_6.csv");
no_forcing_Re_1000_1e_6 = readmatrix("no_forcing_Re_1000_1e_6.csv");
no_forcing_Re_1900_1e_6 = readmatrix("no_forcing_Re_1900_1e_6.csv");
no_forcing_Re_200_1e_6 = readmatrix("no_forcing_Re_200_1e_6.csv");
no_initial_Re_200_1e_6 = readmatrix("no_initial_Re_200_1e_6.csv");
no_initial_Re_1000_1e_6 = readmatrix("no_initial_Re_1000_1e_6.csv");
no_initial_Re_1900_1e_6 = readmatrix("no_initial_Re_1900_1e_6.csv");
t = readmatrix("t.csv");



% Figure 4a
line_styles = {"-","--","-.",":",""};
markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
T = 70000;
t = linspace(0, T, 1000000);
figure;
loglog(t, no_forcing_Re_200_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{1}); hold on;
loglog(t, no_forcing_Re_1000_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{2});
loglog(t, no_forcing_Re_1900_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{3});

% Plot the horizontal line for the upper bound
loglog(t,Upper_bound_1e_6, 'LineWidth', 1.5, 'DisplayName', 'upper bound','LineStyle', line_styles{4});
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|\mathbf{y}_{\tau}(t)\|$', 'Interpreter', 'latex');
set(gca, 'FontSize', 13);
legend('Re$ = 200$','Re $= 1000$','Re $= 1900$','$\xi$','Interpreter', 'latex');
grid on;
hold off;

disp('Figure 4a');


% Figure 4b
figure;
loglog(t, no_initial_Re_200_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{1}); hold on;
loglog(t, no_initial_Re_1000_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{2});
loglog(t, no_initial_Re_1900_1e_6, 'LineWidth', 1.5,'LineStyle', line_styles{3});
 
% Plot the horizontal line for the upper bound
loglog(t,Upper_bound_1e_6, 'LineWidth', 1.5, 'DisplayName', 'upper bound','LineStyle', line_styles{4});
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|\mathbf{y}_{\tau}(t)\|$', 'Interpreter', 'latex');
set(gca, 'FontSize', 13);
legend('Re$ = 200$','Re $= 1000$','Re $= 1900$','$\xi$','Interpreter', 'latex');
grid on;
hold off;


disp('Figure 4b');



