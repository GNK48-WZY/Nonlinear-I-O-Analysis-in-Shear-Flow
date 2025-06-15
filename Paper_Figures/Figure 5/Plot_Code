beta_theorem_inf = readmatrix("beta_theorem_inf.csv");
local_Gamma_theorem = readmatrix("local_Gamma_theorem.csv");
local_beta_theorem_p = readmatrix("local_beta_theorem_p.csv");
local_u_tau_Lp = readmatrix("local_u_tau_Lp.csv");
local_y_tau_Lp = readmatrix("local_y_tau_Lp.csv");
p_list = readmatrix("p_list.csv");
u_tau_Linf = readmatrix("u_tau_Linf.csv");
y_tau_Linf = readmatrix("y_tau_Linf.csv");


colors = lines(16);
line_styles = ["-","--","-.",":",""];
markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
figure;
h1 = loglog(p_list, local_y_tau_Lp, '--r', 'LineWidth',1.5,'Marker',markers{1},'Color', colors(1, :)); hold on;
h2 = loglog(p_list, local_Gamma_theorem(ind_delta) .*local_u_tau_Lp + local_beta_theorem_p, '-.', 'LineWidth',1.5,'Color', colors(2, :), 'Marker',markers{2});
h3 = loglog(p_list, y_tau_Linf.*ones(size(p_list)),'-r','Marker',markers{3}, 'LineWidth',1.5,'Color', colors(3, :));
h4 = loglog(p_list, (local_Gamma_theorem(ind_delta) .* u_tau_Linf + beta_theorem_inf).*ones(size(p_list)),'-k', 'LineWidth',1.5,'Color', colors(4, :),'Marker',markers{4});
xlabel('p', 'Interpreter', 'latex');
ylim([1e-8, 1e3]);
set(gca, 'FontSize', 12);
legend([h1, h2, h3, h4], {...
'$\|\mathbf{y}_{\tau}\|_{\mathcal{L}_p}$', ...
    '$\gamma \|\mathbf{f}_{\tau}\|_{\mathcal{L}_p} + \beta$', ...
    '$\|\mathbf{y}_{\tau}\|_{\mathcal{L}_\infty}$', ...
    '$\gamma \|\mathbf{f}_{\tau}\|_{\mathcal{L}_\infty} + \beta_\infty$'}, ...
    'Interpreter', 'latex');
grid on;
hold off;

