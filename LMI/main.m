%% Using Linear Matrix Inequality to find upper bound of forcing and finite gain
%% Date: 2025/5/8
%% Outputs: u_upper_bound: upper bound of forcing,  Prop_y_u: proportion of output and input, u_tau_Lp: Lp norm of input, 
%% y_tau_Lp: Lp norm of output, beta_theorem_p: beta value for Theorem 5.1, Gamma_theorem: , gamma_cor: finite gain from Cor 5.2, L2 gain: finite gain from Theorem 5.4
%% 


% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/10.2/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))

%% used for norm_analysis for different delta
% Re_list = 200;
% delta_list = 0.5*1e-2;
% Re_list = 200;
% delta_list = 1e-4;
% Re_list = 200;
% delta_list = 1e-6;

%% used for norm_analysis for different Re
% Re_list = 1000;
% delta_list = 1e-6;
% Re_list = 1900;
% delta_list = 1e-6;

% Re_list = logspace(log10(200), log10(2000), 3);  %test
% delta_list = logspace(-6, 0, 5);
% epsilon = 0.01;
Re_list = logspace(log10(200), log10(2000), 16);
delta_list = logspace(-6, 0, 200);
B = eye(9);
C = eye(9);
D = zeros(9,9);
n = [1;0;0;0;0;0;0;0;-1];
%     gamma_results = zeros(length(Re_list), length(delta_list));
%     test_u_delta = zeros(length(Re_list), length(delta_list));
%     test_u_d_gamma = zeros(length(Re_list), length(delta_list));
%     reachable_delta = zeros(length(Re_list), length(delta_list));
%     lambda_min_results = zeros(length(Re_list), length(delta_list));
Lx = 1.75*pi;
Lz = 1.2*pi;
alpha = (2*pi)/Lx;
Beta = pi/2;
Gamma = 2*pi/Lz;
KBG = sqrt(Beta^2+Gamma^2);
KAG = sqrt(alpha^2+Gamma^2);
KABG = sqrt(alpha^2+Beta^2+Gamma^2);
deltaf = zeros(length(Re_list), length(delta_list));
u_upper_bound = zeros(length(Re_list), length(delta_list));
delta_max = zeros(1, length(Re_list));
test_u_delta = zeros(length(Re_list),length(delta_list));
u_upper_bound_max = zeros(1, length(Re_list));
prop_y_u = zeros(length(Re_list), 16);
u_tau_Lp = zeros(length(Re_list), 16);
y_tau_Lp = zeros(length(Re_list), 16);
beta_theorem_p = zeros(length(Re_list), 16);
Gamma_theorem = zeros(length(Re_list), length(delta_list));
Gamma_theorem_max = zeros(1, length(Re_list));
Gamma_theorem_min = zeros(1, length(Re_list));
gamma_cor = zeros(length(Re_list),1);
L2_gain = zeros(length(Re_list), 1);
ind_u_upper_bound_max = zeros(1,length(Re_list));
ind_delta_max = zeros(1,length(Re_list));
ind_Gamma_theorem_max = zeros(1,length(Re_list));
ind_Gamma_theorem_min = zeros(1,length(Re_list));
[RHS_J_mean_shear, nonlinear, u] = nonliner(alpha,Beta,Gamma,KBG, KAG,KABG);
delete(gcp('nocreate'));
% parpool(4);

%% Theorem 5.1
% parfor ind_Re = 1:16
%     Re = Re_list(ind_Re);
%     [local_deltaf, local_u_upper_bound, local_test_u_delta, local_prop_y_u, local_u_tau_Lp, local_y_tau_Lp, local_forcing, local_forcing_amp]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear,u,B,C);
%     deltaf{ind_Re} = local_deltaf;
%     [delta_max{ind_Re},ind_delta_max{ind_Re}] = max(local_deltaf);
%     u_upper_bound{ind_Re} = local_u_upper_bound;
%     [u_upper_bound_max{ind_Re},ind_u_upper_bound_max{ind_Re}] = max(local_u_upper_bound);
%     test_u_delta{ind_Re} = local_test_u_delta;
%     prop_y_u{ind_Re} = local_prop_y_u;
%     u_tau_Lp{ind_Re} = local_u_tau_Lp;
%     y_tau_Lp{ind_Re} = local_y_tau_Lp;
%     forcing_LMI{ind_Re} = local_forcing;
%     forcing_bi{ind_Re} = local_forcing_amp;
% end

for ind_Re = 1:length(Re_list)
    Re = Re_list(ind_Re);
    [local_deltaf, local_u_upper_bound, local_test_u_delta, local_prop_y_u, local_u_tau_Lp, local_y_tau_Lp, local_forcing, local_beta_theorem_p, local_Gamma_theorem, local_gamma_cor, local_L2_gain]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear,u,B,C, alpha,Beta,Gamma,KBG, KAG,KABG);
%     LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear,u,B,C);
    deltaf(ind_Re, :) = local_deltaf;
    [delta_max(ind_Re),ind_delta_max(ind_Re)] = max(local_deltaf);
    u_upper_bound(ind_Re, :) = local_u_upper_bound;
    [u_upper_bound_max(ind_Re),ind_u_upper_bound_max(ind_Re)] = max(local_u_upper_bound);
    test_u_delta(ind_Re,:) = local_test_u_delta;
    prop_y_u(ind_Re,:) = local_prop_y_u;
    u_tau_Lp(ind_Re,:) = local_u_tau_Lp;
    y_tau_Lp(ind_Re,:) = local_y_tau_Lp;
%     forcing_LMI{ind_Re} = local_forcing;
    % forcing_bi(ind_Re,:) = local_forcing_amp;
    beta_theorem_p(ind_Re,:) = local_beta_theorem_p;
    Gamma_theorem(ind_Re,:)= local_Gamma_theorem;
    non_nan_mask = ~isnan(local_Gamma_theorem);
    non_nan_indices = find(non_nan_mask);
    
    if ~isempty(non_nan_indices)
        % Get the last non-NaN index and value
        ind_last_non_nan = non_nan_indices(end);
        last_non_nan_Gamma = local_Gamma_theorem(ind_last_non_nan);
    else
        % If all values are NaN, assign NaN
        ind_last_non_nan = NaN;
        last_non_nan_Gamma = NaN;
    end
    
    % Store results
    Gamma_theorem_last(ind_Re) = last_non_nan_Gamma;
    ind_Gamma_theorem_last(ind_Re) = ind_last_non_nan;
    [Gamma_theorem_max(ind_Re), ind_Gamma_theorem_max(ind_Re)]= max(local_Gamma_theorem);
    [Gamma_theorem_min(ind_Re), ind_Gamma_theorem_min(ind_Re)]= min(local_Gamma_theorem);
    gamma_cor(ind_Re,:) = local_gamma_cor; 
    L2_gain(ind_Re,:) = local_L2_gain;
end


%     if any(delta_max == 0)
%         ind_Re = find(delta_max);
%     end

% log_Re_list = log10(Re_list);
% log_delta_max = log10(delta_max);
% 
% % line fitting to the data
% coeffs = polyfit(log_Re_list, log_delta_max, 1);
% sigma = coeffs(1);
% A_log = coeffs(2);
save ODE45_9D_oldLMI_function_NEW_51_2_5_test.mat
% figure;
% loglog(Re_list, delta_max, 'o'); hold on;
% fit_line = 10.^polyval(coeffs, log_Re_list);
% loglog(Re_list, fit_line, '-');
% xlabel('Re');
% ylabel('\delta_{p}');
% title(['Scaling exponent \sigma = ', num2str(sigma)]);
% 
% disp(['Scaling exponent (sigma): ', num2str(sigma)]);
% disp(['Intercept (A): ', num2str(A_log)]);
% saveas(gcf,'figure.fig');

