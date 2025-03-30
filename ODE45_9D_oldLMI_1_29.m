% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/10.2/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))


% Re_list = 200;
% delta_list = 0.5*1e-2;


% Re_list = 200;
% delta_list = 1e-4;

% Re_list = 200;
% delta_list = 1e-6;

% Re_list = 1000;
% delta_list = 1e-6;
% 
% Re_list = 2000;
% delta_list = 1e-6;
% Re_list = logspace(log10(200), log10(2000), 3);
% delta_list = logspace(-6, 0, 5);
% Re_list= Re_list(1:2);
% delta_list=delta_list(1);
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


function [local_deltaf, local_u_upper_bound, local_test_u_delta, local_prop_y_u, local_u_tau_Lp, local_y_tau_Lp, local_forcing, local_beta_theorem_p, local_Gamma_theorem, local_gamma_cor, local_L2_gain]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear, u,B,C, alpha, Beta, Gamma,KBG, KAG,KABG)
    local_u_upper_bound = NaN*ones(1, length(delta_list));
    local_deltaf = NaN*ones(1, length(delta_list));
    local_test_u_delta = NaN*ones(1,length(delta_list));
    local_Gamma_theorem = NaN*ones(1,length(delta_list));
    linear_term = [(Beta^2)/Re;
        ((4*Beta^2)/3 + Gamma^2)/Re;
        (Beta^2+Gamma^2)/Re;
        (3*alpha^2+4*Beta^2)/(3*Re);
        (alpha^2+Beta^2)/Re;
        (3*alpha^2+4*Beta^2+3*Gamma^2)/(3*Re);
        (alpha^2+Beta^2+Gamma^2)/Re;
        (alpha^2+Beta^2+Gamma^2)/Re;
        (9*Beta^2)/Re];
    RHS_R_viscous = diag(double(linear_term));
    A = RHS_J_mean_shear - RHS_R_viscous;

   %% Theorem 5.4
    D = zeros(9,9);
    I = eye(size(A));
%     s = tf('s'); % Laplace variable
%     G = C * inv(s * I - A) * B + D;   
%     % Compute the L2 gain (H-infinity norm)
%     local_L2_gain = norm(G, Inf);
    sys = ss(A, B, C, D); 
    local_L2_gain = hinfnorm(sys);

    % Display the result
    disp(['The L2 gain of the system is: ', num2str(local_L2_gain)]);


for ind_delta = 1:length(delta_list)
    yalmip('clear');
    delta = delta_list(ind_delta);
    delta2 = delta^2;

    %             if ~isempty(n)  %we have a non-zero null space of nonlinear term.
    %                 n = double(n);
    %                 lambda=sdpvar(length(A),1);
    %                 kappa=sdpvar(length(A),1);
    %                 %n=double(n);
    %                 %implementing the inequality constraint like the Lagrange
    %                 %multiplier
    %                 for ind_e=1:length(A)
    %                       e=zeros(length(A),1);
    %                       e(ind_e)=1;
    %                       lambda_af=lambda_af+lambda(ind_e,1)*e*n';
    %                       kappa_ff=kappa_ff+kappa(ind_e,1)*(e*n'+n*e');
    %                 end
    %             else
    %                 lambda=sdpvar(1,1);
    %                 lambda_af=zeros(size(A));
    %                 kappa=sdpvar(1,1);
    %                 kappa_ff=zeros(size(A));
    %             end
    [~, ~, lambda_af, kappa_ff] = check(A, C, B);
    [F_square] = F__square(nonlinear, u);

    %             s = sdpvar(length(F_square), 1);
    %             s_bound = zeros(size(A));
    %             for m_ind = 1:length(F_square)
    %                  s_bound = s_bound + s(m_ind) * delta2 * double(F_square{m_ind});
    %             end
    %
    %             if length(s) <= length(A)
    %                 diag_s = diag(s);
    %             else
    %                 diag_s = diag(s(1:length(A))) + eye(size(A)) * s(length(A) + 1);
    %             end
    [s, diag_s, s_bound] = s_diag_s(F_square, A, delta2);
    %% Corollary 5.2
    % Initialize variables and parameters
    P = lyap(transpose(A), I);
    D = zeros(9,9);
    
    % Calculate eigenvalues of P
    lambda_max_P = max(eig(P)); % Maximum eigenvalue of P
    lambda_min_P = min(eig(P)); % Minimum eigenvalue of P
    
    % Calculate L2 norms
    norm_B = norm(B, 2); % ||B||_2
    norm_C = norm(C, 2); % ||C||_2
    norm_D = norm(D, 2); % ||D||_2
    
    % Compute gamma
    local_gamma_cor = norm_D + (2 * lambda_max_P^2 * norm_B * norm_C) / lambda_min_P;
    
    % Initialize a range of p values to test
    p_values = logspace(0, log10(10000), 16);
    beta_cor_p_values = zeros(size(p_values)); % To store beta for each p
    cor2_deltaf = sqrt(value(delta2)) * sqrt(lambda_min_P / lambda_max_P);
    % p=inf, rho =1
    beta_coro_inf = 1 * norm_C * cor2_deltaf * sqrt(lambda_max_P / lambda_min_P);

    % Loop over p values
    for i = 1:length(p_values)
        p = p_values(i);
        
 
        rho = ((2 * lambda_max_P) / p)^(1 / p);
        rho = double(rho);
        
        % Compute delta and beta (assuming co2_deltaf is defined)
        cor2_deltaf = sqrt(value(delta2)) * sqrt(lambda_min_P / lambda_max_P);
        beta_cor_p = rho * norm_C * cor2_deltaf * sqrt(lambda_max_P / lambda_min_P);
        
        % Store beta
        beta_cor_p_values(i) = beta_cor_p;
    end


    
    % Construct LMI conditions
    epsilon=0.01;
    sdp_options = sdpsettings('solver', 'mosek');
    options.bisection.absgaptol = 1e-16;
    [result_yalmip, P, lambda_af, kappa_ff] = LMI(A,s_bound, epsilon, lambda_af, kappa_ff, diag_s, B, C, s, sdp_options);


    if result_yalmip.problem == 0

        lambda_af_value = value(lambda_af); % Lagrange mutiplier
        kappa_ff_value = value(kappa_ff); % Lagrange mutiplier
        optp = value(P);
        umax = max(eig(optp)); % c2
        umin = min(eig(optp)); % c1
        c3 = epsilon;
        c2 = umax;
        c1 = umin;
        L = norm(B);
        c4 = norm(2*value(P));
        local_deltaf(ind_delta) = sqrt(value(delta2)) * sqrt(umin / umax);
        eta_1 = norm(C);
        eta_2 = 0;        
        local_Gamma_theorem(ind_delta) = eta_2 + (eta_1*c2*c4*L)/(c1*c3);
        % For Lp=inf
        beta_theorem_inf = eta_1*local_deltaf(ind_delta)*sqrt(c2/c1)*1;
        
        local_u_upper_bound(ind_delta) = (c1*c3*delta)/(c2*c4*L);
        disp(['Feasible solution found for Re = ', num2str(Re), ', delta = ', num2str(delta_list(ind_delta))]);
        disp(['P: ', mat2str(optp)]);
        disp(['umax: ', num2str(umax)]);
        disp(['umin: ', num2str(umin)]);

        T = 20000;
        omega=1;
        [norm_u_simulations(ind_delta), x, ~, ~,~, t, local_forcing, y, dt, ~] = compute_norm_u(Re, A, B, T, omega,local_deltaf(ind_delta), local_u_upper_bound(ind_delta),alpha,Beta,Gamma,KBG, KAG,KABG);
        
        %% Theorem 5.1          
        y_norm = sqrt(sum(y.^2,2));
        forcing_norm=sqrt(sum(local_forcing.^2,2));
        u_tau_Linf = max(abs(forcing_norm));  % Linf norm of u 
        y_tau_Linf = max(abs(y_norm));  % Linf norm of y 
        inequality_holds_Linf = y_tau_Linf <= local_Gamma_theorem(ind_delta) * u_tau_Linf + beta_theorem_inf;
        if inequality_holds_Linf
                fprintf('The inequality holds: ||y_tau||_Linf <= gamma * ||u_tau||_Linf + beta.\n');
        else
                fprintf('The inequality does NOT hold.\n');
        end
        
%         p_list = [1,2,3,4];
        p_list = logspace(log10(1), log10(10000), 16);% specify the norm type, e.g., p=2 for L2 norm
        %p_list = [1,2,10,20,40,80,100,200,400,800,1000,2000,4000,5000,8000,10000];
       % add for loop to get Lp stable for different p, and then we need get the norm_u_simulation over t and norm_u_theorem over t to show the upper bound
        for ind_p = 1:length(p_list)
            p = p_list(ind_p);
            local_beta_theorem_p(ind_p) = eta_1*local_deltaf(ind_delta)*(sqrt(c2/c1))*(((2*c2)/(c3*p))^(1/p));
            local_u_tau_Lp(ind_p) = double((sum(abs(sym(forcing_norm)).^p) * dt)^(1/p));  % Lp norm of u over the interval
            local_y_tau_Lp(ind_p) = double((sum(abs(sym(y_norm)).^p) * dt)^(1/p));  % Lp norm of y over the interval
            
            local_prop_y_u(ind_p) = local_y_tau_Lp(ind_p)/local_u_tau_Lp(ind_p);
            inequality_holds(ind_p) = local_y_tau_Lp(ind_p) <= local_Gamma_theorem(ind_delta) * local_u_tau_Lp(ind_p) + local_beta_theorem_p(ind_p);
            fprintf('Lp norm of input u(t): %.18f\n', local_u_tau_Lp(ind_p));
            fprintf('Lp norm of output y(t): %.18f\n', local_y_tau_Lp(ind_p));
    %         fprintf('Gamma * Lp norm of u + Beta: %.18f\n', Gamma_theorem* u_tau_Linf + beta_theorem_inf);
            fprintf('Gamma * Lp norm of u + Beta: %.18f\n', local_Gamma_theorem(ind_delta)* local_u_tau_Lp(ind_p) + local_beta_theorem_p(ind_p));
    
            if inequality_holds
                fprintf('The inequality holds: ||y_tau||_Lp <= gamma * ||u_tau||_Lp + beta.\n');
            else
                fprintf('The inequality does NOT hold.\n');
            end

            
%             upper_bound_y = local_Gamma_theorem(ind_delta)*local_u_upper_bound(ind_delta) +beta_theorem_inf*ones(size(t));
%             figure;
% %             Plot the norm from simulation
% %             loglog(t, u_square, 'LineWidth', 1.5, 'DisplayName', 'norm usimulation'); hold on;
%             loglog(t, y_norm, 'LineWidth', 1.5, 'DisplayName', 'norm y simulation'); hold on;
% 
%             % Plot the horizontal line for the upper bound
%             loglog(t,upper_bound_y, 'DisplayName', 'upper bound');
%             %loglog(t, local_u_upper_bound(ind_delta) * ones(size(t)), '--r', 'LineWidth', 1.5, 'DisplayName', 'Upper Bound');
%             % Plot Lp norms as horizontal lines
% %             loglog(t, u_tau_Lp(ind_p) * ones(size(t)), '--', 'DisplayName', ['Lp Norm (p=', num2str(p_list(ind_p)), ')']); % Lp_norm equal to zero after p=20, may be is numerical problem
% %             loglog(t, u_tau_Linf * ones(size(t)), '-.', 'DisplayName', 'L∞ Norm'); % Plot Linf norm
%             xlabel('Time (t)');
%             ylabel('Norm Values');
%             title(['Norm Analysis for \delta = ', num2str(delta_list(ind_delta))]);
%             legend('show');
%             grid on;
%             hold off;

 
        end
%% figure p vs. beta(Them 5.1 and Cor 5.2)

        % figure;
        % loglog(p_list, local_beta_theorem_p, '-o', 'LineWidth', 1.5); % First line
        % hold on;
        % loglog(p_values, beta_cor_p_values, '-o', 'LineWidth', 1.5); % Second line
        % hold off;
        % 
        % % Add title and labels
        % title('Effect of p on Beta', 'FontSize', 14);
        % xlabel('p', 'FontSize', 12);
        % ylabel('\beta (Beta)', 'FontSize', 12);
        % 
        % % Customize grid, ticks, and legend
        % grid on;
        % xticks(p_values);
        % legend('\beta_{theorem} vs. p', '\beta_{cor} vs. p', 'FontSize', 10, 'Location', 'best');
        % 
        % % Adjust font size for the axes
        % set(gca, 'FontSize', 12);
%%
%         % 
        colors = lines(16);
        line_styles = {"-","--","-.",":",""};
        markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
        figure;
        loglog(p_list, local_y_tau_Lp, '--r', 'LineWidth',1.5,'Marker',markers{1},'Color', colors(1, :),'DisplayName','LHS of Inequality (p=1-10^4 )'); hold on;
        loglog(p_list, local_Gamma_theorem(ind_delta) .*local_u_tau_Lp + local_beta_theorem_p, '-.', 'LineWidth',1.5,'Color', colors(2, :), 'Marker',markers{2},'DisplayName', 'RHS of Inequality (p=1-10^4)');
        loglog(p_list, y_tau_Linf.*ones(size(p_list)),'-r','Marker',markers{3}, 'LineWidth',1.5,'Color', colors(3, :),'DisplayName', 'LHS of Inequality (p=\infty)');
        loglog(p_list, (local_Gamma_theorem(ind_delta) .* u_tau_Linf + beta_theorem_inf).*ones(size(p_list)),'-k', 'LineWidth',1.5,'Color', colors(4, :),'Marker',markers{4},'DisplayName' ,'RHS of Inequality (p=\infty)');
        xlabel('p', 'Interpreter', 'latex');
%         title('Gap between LHS and RHS change over p');
        set(gca, 'FontSize', 13);
        legend('show');
        grid on;
        hold off;


        if norm_u_simulations(ind_delta) <= delta
            local_test_u_delta(ind_delta) = 1;
        else
            local_test_u_delta(ind_delta) = 0;
        end

    else
        disp(['No feasible solution found for Re = ', num2str(Re), ', delta = ', num2str(delta_list(ind_delta))]);
        disp(['YALMIP error: ', yalmiperror(result_yalmip.problem)]);
    end
end

end

%% 9D Shear flow model
function [RHS_J_mean_shear,nonlinear, u] = nonliner(alpha,Beta,Gamma,KBG, KAG,KABG)

u = sym('u', [9, 1]);
term1 = -sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG+sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG;
term2 = (5/3)*sqrt(2/3)*Gamma^2*u(4)*u(6)/KAG-Gamma^2*u(5)*u(7)/(sqrt(6)*KAG) ...
    -alpha*Beta*Gamma*u(5)*u(8)/(sqrt(6)*KAG*KABG)-sqrt(3/2)*Beta*Gamma*u(1)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(3)*u(9)/KBG;
term3 = 2*alpha*Beta*Gamma*(u(4)*u(7)+u(5)*u(6))/(sqrt(6)*KAG*KBG)+(Beta^2*(3*alpha^2+Gamma^2)-3*Gamma^2*(alpha^2+Gamma^2))*u(4)*u(8)/(sqrt(6)*KAG*KBG*KABG);
term4 = -alpha*u(1)*u(5)/sqrt(6)-10*alpha^2*u(2)*u(6)/(3*sqrt(6)*KAG)  ...
    -sqrt(3/2)*alpha*Beta*Gamma*u(3)*u(7)/KAG*KBG-sqrt(3/2)*alpha^2*Beta^2*u(3)*u(8)/KAG*KBG*KABG-alpha*u(5)*u(9)/sqrt(6);
term5 =  alpha*u(1)*u(4)/sqrt(6)+alpha^2*u(2)*u(7)/(sqrt(6)*KAG)-alpha*Beta*Gamma*u(2)*u(8)/(sqrt(6)*KAG*KABG)+alpha*u(4)*u(9)/sqrt(6)+2*alpha*Beta*Gamma*u(3)*u(6)/(sqrt(6)*KAG*KBG);
term6 =  alpha*u(1)*u(7)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(1)*u(8)/KABG  ...
    +10*(alpha^2-Gamma^2)*u(2)*u(4)/(KAG*3*sqrt(6))-2*sqrt(2/3)*u(3)*u(5)*alpha*Beta*Gamma/(KAG*KBG)+alpha*u(7)*u(9)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(8)*u(9)/KABG;
term7 = -alpha*(u(1)*u(6)+u(6)*u(9))/sqrt(6)+(Gamma^2-alpha^2)*u(2)*u(5)/(sqrt(6)*KAG)+alpha*Beta*Gamma*u(3)*u(4)/(sqrt(6)*KAG*KBG);
term8 = 2*alpha*Beta*Gamma*u(2)*u(5)/(sqrt(6)*KAG*KABG)+Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*u(3)*u(4)/(sqrt(6)*KAG*KBG*KABG);
term9 = sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG;
nonlinear= [term1;
    term2;
    term3;
    term4;
    term5;
    term6;
    term7;
    term8;
    term9];



nonlinear_gradient = sym(zeros(length(u), length(u)));
for i = 1:length(nonlinear)
    nonlinear_gradient(i, :) = gradient(nonlinear(i), u);
end
%nonlinear_gradinet = gradient(nonlinear, u);
a_bar = [1;0;0;0;0;0;0;0;0];
nonlinear_gradient_sub = double(subs(nonlinear_gradient, u, a_bar));
RHS_J_mean_shear = nonlinear_gradient_sub;

end


function [F_square] = F__square(nonlinear, u)

F_square = cell(1, length(nonlinear));
for n_ind = 1:length(nonlinear)
    %initialize F{n_ind}
    F{n_ind} = sym(zeros(length(u), length(u)));
    for x_ind = 1:length(u)
        for y_ind = 1:length(u)
            F{n_ind}(x_ind, y_ind) = 1/2 * diff(diff(nonlinear(n_ind), u(x_ind)), u(y_ind));
        end
    end

    %convert F{n_ind} to double if it's symbolic
    [V, D] = eig(double(F{n_ind}));
    %store F_square as double
    F_square{n_ind} = double(V * D^2 /(V));
end

end


%% LMI
function[result_yalmip, P,lambda_af,kappa_ff] = LMI(A,s_bound, epsilon, lambda_af, kappa_ff, diag_s, B, C, s, sdp_options)

P = sdpvar(size(A, 1), size(A, 2));
G11 = A' * P + P * A + s_bound + epsilon*eye(size(A));
G12 = P*B + C'*lambda_af;
G22 = -diag_s + kappa_ff;

% LMI
G = [G11, G12;
    G12', G22;];

V_ineq = P - epsilon * eye(size(A));

% Solve the LMI
F = [V_ineq >= 0, G <= 0, s >= 0];
result_yalmip = optimize(F, [], sdp_options);
end

function[s, diag_s, s_bound] = s_diag_s(F_square, A, delta2)
s = sdpvar(length(F_square), 1);
s_bound = zeros(size(A));
for m_ind = 1:length(F_square)
    s_bound = s_bound + s(m_ind) * delta2 * double(F_square{m_ind});
end

if length(s) <= length(A)
    diag_s = diag(s);
else
    diag_s = diag(s(1:length(A))) + eye(size(A)) * s(length(A) + 1);
end
end


function[lambda, kappa, lambda_af, kappa_ff] = check(A, C, B)
n = [1;0;0;0;0;0;0;0;-1];
lambda_af=zeros(size(C,1),size(C,1));
kappa_ff=zeros(size(B,2),size(B,2));
if ~isempty(n)  %%we have a non-zero null space of nonlinear term.
    n = double(n);
    lambda=sdpvar(length(A),1);
    kappa=sdpvar(length(A),1);
    %n=double(n);
    %implementing the inequality constraint like the Lagrange
    %multiplier
    for ind_e=1:length(A)
        e=zeros(length(A),1);
        e(ind_e)=1;
        lambda_af=lambda_af+lambda(ind_e,1)*e*n';
        kappa_ff=kappa_ff+kappa(ind_e,1)*(e*n'+n*e');
    end
else
    lambda=sdpvar(1,1);
    lambda_af=zeros(size(A));
    kappa=sdpvar(1,1);
    kappa_ff=zeros(size(A));
end
end

%% forcing
%function d_forcing_func = generate_d_forcing(delta, T, omega)
%phi = rand(1,1)*2*pi;
% force=randn(9,1);
% d_forcing_func = @(t) (delta)*force/norm(force);
%[sin(omega*t)*cos(phi);sin(omega*t)*sin(phi);sin(omega*t)*cos(phi); sin(omega*t)*sin(phi);sin(omega*t)*cos(phi); sin(omega*t)*sin(phi); sin(omega*t)*cos(phi); sin(omega*t)*sin(phi);sin(omega*t)*cos(phi)];
%end
% 
% function forcing =d_forcing_random(t,u_upper_bound)
% force=randn(9,1);
% forcing=u_upper_bound*force/norm(force);
% end


%% Norm_u computation
function [norm_u, u, G, sys, max_SVD, t, local_forcing, y, dt, a_square] = compute_norm_u(Re, A, B, T, omega, deltaf, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG)


ini=randn(9,1);
u0 = deltaf*ini/norm(ini);
% u0 =zeros(9,1);

T = 70000;
tspan = linspace(0, T, 10000);

I = eye(size(A));
C = eye(9);
D = zeros(9, 9);

G = C * inv(1i * omega * I - A) * B;
max_SVD = max(svd(G));
sys = ss(A, B, C, D);
% u0 = zeros(9,1);
% tol = 1e-9;  % Tolerance for convergence
% max_iter = 1000;  % Maximum fixed-point iterations per step

% dt = diff(tspan);
% dt = dt(1);
%[t, u, forcing] = ode45(@(t,u)  ode_system_123(t, u, A, B, u_upper_bound), tspan, u0);
[t, u, local_forcing] = rk5_solver(@(t,u)  ode_system_123(t, u, A, B, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG), tspan, u0);

% Linear Equation analytical solution
% for t_ind = 1:length(tspan)
%     syms t_sym real; % Create a symbolic variable for time
%     u_sym = expm(A *tspan(t_ind)) * u0; % Analytical solution: u(t)=exp(At)*u0
%     u_analytical(:,t_ind) = double(u_sym); % Store symbolic expression
% end 

%old norm integrate over time
dt=diff(tspan);
dt=dt(1);
% u_square=sum(u.^2,2);
% norm_u=sqrt(sum(u_square(2:end).*dt)/max(tspan));

%1 norm maximal over time
a_square=sqrt(sum(u.^2,2));
%norm_u = sqrt(sum(u_square.^2)*dt);
y =C*u';
y = y'; % output
norm_u=max(a_square);
% 
% 
% % %Bisection
% max_u = 0.01; 
% min_u = 0;  
% % Tolerance and maximum iterations
% tol = 1e-14;
% maxIter = 100;
% 
% % Call the bisection method
% % u0=zeros(9,1);
% % [forcing_amp, iterations] = bisectionMethod(min_u, max_u, tol, maxIter, A, B,tspan, u0);
% num_simulations = 20;  % Number of simulations per step
% [forcing_amp, iterations] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations);
% % Display results
% fprintf('Approximated root: %.20f\n', forcing_amp);
% fprintf('Number of iterations: %d\n', iterations);
% % 

end



function [du_dt, forcing] = ode_system_123(t, u, A, B, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG)


force=randn(9,1);
forcing=u_upper_bound*force/norm(force);

term1 = -sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG+sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG;
term2 = (5/3)*sqrt(2/3)*Gamma^2*u(4)*u(6)/KAG-Gamma^2*u(5)*u(7)/(sqrt(6)*KAG) ...
    -alpha*Beta*Gamma*u(5)*u(8)/(sqrt(6)*KAG*KABG)-sqrt(3/2)*Beta*Gamma*u(1)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(3)*u(9)/KBG;
term3 = 2*alpha*Beta*Gamma*(u(4)*u(7)+u(5)*u(6))/(sqrt(6)*KAG*KBG)+(Beta^2*(3*alpha^2+Gamma^2)-3*Gamma^2*(alpha^2+Gamma^2))*u(4)*u(8)/(sqrt(6)*KAG*KBG*KABG);
term4 = -alpha*u(1)*u(5)/sqrt(6)-10*alpha^2*u(2)*u(6)/(3*sqrt(6)*KAG)  ...
    -sqrt(3/2)*alpha*Beta*Gamma*u(3)*u(7)/KAG*KBG-sqrt(3/2)*alpha^2*Beta^2*u(3)*u(8)/KAG*KBG*KABG-alpha*u(5)*u(9)/sqrt(6);
term5 =  alpha*u(1)*u(4)/sqrt(6)+alpha^2*u(2)*u(7)/(sqrt(6)*KAG)-alpha*Beta*Gamma*u(2)*u(8)/(sqrt(6)*KAG*KABG)+alpha*u(4)*u(9)/sqrt(6)+2*alpha*Beta*Gamma*u(3)*u(6)/(sqrt(6)*KAG*KBG);
term6 =  alpha*u(1)*u(7)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(1)*u(8)/KABG  ...
    +10*(alpha^2-Gamma^2)*u(2)*u(4)/(KAG*3*sqrt(6))-2*sqrt(2/3)*u(3)*u(5)*alpha*Beta*Gamma/(KAG*KBG)+alpha*u(7)*u(9)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(8)*u(9)/KABG;
term7 = -alpha*(u(1)*u(6)+u(6)*u(9))/sqrt(6)+(Gamma^2-alpha^2)*u(2)*u(5)/(sqrt(6)*KAG)+alpha*Beta*Gamma*u(3)*u(4)/(sqrt(6)*KAG*KBG);
term8 = 2*alpha*Beta*Gamma*u(2)*u(5)/(sqrt(6)*KAG*KABG)+Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*u(3)*u(4)/(sqrt(6)*KAG*KBG*KABG);
term9 = sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG;
nonlinear= [term1;
    term2;
    term3;
    term4;
    term5;
    term6;
    term7;
    term8;
    term9];

du_dt = A * u + B * forcing+ nonlinear;
end

function [t, x, forcing] = rk5_solver(ode_func, tspan, u0)
    % RK5 Solver for ODEs
    % Inputs:
    %   ode_func - Function handle for the ODE system (f(t, u))
    %   tspan - Time vector [t0, t1, ..., tn]
    %   u0 - Initial condition (row vector)
    % Outputs:
    %   t - Time vector (same as tspan)
    %   x - Solution matrix (size: length(tspan) x length(u0))
    %   forcing - Forcing term evaluated at each step (optional)
    
    % Number of time steps
    N = length(tspan);
    
    % Preallocate solution and forcing arrays
    x = zeros(N, length(u0));
    forcing = zeros(N, length(u0));
    x(1, :) = u0; % Set initial condition
    
    % Coefficients for RK5 (Dormand-Prince example)
    c = [0, 1/4, 3/8, 12/13, 1, 1/2];
    b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    a = [
        0,      0,      0,       0,       0,      0;
        1/4,    0,      0,       0,       0,      0;
        3/32,   9/32,   0,       0,       0,      0;
        1932/2197, -7200/2197, 7296/2197, 0,       0,  0;
        439/216, -8,    3680/513, -845/4104, 0,    0;
        -8/27,   2,    -3544/2565, 1859/4104, -11/40, 0
    ];
    
    % Time-stepping loop
    for i = 1:N-1
        h = tspan(i+1) - tspan(i); % Time step size
        u_current = x(i, :)'; % Current state
        
        % Compute RK5 coefficients (k1 to k6)
        k1 = h * ode_func(tspan(i), u_current);
        k2 = h * ode_func(tspan(i) + c(2)*h, u_current + a(2,1)*k1);
        k3 = h * ode_func(tspan(i) + c(3)*h, u_current + a(3,1)*k1 + a(3,2)*k2);
        k4 = h * ode_func(tspan(i) + c(4)*h, u_current + a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3);
        k5 = h * ode_func(tspan(i) + c(5)*h, u_current + a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4);
        k6 = h * ode_func(tspan(i) + c(6)*h, u_current + a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5);
        
        % Update state
        x(i+1, :) = (u_current + b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4 + b(5)*k5 + b(6)*k6)';
        
        % Optionally, store forcing term for analysis
        [~, forcing(i+1, :)] = ode_func(tspan(i), u_current);
    end
    
    t = tspan; % Output time vector
end


% % Runge-Kutta method
% function [t, x, forcing] = rk4_solver(ode_func, tspan, u0)
%     % Number of time steps
%     N = length(tspan);
%     % Preallocate u matrix (number of time steps by 9 states)
%     x = zeros(N, length(u0));
%     x(1, :) = u0;  % Set initial condition
%     forcing = zeros(N, length(u0));
%     % Time stepping
%     for i = 1:N-1
%         h = tspan(i+1) - tspan(i);  % Calculate step size
%         % Perform one RK4 step and store result in the next row
%         [x(i+1, :),forcing(i+1,:)] = rk4_step(x(i, :), h, ode_func);
% 
%     end
% 
%     t = tspan;  % Output the time vector
% end
% 
% function [u_next, forcing] = rk4_step(u_current, h, ode_func)
%     % Calculate the RK4 coefficients
%     k1 = h * ode_func(0, u_current')';
%     k2 = h * ode_func(0, (u_current' + 0.5 * k1'))';
%     k3 = h * ode_func(0, (u_current' + 0.5 * k2'))';
%     k4 = h * ode_func(0, (u_current' + k3'))';
%     [~,forcing] = ode_func(0, u_current');
%     forcing = forcing';
%     % Update the state using RK4 formula
%     u_next = u_current + (k1 + 2*k2 + 2*k3 + k4) / 6;
% end

% 
% % Bisection Method
% function [forcing_amp, iterations] = bisectionMethod(min_u, max_u, tol, maxIter, A, B,tspan, u0)
%     % Bisection Method to find the root of a function
%     % f: function handle
%     % a, b: interval [a, b] where f(a) * f(b) < 0
%     % tol: tolerance for convergence
%     % maxIter: maximum number of iterations
%     % Returns:
%     %   root: the approximated root
%     %   iterations: number of iterations performed
% 
%     % Initialize variables
%     iterations = 0;
%     forcing_amp = (max_u + min_u) / 2; % Midpoint
% 
%     % Start iteration loop
%     while (max_u - min_u) / 2 > tol && iterations < maxIter
%         iterations = iterations + 1;
% 
%         % Update the midpoint
%         forcing_amp = (max_u + min_u) / 2;
%         [t, u, forcing] = rk5_solver(@(t,u)  ode_system_123(t, u, A, B, forcing_amp), tspan, u0);
%         u_square=sqrt(sum(u.^2,2));
%         max_u_square = max(u_square(end-10:end));
% 
% 
%         % Update the interval
%         if max_u_square < 10.^-8
%             min_u = forcing_amp; % Root is in [a, root]
%         else
%             max_u = forcing_amp; % Root is in [root, b]
%         end
%     end
% 
%     % Final check for convergence
%     if (max_u - min_u) > tol
%         warning('Bisection method did not converge within the maximum iterations.');
%     end
% end
% 
% 
% function [forcing_amp, iterations] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations)
%     % Modified Bisection Method with multiple simulations per step
%     for iter = 1:maxIter
%         mid_u = (min_u + max_u) / 2;
%         norm_values = zeros(num_simulations, 1);
% 
%         % Run multiple simulations
%         for sim = 1:num_simulations
%             ini = randn(9, 1);
%             u0_sim = mid_u * ini / norm(ini);
%             [~, u, ~] = rk5_solver(@(t, u) ode_system_123(t, u, A, B, mid_u), tspan, u0_sim);
%             norm_values(sim) = max(sqrt(sum(u.^2, 2)));  % Max norm
%         end
% 
%         % Evaluate average norm
%         avg_norm = mean(norm_values);
% 
%         % Update bounds
%         if avg_norm > 1
%             max_u = mid_u;
%         else
%             min_u = mid_u;
%         end
% 
%         % Check convergence
%         if abs(max_u - min_u) < tol
%             forcing_amp = mid_u;
%             iterations = iter;
%             return;
%         end
%     end
% 
%     % If max iterations reached
%     forcing_amp = (min_u + max_u) / 2;
%     iterations = maxIter;
% end