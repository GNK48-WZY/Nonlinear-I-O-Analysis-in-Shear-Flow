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

            
             upper_bound_y = local_Gamma_theorem(ind_delta)*local_u_upper_bound(ind_delta) +beta_theorem_inf*ones(size(t));
            figure;
%             Plot the norm from simulation
%             loglog(t, u_square, 'LineWidth', 1.5, 'DisplayName', 'norm usimulation'); hold on;
            loglog(t, y_norm, 'LineWidth', 1.5, 'DisplayName', 'norm y simulation'); hold on;

            % Plot the horizontal line for the upper bound
            loglog(t,upper_bound_y, 'DisplayName', 'upper bound');
            %loglog(t, local_u_upper_bound(ind_delta) * ones(size(t)), '--r', 'LineWidth', 1.5, 'DisplayName', 'Upper Bound');
            % Plot Lp norms as horizontal lines
%             loglog(t, u_tau_Lp(ind_p) * ones(size(t)), '--', 'DisplayName', ['Lp Norm (p=', num2str(p_list(ind_p)), ')']); % Lp_norm equal to zero after p=20, may be is numerical problem
%             loglog(t, u_tau_Linf * ones(size(t)), '-.', 'DisplayName', 'L∞ Norm'); % Plot Linf norm
            xlabel('Time (t)');
            ylabel('Norm Values');
            title(['Norm Analysis for \delta = ', num2str(delta_list(ind_delta))]);
            legend('show');
            grid on;
            hold off;

 
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
%         colors = lines(16);
%         line_styles = {"-","--","-.",":",""};
%         markers = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
%         figure;
%         loglog(p_list, local_y_tau_Lp, '--r', 'LineWidth',1.5,'Marker',markers{1},'Color', colors(1, :),'DisplayName','LHS of Inequality (p=1-10^4)'); hold on;
%         loglog(p_list, local_Gamma_theorem(ind_delta) .*local_u_tau_Lp + local_beta_theorem_p, '-.', 'LineWidth',1.5,'Color', colors(2, :), 'Marker',markers{2},'DisplayName', 'RHS of Inequality (p=1-10^4)');
%         loglog(p_list, y_tau_Linf.*ones(size(p_list)),'-r','Marker',markers{3}, 'LineWidth',1.5,'Color', colors(3, :),'DisplayName', 'LHS of Inequality (p=\infty)');
%         loglog(p_list, (local_Gamma_theorem(ind_delta) .* u_tau_Linf + beta_theorem_inf).*ones(size(p_list)),'-k', 'LineWidth',1.5,'Color', colors(4, :),'Marker',markers{4},'DisplayName' ,'RHS of Inequality (p=\infty)');
%         xlabel('p', 'Interpreter', 'latex');
% %         title('Gap between LHS and RHS change over p');
%         set(gca, 'FontSize', 13);
%         legend('show');
%         grid on;
%         hold off;


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