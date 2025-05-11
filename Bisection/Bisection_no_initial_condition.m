%% Bisection method to get the critical forcing for no initial condition
%% Date: 2025/5/8
%% Outputs: local_forcing_amp: critical forcing without initial condition
%% Important condition: for in figure 9, the different \Pi means that change the avg_norm > 1e-10(e-9, e-8)


% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/10.2/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))
Re_list = logspace(log10(200), log10(2000), 16);
% delta_list = logspace(-6, 0, 200);
% delta_list = 1e-6;
% Re_list = logspace(log10(200), log10(2000), 4);
% delta_list = logspace(-6, 0, 2);
% Re_list= Re_list(1:5);
% delta_list=delta_list(1:2);
% epsilon = 0.01;
% Re_list = 200;
% delta_list = 1e-6;
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
[RHS_J_mean_shear, nonlinear, u] = nonliner(alpha,Beta,Gamma,KBG, KAG,KABG);
delete(gcp('nocreate'));
parpool(4);
% forcing_bi = zeros(length(Re_list), length(delta_list));


parfor ind_Re = 1:length(Re_list)
    Re = Re_list(ind_Re);
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
    [local_forcing_amp] = compute_norm_u(Re, A, B, alpha,Beta,Gamma,KBG, KAG,KABG);   
end

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

function [local_forcing_amp] = compute_norm_u(Re, A, B, alpha,Beta,Gamma,KBG, KAG,KABG)


ini=randn(9,1);
u0 = deltaf*ini/norm(ini);
% u0 =zeros(9,1);

T = 70000;
tspan = linspace(0, T, 1000000);

I = eye(size(A));
C = eye(9);
D = zeros(9, 9);

% G = C * inv(1i * omega * I - A) * B;
% max_SVD = max(svd(G));
% sys = ss(A, B, C, D);
% u0 = zeros(9,1);
% tol = 1e-9;  % Tolerance for convergence
% max_iter = 1000;  % Maximum fixed-point iterations per step
% dt = diff(tspan);
% dt = dt(1);
%[t, u, forcing] = ode45(@(t,u)  ode_system_123(t, u, A, B, u_upper_bound), tspan, u0);
% [t, u, forcing] = rk5_solver(@(t,u)  ode_system_123(t, u, A, B, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG), tspan, u0);
% Linear Equation analytical solution
% for t_ind = 1:length(tspan)
%     syms t_sym real; % Create a symbolic variable for time
%     u_sym = expm(A *tspan(t_ind)) * u0; % Analytical solution: u(t)=exp(At)*u0
%     u_analytical(:,t_ind) = double(u_sym); % Store symbolic expression
% end 
% %old norm integrate over time
% dt=diff(tspan);
% dt=dt(1);
% % u_square=sum(u.^2,2);
% % norm_u=sqrt(sum(u_square(2:end).*dt)/max(tspan));
% 
% %1 norm maximal over time
% a_square=sqrt(sum(u.^2,2));
% %norm_u = sqrt(sum(u_square.^2)*dt);
% y =C*u';
% y = y'; % output
% norm_u=max(a_square);


%% Bisection method
max_u = 1e-3; 
min_u = 0;  
% Tolerance and maximum iterations
tol = 1e-20;
maxIter = 300;

% Call the bisection method
% u0=zeros(9,1);
% [forcing_amp, iterations] = bisectionMethod(min_u, max_u, tol, maxIter, A, B,tspan, u0);
num_simulations = 5;  % Number of simulations per step
[local_forcing_amp, iterations,norm_history, mid_u_history] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations,alpha,Beta,Gamma,KBG, KAG,KABG);
% Display results
fprintf('Approximated root: %.20f\n', local_forcing_amp);
fprintf('Number of iterations: %d\n', iterations);
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

function [local_forcing_amp, iterations, norm_history, mid_u_history] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, ~, num_simulations, alpha,Beta,Gamma,KBG, KAG,KABG)
    local_u_bisection = cell(maxIter, num_simulations); % Preallocate cell array
    mid_u_history = NaN * ones(maxIter, 1); % Preallocate array for mid_u values

    % Modified Bisection Method with multiple simulations per step
    for iter = 1:maxIter
        mid_u = (min_u + max_u) / 2;
        norm_values = zeros(num_simulations, 1);

        % Run multiple simulations
        for sim = 1:num_simulations
            ini = randn(9, 1);
%             u0_sim = deltaf * ini / norm(ini);
            u0_sim = zeros(9, 1); % Reset initial condition
            [~, u, ~] = rk5_solver(@(t, u) ode_system_123(t, u, A, B, mid_u, alpha, Beta, Gamma, KBG, KAG, KABG), tspan, u0_sim);
            
            % Store u for this simulation and iteration
            local_u_bisection{iter, sim} = u;

            % Calculate norm only for late time
            late_time_indices = round(0.9 * length(tspan)):length(tspan);
            late_time_indices = 1:length(tspan);

            norm_values(sim) = mean(sqrt(sum(u(late_time_indices, :).^2, 2)));  % Average norm over late time
%             norm_values_timehistory = sqrt(sum(u(late_time_indices, :).^2, 2));  % t over norm_values_timehistory

        end

        % Store mid_u for this iteration
        mid_u_history(iter) = mid_u;
        avg_norm = mean(norm_values);
        norm_history(iter) = avg_norm;
        % Evaluate average norm
        fprintf('Iter %3d: mid forcing = %.3e, avg norm = %.3e\n', iter, mid_u, avg_norm);

        % Update bounds
        if avg_norm > 1e-10
            max_u = mid_u;
        else
            min_u = mid_u;
        end

        % Check convergence
        if abs(max_u - min_u) < tol
            local_forcing_amp = mid_u;
            iterations = iter;
            break;
        end

        
    end

%     % If max iterations reached
%     if iter == maxIter
%             local_forcing_amp = (min_u + max_u) / 2;
%             iterations = maxIter;
%     end

    local_forcing_amp = (min_u + max_u) / 2;
%     iterations = maxIter;

%     % Save mid_u values for analysis
%     save('local_mid_u_bisection.mat', 'local_mid_u_bisection');
% 
%     % Separate data into two groups based on local_forcing_amp
     above_threshold = mid_u_history > local_forcing_amp;
     below_threshold = mid_u_history <= local_forcing_amp;
%     
%     save('local_u_bisection.mat',"local_u_bisection", "below_threshold", "above_threshold");
% 
%     Plot |u| as a function of time for both groups

% % Above threshold
%     figure;
%     axes('XScale', 'log', 'YScale', 'log')
%     hold on;
%     for iter = 1:maxIter
%         for sim = 1:num_simulations
%             u = local_u_bisection{iter, sim};
%             if above_threshold(iter)
%             loglog(tspan(late_time_indices), sqrt(sum(u(late_time_indices, :).^2, 2))); % Above threshold in red
%             hold on;
%             end
%         end
%     end
% %     xlim([1, 10000])
%     xlabel('Time($\mathbf{t}$)','Interpreter', 'latex');
%     ylabel('$\|\mathbf{a}\|$','Interpreter', 'latex');
%     hold off;
% 
% % Below threshold
%     figure;
%     axes('XScale', 'log', 'YScale', 'log')
%     hold on;
%     for iter = 1:maxIter
%         for sim = 1:num_simulations
%             u = local_u_bisection{iter, sim};
%             if below_threshold(iter)
%             loglog(tspan(late_time_indices), sqrt(sum(u(late_time_indices, :).^2, 2))); % below threshold
%             hold on;
%             end
%         end
%     end
% %     xlim([1, 10000])
%     xlabel('Time($\mathbf{t}$)','Interpreter', 'latex');
%     ylabel('$\|\mathbf{a}\|$','Interpreter', 'latex');
%     hold off;
end