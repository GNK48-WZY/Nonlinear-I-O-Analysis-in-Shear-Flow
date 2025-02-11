% addpath(genpath('/home/zhw22003/PRE-Model/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/PRE-Model/mosek/10.2/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/PRE-Model/mosek'))
Re_list = logspace(log10(200), log10(2000), 16);
delta_list = logspace(-6, 0, 200);
% Re_list = logspace(log10(200), log10(2000), 4);
% delta_list = logspace(-6, 0, 2);
% Re_list= Re_list(1:5);
% delta_list=delta_list(1:2);
% epsilon = 0.01;
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
% parpool(32);
forcing_bi = zeros(length(Re_list), length(delta_list));


for ind_Re = 1:length(Re_list)
    Re = Re_list(ind_Re);
    [local_forcing_amp]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear,u,B,C,alpha,Beta,Gamma,KBG, KAG,KABG);
    forcing_bi(ind_Re,:) = local_forcing_amp;
   
end

save Bisection_2_5.mat
function [local_forcing_amp, local_u_bisection, local_mid_u_bisection]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear, u,B,C, alpha,Beta,Gamma,KBG, KAG,KABG)
    local_forcing_amp = NaN*ones(1, length(delta_list)); % 1 x length(delta_list)
    local_u_bisection = cell(1, length(delta_list)); % Cell array of size 1 x length(delta_list)
    local_mid_u_bisection = NaN*ones(1, length(delta_list)); % 1 x length(delta_list)

    % Initialize other variables
    local_u_upper_bound = NaN*ones(1, length(delta_list));
    local_deltaf = NaN*ones(1, length(delta_list));
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
    [lambda, kappa, lambda_af, kappa_ff] = check(A, C, B);
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
        local_u_upper_bound(ind_delta) = (c1*c3*delta)/(c2*c4*L);
        disp(['Feasible solution found for Re = ', num2str(Re), ', delta = ', num2str(delta_list(ind_delta))]);
        disp(['P: ', mat2str(optp)]);
        disp(['umax: ', num2str(umax)]);
        disp(['umin: ', num2str(umin)]);

        T = 20000;
        omega=1;
        [local_forcing_amp, local_u_bisection, local_mid_u_bisection] = compute_norm_u(Re, A, B, T, omega,local_deltaf(ind_delta), local_u_upper_bound(ind_delta), delta_list, alpha,Beta,Gamma,KBG, KAG,KABG);
    else
        disp(['No feasible solution found for Re = ', num2str(Re), ', delta = ', num2str(delta_list(ind_delta))]);
        disp(['YALMIP error: ', yalmiperror(result_yalmip.problem)]);
    end
end

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

function [local_forcing_amp, local_u_bisection, local_mid_u_bisection] = compute_norm_u(Re, A, B, T, omega, deltaf, u_upper_bound, delta_list, alpha,Beta,Gamma,KBG, KAG,KABG)


ini=randn(9,1);
u0 = deltaf*ini/norm(ini);
% u0 =zeros(9,1);

T = 2000;
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
[t, u, forcing] = rk5_solver(@(t,u)  ode_system_123(t, u, A, B, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG), tspan, u0);

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


% %Bisection
max_u = 0.01; 
min_u = 0;  
% Tolerance and maximum iterations
tol = 1e-14;
maxIter = 100;

% Call the bisection method
% u0=zeros(9,1);
% [forcing_amp, iterations] = bisectionMethod(min_u, max_u, tol, maxIter, A, B,tspan, u0);
num_simulations = 20;  % Number of simulations per step
[local_forcing_amp, iterations, local_u_bisection, local_mid_u_bisection] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations,deltaf,delta_list,alpha,Beta,Gamma,KBG, KAG,KABG);
% Display results
fprintf('Approximated root: %.20f\n', local_forcing_amp);
fprintf('Number of iterations: %d\n', iterations);
% 
% %% t vs local_mid_u_bisection
% t = linspace(0, 10, 40); % Creates a 1x40 time vector from 0 to 10
% loglog(t, local_mid_u_bisection);
% xlabel('Time');
% ylabel('Mid-u Values');
% title('Log-Log Plot of Mid-u vs Time');
% 
% %% t vs local_u_bisection
% % Initialize a vector to store the maximum values
% max_values = zeros(40, 1); % 40 iterations
% 
% % Loop through each iteration
% for iter = 1:40
%     % Initialize a temporary vector to store max values for this iteration
%     temp_max = zeros(20, 1); % 20 simulations per iteration
%     
%     % Loop through each simulation
%     for sim = 1:20
%         data = local_u_bisection{iter, sim}; % Extract data from the cell
%         
%         % Compute the maximum value for this simulation
%         % Use max(data(:)) to ensure a scalar is returned
%         temp_max(sim) = max(data(:));
%     end
%     
%     % Store the maximum value across all simulations for this iteration
%     max_values(iter) = max(temp_max);
% end
% 
% % Create a time vector (iteration numbers)
% t = 1:40;
% 
% % Plot the maximum values
% loglog(t, max_values, '-o'); % Plot with circles for each point
% xlabel('Iteration');
% ylabel('Maximum Value');
% title('Log-Log Plot of Maximum Values Across Simulations');
% grid on;

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
function [local_forcing_amp, iterations, local_u_bisection, local_mid_u_bisection] = bisectionMethodMulti(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations, deltaf, delta_list, alpha,Beta,Gamma,KBG, KAG,KABG)
    local_u_bisection = cell(maxIter, num_simulations); % Preallocate cell array
    local_mid_u_bisection = NaN * ones(maxIter, 1); % Preallocate array for mid_u values

    % Modified Bisection Method with multiple simulations per step
    for iter = 1:maxIter
        mid_u = (min_u + max_u) / 2;
        norm_values = zeros(num_simulations, 1);

        % Run multiple simulations
        for sim = 1:num_simulations
            ini = randn(9, 1);
            u0_sim = deltaf * ini / norm(ini);
%             u0_sim = zeros(9, 1); % Reset initial condition
            [~, u, ~] = rk5_solver(@(t, u) ode_system_123(t, u, A, B, mid_u, alpha,Beta,Gamma,KBG, KAG,KABG), tspan, u0_sim);
            
            % Store u for this simulation and iteration
            local_u_bisection{iter, sim} = u;

            % Calculate norm only for late time (last 50% of the time domain)
            late_time_indices = round(0.9 * length(tspan)):length(tspan);
            norm_values(sim) = mean(sqrt(sum(u(late_time_indices, :).^2, 2)));  % Average norm over late time
%             norm_values_timehistory = sqrt(sum(u(late_time_indices, :).^2, 2));  % t over norm_values_timehistory

        end

        % Store mid_u for this iteration
        local_mid_u_bisection(iter) = mid_u;

        % Evaluate average norm
        avg_norm = mean(norm_values);

        % Update bounds
        if avg_norm > 10.^-8
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

    % If max iterations reached
    if iter == maxIter
        local_forcing_amp = (min_u + max_u) / 2;
        iterations = maxIter;
    end

%     % Save mid_u values for analysis
%     save('local_mid_u_bisection.mat', 'local_mid_u_bisection');
% 
%     % Separate data into two groups based on local_forcing_amp
     above_threshold = local_mid_u_bisection > local_forcing_amp;
     below_threshold = local_mid_u_bisection <= local_forcing_amp;
%     
%     save('local_u_bisection.mat',"local_u_bisection", "below_threshold", "above_threshold");
% 
    % Plot |u| as a function of time for both groups
    colors = lines(16);

    figure;
    hold on;
    for iter = 1:maxIter
        for sim = 1:num_simulations
            u = local_u_bisection{iter, sim};
            if above_threshold(iter)
            plot(tspan(late_time_indices), sqrt(sum(u(late_time_indices, :).^2, 2)), 'Color', colors(1, :)); % Above threshold in red
            end
        end
    end
    xlabel('Time');
    ylabel('|u|');
    title('|u| as a Function of Time');
    legend('Above Threshold');
    hold off;
    

    figure;
    hold on;
    for iter = 1:maxIter
        for sim = 1:num_simulations
            u = local_u_bisection{iter, sim};
            if below_threshold(iter)
            plot(tspan(late_time_indices), sqrt(sum(u(late_time_indices, :).^2, 2)), 'Color', colors(1, :)); % Above threshold in red
            end
        end
    end
    xlabel('Time');
    ylabel('|u|');
    title('|u| as a Function of Time');
    legend('Below Threshold');
    hold off;
    
end

function [local_forcing_amp, iterations, local_u_bisection, local_mid_u_bisection] = bisectionMethodMulti2(min_u, max_u, tol, maxIter, A, B, tspan, u0, num_simulations, deltaf, delta_list, alpha,Beta,Gamma,KBG, KAG,KABG)
    local_u_bisection = cell(length(maxIter), length(num_simulations)); % Preallocate cell array
    local_mid_u_bisection = NaN * ones(length(maxIter), 1); % Preallocate array for mid_u values

    % Modified Bisection Method with multiple simulations per step
    for iter = 1:maxIter
        mid_u = (min_u + max_u) / 2;
        norm_values = zeros(num_simulations, 1);

        % Run multiple simulations
        for sim = 1:num_simulations
            ini = randn(9, 1);
            u0_sim = deltaf * ini / norm(ini);
            u0_sim = 0*zeros(9,1);
            [~, u, ~] = rk5_solver(@(t, u) ode_system_123(t, u, A, B, mid_u, alpha,Beta,Gamma,KBG, KAG,KABG), tspan, u0_sim);
            norm_values(sim) = max(sqrt(sum(u.^2, 2)));  % Max norm

            % Store u for this simulation and iteration
            local_u_bisection{iter, sim} = u;
        end

        % Store mid_u for this iteration
        local_mid_u_bisection(iter) = mid_u;

        % Evaluate average norm
        avg_norm = mean(norm_values);

        % Update bounds
        if avg_norm > 10.^-8
            max_u = mid_u;
        else
            min_u = mid_u;
        end

        % Check convergence
        if abs(max_u - min_u) < tol
            local_forcing_amp = mid_u;
            iterations = iter;
            return;
        end
    end

    % If max iterations reached
    local_forcing_amp = (min_u + max_u) / 2;
    iterations = maxIter;
    save('local_mid_u_bisection.mat', 'local_mid_u_bisection');

end