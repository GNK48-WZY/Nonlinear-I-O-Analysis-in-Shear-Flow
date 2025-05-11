function [norm_u, u, G, sys, max_SVD, t, local_forcing, y, dt, a_square] = compute_norm_u(Re, A, B, T, omega, deltaf, u_upper_bound, alpha,Beta,Gamma,KBG, KAG,KABG)


ini=randn(9,1);
u0 = deltaf*ini/norm(ini);
% u0 =zeros(9,1);

T = 70000;
tspan = linspace(0, T, 1000000);

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