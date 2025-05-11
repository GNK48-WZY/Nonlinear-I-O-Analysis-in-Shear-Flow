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