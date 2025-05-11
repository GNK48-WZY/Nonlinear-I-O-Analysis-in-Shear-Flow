function[result_yalmip, P,lambda_af,kappa_ff] = LMI(A,s_bound, epsilon, lambda_af, kappa_ff, diag_s, B, C, s, sdp_options)
%% Solve the LMI
%% diag_s and s are from the s_diag_s.m, lambda_af and kappa_ff from check.m


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