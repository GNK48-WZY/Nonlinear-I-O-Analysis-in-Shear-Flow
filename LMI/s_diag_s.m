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