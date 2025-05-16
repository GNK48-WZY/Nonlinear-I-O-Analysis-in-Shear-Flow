
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
