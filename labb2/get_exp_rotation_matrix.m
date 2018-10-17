function R = get_exp_rotation_matrix(n, alpha)
R = expm(alpha * liu_crossop(n));
end

