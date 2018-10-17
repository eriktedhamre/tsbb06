function R = get_rod_rotation_matrix(n, alpha)
R = eye(3) + (1 - cos(alpha))*liu_crossop(n)^2+sin(alpha)*liu_crossop(n);
end

