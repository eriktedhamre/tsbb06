function [n_hat, alpha] = get_rotation_angle_and_axis(R)
v = liu_crossop((R-transpose(R))/2);
c = (trace(R) - 1)/2;
s = norm(v);
if s == 0
    if c > 0 %  cos(alpha) = 1 && sin(alpha) = 0 => alpha = 0
        arb_rotation_axis = transpose([1, 1, 1]);
        n_hat = arb_rotation_axis/ norm(arb_rotation_axis);
        alpha = 0;
    else % cos(alpha) = -1 && sin(alpha) = 0 => alpha = pi
        max_norm_n = [0, 0, 0];
        for col = [1:3]
            n_k = R + eye(3);
            if norm(n_k) > norm(max_norm_n)
                max_norm_n = n_k;
            end
        end
        n_hat = norm(max_norm_n);
        alpha = pi;
    end
else
    n_hat = v/s;
    alpha = atan(s/c);
end

