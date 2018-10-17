function q = get_unit_quaternion(n, alpha)
s_part = cos(alpha/2);
v_part = sin(alpha/2) * n;
q = [s_part; v_part];

end

