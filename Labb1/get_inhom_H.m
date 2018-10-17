function H = get_inhom_H(A)
A_last = A(:, end);
A_rest = A(:, 1:end-1);

z_line = -inv((transpose(A_rest)*A_rest))*transpose(A_rest)*A_last;
z = [z_line; 1];


H = reshape(z,3,3);

end

