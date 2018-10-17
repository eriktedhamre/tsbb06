function R = get_rotation_opp(A,B)
[U, S, V] = svd(A*B');
R1 = V*U';
if det(R1) < 0
    disp('---- R1 not in SO(3)! ----')
    tau = sign(det(U)*det(V));
    sigma = [1 ,0, 0; 0, 1, 0; 0, 0 , tau];
    R = V*sigma*U';
else
    R = R1;
end
end

