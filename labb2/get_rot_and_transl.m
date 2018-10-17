function [ Rest, Test ] = get_rot_and_transl(x1, x2)
[M, N] = size(x1)
a0 = mean(x1,2);
b0 = mean(x2,2);
A = x1 - a0*ones(1,N);
B = x2 - b0*ones(1,N);
Rest = get_rotation_opp(A, B);
Test = b0 - Rest*a0;


end

