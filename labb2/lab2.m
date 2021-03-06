%% Variables

n = transpose([1, 2, 0]);
n = n/norm(n)
alpha = 0.3


%% 4.1 Rotation axis, rotation angle, and rotation matrix



R_rodr1 = get_rod_rotation_matrix(n, alpha)

rotated_n = R*n
[extracted_n extracted_alpha] = get_rotation_angle_and_axis(R_rodr1)


%% 4.2 Matrix exponentials
R_exp = get_exp_rotation_matrix(n, alpha)
R_rodr2 = liu_rodrigues(n, alpha)

[vectors1 values1] = eig(R_rodr2);
[vectors2 values2] = eig(alpha*liu_crossop(n));

vectors1
vectors2

values1
expm(values2)


%% 4.3 Unit quaternions

rotation_q = get_unit_quaternion(n, alpha);

x0 = randn(3,1)
p = [0;x0];
rotated1_x0 = R_rodr2*x0
temp = liu_qmult(rotation_q, p);
rotated2_x0 = liu_qmult(temp, liu_qconj(rotation_q))

R_liu = liu_R_from_q(rotation_q)
R_rodr2

%% 5.1 Creating Synthetic data

N = 9;
x1 = 2*rand(3, N)-1;
t = 2*rand(3,1)-1;      % translation vector
n = 2*rand(3,1)-1;

n = n/norm(n);          % normalized rotation axis
a = 2*pi*rand(1,1);     % rotation angle
R = liu_rodrigues(n,a); % R from (n,a)

x2 = R*x1+t*ones(1,N);

s = 0.01;
x1n = x1 + s*randn(3,N);
x2n = x2 + s*randn(3,N);
%% 5.2 Estimating R and t


[Rest, Test] = get_rot_and_transl(x1n, x2n)
determinant_of_R = det(Rest)
R
t
x2e = Rest*x1n + Test*ones(1, N);
err = norm(x2e - x2n, 'fro')

%% 5.2 part 2
%load rigiddataA
[RdataA, TdataA] = get_rot_and_transl(data1, data2) 
determinant_of_R = det(RdataA)

[M, N] = size(data1)
[ndataA, alphadataA] = get_rotation_angle_and_axis(RdataA) 

est_data2A = RdataA*data1 + TdataA*ones(1, N);
err = norm(est_data2A - data2, 'fro')
%% Part 2 b
%load rigiddataB
[RdataB, TdataB] = get_rot_and_transl(data1, data2) 
determinant_of_R = det(RdataB)

[M, N] = size(data1);
[ndataB, alphadataB] = get_rotation_angle_and_axis(RdataB);

est_data2B = RdataB*data1 + TdataB*ones(1, N);
err = norm(est_data2B - data2, 'fro')
