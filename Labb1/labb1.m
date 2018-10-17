

draw_images()

y1 = vgg_get_homg([516 344; 780 36; 43 472; 236 20; 456 490; 217 239; 44 283; 266 448]');
y2 = vgg_get_homg([476 347; 570 41; 146 588; 168 123; 476 490; 217 325; 87 409; 319 503]');

figure(1);hold('on');plot(y1(1,:),y1(2,:),'go')
figure(2);hold('on');plot(y2(1,:),y2(2,:),'go')

A = make_A_matrix(8, y1, y2);

A0 = make_A_matrix(4, y1, y2)

H1 = get_inhom_H(A)
[U S V] = svd(A);
H2 = reshape(V(:,end),3,3)

%e1_hom = get_geometric_error(H2, y1, y2)

%img2t=image_resample(double(img1),H1,640,800);
%figure(3);imagesc(uint8(img2t))
%imagesc(uint8(img2t)-img2)

diag(S)'
figure(4);plot(log(diag(S)),'o');

[y1tilde T1]= liu_preconditioning(y1)
[y2tilde T2]= liu_preconditioning(y2)

Atilde = make_A_matrix(8, y1tilde, y2tilde);
[Utilde Stilde Vtilde] = svd(Atilde);
Htilde = reshape(Vtilde(:,end),3,3)

H3 = inv(T2)*Htilde*T1;

diag(Stilde)'
figure(4);plot(log(diag(Stilde)),'o');
%e3 = get_geometric_error(H3, y1, y2)

img2t=image_resample(double(img1),H3,640,800);
figure(3);imagesc(uint8(img2t))

e4 = get_H_diffs(H1to2p, H1, y1, y2)
e5 = get_H_diffs(H1to2p, H2, y1, y2)
e6 = get_H_diffs(H1to2p, H3, y1, y2)

l1 = cross(y1(:,1), y1(:,2));
figure(1);drawline(l1,'axis','xy');

l2 = inv(H1')*l1;
figure(2);drawline(l2,'axis','xy');


yl1 = vgg_get_homg([516 344; 780 36; 691 140; 362 524; 321 572 ]');
yl2 = vgg_get_homg([476 347; 570 41; 540 140; 417 544; 401 597]');

[y1ltilde Tl1]= liu_preconditioning(yl1)
[y2ltilde Tl2]= liu_preconditioning(yl2)

Al = make_A_matrix(5, yl1, yl2);
[Ul Sl Vl] = svd(Al);
Hl = reshape(Vl(:,end),3,3)

Altilde = make_A_matrix(5, y1ltilde, y2ltilde);
[Ultilde Sltilde Vltilde] = svd(Altilde);
Hltilde = reshape(Vltilde(:,end),3,3);
Hltilde1 = inv(Tl2)*Hltilde*Tl1
img2t=image_resample(double(img1),Hl,640,800);
figure(5);imagesc(uint8(img2t));

img3t=image_resample(double(img1),Hltilde1,640,800);
figure(6);imagesc(uint8(img2t));


%%
A_last = curr_A(:, end);
A_rest = curr_A(:, 1:end-1);

z_line = -inv((transpose(A_rest)*A_rest))*transpose(A_rest)*A_last;
z = [z_line; 1];


H1 = reshape(z,3,3);

img2t=image_resample(double(img1),H1,640,800);
figure(3);imagesc(uint8(img2t))
%imagesc(uint8(img2t)-img2)

%%
y2b = vgg_get_nonhomg(H1*y1);
y1b = vgg_get_nonhomg(inv(H1)*y2);
figure(1);plot(y1b(1,:),y1b(2,:),'rx');
figure(2);plot(y2b(1,:),y2b(2,:),'rx');

e1 = 0;
for k = 1:length(y1),
    e1 = e1 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2 + ...
    + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end
e1 = sqrt(e1)

%%  Unsymmetric error
e2 = 0;
for k = 1:length(y1),
    e2 = e2 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2;
end
e2 = sqrt(e2)

e3 = 0;
for k = 1:length(y1),
    e3 = e3 + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end
e3 = sqrt(e3)
%% Setting first elem in z to 1

img1=imread('img1.ppm');
img2=imread('img2.ppm');
figure(1);imagesc(img1);
figure(2);imagesc(img2);

y1 = vgg_get_homg([516 344; 780 36; 43 472; 236 20; 456 490; 217 239; 44 283; 266 448]');
y2 = vgg_get_homg([476 347; 570 41; 146 588; 168 123; 476 490; 217 325; 87 409; 319 503]');

figure(1);hold('on');plot(y1(1,:),y1(2,:),'go')
figure(2);hold('on');plot(y2(1,:),y2(2,:),'go')



A = [];
N = 8;
for cnt=1:N,
    coord1 = y1(:, cnt);
    coord2 = y2(:, cnt);
    A = [A;
        [coord1(1) 0 -coord1(1)*coord2(1)  coord1(2) 0 -coord1(2)*coord2(1) 1 0 -coord2(1)];
        [0  coord1(1)  -coord1(1)*coord2(2) 0 coord1(2) -coord1(2)*coord2(2) 0 1 -coord2(2)]];
end

A0 = [];
M = 4;
for cnt=1:M,
    coord1 = y1(:, cnt);
    coord2 = y2(:, cnt);
    A0 = [A0;
        [coord1(1) 0 -coord1(1)*coord2(1)  coord1(2) 0 -coord1(2)*coord2(1) 1 0 -coord2(1)];
        [0  coord1(1)  -coord1(1)*coord2(2) 0 coord1(2) -coord1(2)*coord2(2) 0 1 -coord2(2)]];
end

curr_A = A;

A_first= curr_A(:, 1);
A_rest = curr_A(:, 2:end);

z_line = -inv((transpose(A_rest)*A_rest))*transpose(A_rest)*A_first;
z = [1; z_line];


H1 = reshape(z,3,3);

img2t=image_resample(double(img1),H1,640,800);
figure(3);imagesc(uint8(img2t))
%imagesc(uint8(img2t)-img2)


