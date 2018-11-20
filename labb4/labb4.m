im1=imread('dinosaur0.png');
im2=imread('dinosaur5.png');
load('dino_Ps','P');
figure(1);image([im1 im2]);
C1=P{1};       % The camera matrix corresponding to image 1
C2=P{6};       % The camera matrix corresponding to image 2

n1 = null(C1);
n2 = null(C2)
e12 = C1*n2
e21 = C2*n1;
F21 = liu_crossop(e21)*C2*pinv(C1)
%%
[y1,y2]=correspondences_select(im1,im2);
y1(end+1,:)=1;
y2(end+1,:)=1;

%%

l2=F21*y1;
l1=F21'*y2;
for ix=1:size(y1,2),
    l1(:,ix)=-l1(:,ix)/norm(l1(1:2,ix))*sign(l1(3,ix));
    l2(:,ix)=-l2(:,ix)/norm(l2(1:2,ix))*sign(l2(3,ix));
end

figure(2);clf;
subplot(1,2,1);image(im1);hold on;
subplot(1,2,2);image(im2);hold on;
for k=1:size(y1,2),
    subplot(1,2,1);plot(y1(1,k),y1(2,k), 'r*');drawline(l1(:,k));
    subplot(1,2,2);plot(y2(1,k),y2(2,k), 'r*');drawline(l2(:,k));
end

%%
e12_norm = e12/e12(3)
e21_norm = e21/e21(3)

abs(sum(y1.*l1))
abs(sum(y2.*l2))

e12_line = cross(l1(:,1),l1(:,2))
e12_line = e12_line/e12_line(3)

e21_line = cross(l2(:,1),l2(:,2))
e21_line = e21_line/e21_line(3)

%%
w=720;h=576; % Image width and height
lambda = -5:0.01:5; %Replace with the interval you choose
dist = LoopZhangDistortion(e12,e21,w,h,lambda);
figure(3);plot(lambda,dist);
set(gca,'YScale','log');grid on;
lmin = lambda(dist==min(dist)) %Print minimum lambda

%%
lambda = -1.25   % Insert your value here
w1=liu_crossop(e12)*[lambda 1 0]';
w1=w1/w1(3);
w2=F21*[lambda 1 0]';
w2=w2/w2(3);
Hp1=[1 0 0;0 1 0;w1'];
Hp2=[1 0 0;0 1 0;w2'];

%%
vcp=0;
Hr1=[F21(3,2)-w1(2)*F21(3,3) w1(1)*F21(3,3)-F21(3,1) 0;
F21(3,1)-w1(1)*F21(3,3) F21(3,2)-w1(2)*F21(3,3) F21(3,3)+vcp;
0 0 1];
Hr2=[w2(2)*F21(3,3)-F21(2,3) F21(1,3)-w2(1)*F21(3,3) 0;
w2(1)*F21(3,3)-F21(1,3) w2(2)*F21(3,3)-F21(2,3) vcp;
0 0 1];

flip_matrix = [ -1 0 0;
                0 -1 0;
                0 0 1]

H1 = flip_matrix*Hr1*Hp1;
H2 = flip_matrix*Hr2*Hp2;

Fr = inv(H2)'*F21*inv(H1)

%%

oldcorners=[1 w w 1;1 1 h h];
newcorners1=map_points(H1,oldcorners);
newcorners2=map_points(H2,oldcorners);
mincol=min([newcorners1(1,:) newcorners2(1,:)]);
minrow=min([newcorners1(2,:) newcorners2(2,:)]);
T=[1 0 -mincol+1;0 1 -minrow+1;0 0 1]
newcorners1=map_points(T*H1,oldcorners);
newcorners2=map_points(T*H2,oldcorners);
maxcol=max([newcorners1(1,:) newcorners2(1,:)]);
maxrow=max([newcorners1(2,:) newcorners2(2,:)]);
inv_T=inv(diag([maxcol/w maxrow/h 1]))
T=inv(diag([maxcol/w maxrow/h 1]))*T

%%

imr1=image_resample(double(im1),T*H1,h,w);
imr2=image_resample(double(im2),T*H2,h,w);
figure(5);image(uint8([imr1 imr2]));hold on;

y2_trans = map_points(T*H2,y2(1:2,:))
y1_trans = map_points(T*H1,y1(1:2,:))

for k=1:size(y1_trans,2),
    plot(y1_trans(1,k),y1_trans(2,k), 'r*');
    plot(y2_trans(1,k)+720,y2_trans(2,k), 'r*');
    plot([y1_trans(1,k) y2_trans(1,k)+720], [y1_trans(2,k) y2_trans(2,k)])
end
%%

[y1_v1,y2_v1]=correspondences_select(im1,im2);
y1_v1(end+1,:)=1;
y2_v1(end+1,:)=1;
%%
F = fmatrix_n8pa(y1_v1,y2_v1)
%%
figure(6);clf;
subplot(1,2,1);image(im1);hold on;
subplot(1,2,2);image(im2);hold on;

l2_8p=F*y1_v1;
l1_8p=F'*y2_v1;

for k=1:size(y1_v1,2),
    l1_8p(:,k)=-l1_8p(:,k)/norm(l1_8p(1:2,k))*sign(l1_8p(3,k)); %Normalise dual
    l2_8p(:,k)=-l2_8p(:,k)/norm(l2_8p(1:2,k))*sign(l2_8p(3,k)); %homog. coord.
end

for k=1:size(y1_v1,2),
    subplot(1,2,1);plot(y1_v1(1,k),y1_v1(2,k),'or');drawline(l1_8p(:,k));
    subplot(1,2,2);plot(y2_v1(1,k),y2_v1(2,k),'or');drawline(l2_8p(:,k));
end
dist_im1_8p = mean(abs(sum(y1_v1.*l1_8p)))
dist_im2_8p = mean(abs(sum(y2_v1.*l2_8p)))


dist_im1 = mean(abs(sum(y1.*l1)))
dist_im2 = mean(abs(sum(y2.*l2)))

%% Part 5 Triangulation

load('tridata','x','y1','y2','im1','im2','C1','C2');

figure(1);image([im1 im2]);
colormap(gray(256));hold on;
for k=1:size(y1,2),
    plot(y1(2,k),y1(1,k),'ro');
    plot(y2(2,k)+size(im1,2),y2(1,k),'ro');
end
hold off;

xrec=[];
for ix=1:length(y1),
    
    [U,S,V] = svd([liu_crossop(y1(:,ix))*C1; liu_crossop(y2(:,ix))*C2]);
    V_last = V(:,end);
    recpoint = V_last/V_last(4);
    xrec=[xrec recpoint];
end

figure(3);plot3(xrec(1,:),xrec(2,:),xrec(3,:),'o');

er = [];
for ind=1:size(xrec,2),
    curr_er = sqrt(sum((x(:, ind)-xrec(:, ind)).^2));
    er = [er curr_er];
end
x_mean = sum(er)/size(er,2)
x_std = std(er)

%%
y1_rec = C1*xrec;
y2_rec = C2*xrec;

y1_n_rec = [];
y2_n_rec = [];
for i=1:size(xrec,2),
    y1_n = C1*xrec(:, i);
    y1_n_rec = [y1_n_rec y1_n/y1_n(3)];
    y2_n = C2*xrec(:, i);
    y2_n_rec = [y2_n_rec y2_n/y2_n(3)];
end

er= [];
for ind=1:size(y2_n_rec,2),
    curr_er = sqrt(sum( (y2(:, ind)-y2_n_rec(:, ind)).^2));
    er = [er curr_er];
end

proj_mean_y2 = sum(er)/size(er,2)
proj_std_y2 = std(er)


er= [];
for ind=1:size(y1_n_rec,2),
    curr_er = sqrt(sum( (y1(:, ind)-y1_n_rec(:, ind)).^2));
    er = [er curr_er];
end

proj_mean_y1 = sum(er)/size(er,2)
proj_std_y1 = std(er)