load A
size(A)

% C = A*A';                         %Compute the correlation matrix C
% [e l] = eig(C);                   %Compute EVD of C
% [PM p] = sort(diag(l),'decend');  %Sort the eigenvalues: largest first
% PC = e(:,p);                      %Sort the eigenvectors in the same way

[PC S] = svd(A);                  %Compute SVD of A
PM = diag(S);                     %magnitudes are given by the
                                  %singular values
figure(1);
subplot(2,1,1)
plot(PM,'o');
subplot(2,1,2)
plot(log(PM),'o');

M = 4;
figure(2);plot(diff(PM((M + 1):end)),'o');

figure(3);plot(PC(:,1)'*A,PC(:,2)'*A,'o');axis('equal');

%%
singular_squared = PM.^2
e = cumsum(flip(singular_squared))
figure(9);plot([0:99],flip(e),'o');

%%

m=mean(A,2);
A0 = A - m*ones(1,size(A,2));
[PC S] = svd(A);                  
PM = diag(S);
figure(4);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');
figure(5);
subplot(4,1,1);plot(PC(:,1)'*A0,PC(:,2)'*A0,'o');axis('equal')
subplot(4,1,2);plot(PC(:,2)'*A0,PC(:,3)'*A0,'o');axis('equal')
subplot(4,1,3);plot(PC(:,3)'*A0,PC(:,4)'*A0,'o');axis('equal')
subplot(4,1,4);plot(PC(:,4)'*A0,PC(:,5)'*A0,'o');axis('equal')

singular_squared = PM.^2
e = cumsum(flip(singular_squared))
figure(6);plot([0:99],flip(e),'o');

%%
PM

figure(8);plot3(PC(:,1)'*A0, PC(:,2)'*A0, PC(:,3)'*A0, 'o');axis('equal');

%%

im = double(imread('middlebury.png')); % Choose you image here!
size(im)
figure(10);colormap('gray');imagesc(im);

N = 8;
A = im2col(im,[N N],'distinct');
size(A)

[PC S] = svd(A);                 %Compute SVD of A
PM = diag(S.^2);                  %magnitudes are given by the
                                  %singular values
figure(11);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');
%%
figure(12);colormap('gray');
for M=1:6,
    c=PC(:,1:M)'*A;    %Compute coordinates from blocks
    Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
    imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
    subplot(2,3,M);imagesc(imrec);axis('off');     %Display image
    title(sprintf('%d principal components',M));   %Set title
end

%%

figure(14);colormap('gray');
for cnt=1:6,
    subplot(2,3,cnt);h=mesh(reshape(PC(:,cnt),N,N));
    set(h,'edgecolor','black');axis([1 N 1 N -0.5 0.5]);
    title(sprintf('principal component %d',cnt));
end
