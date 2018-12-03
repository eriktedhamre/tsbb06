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
