%% Signal for part 1

k = 0:1:100;
s = sin(k/10);
figure(1);plot(s);
%% Signal for part 4

cert = double(rand(1,101)>0.2);
scert = s.*cert;
figure(6);plot(scert);
s = scert

%%
figure(6);plot(s);
x=(-3:3)';
b0=ones(7,1);
b1=x;
b2=x.^2;
figure(2);
subplot(4,1,1);plot(b0,'-o');
subplot(4,1,2);plot(b1,'-o');
subplot(4,1,3);plot(b2,'-o');

a = exp(-x.^2/4);
figure(2);
subplot(4,1,4);plot(a,'-o');

f0 = b0.*a; f0 = f0(end:-1:1);
f1 = b1.*a; f1 = f1(end:-1:1);
f2 = b2.*a; f2 = f2(end:-1:1);
figure(3);
subplot(3,1,1);plot(f0,'-o');
subplot(3,1,2);plot(f1,'-o');
subplot(3,1,3);plot(f2,'-o');

h0 = conv(s,f0,'same');
h1 = conv(s,f1,'same');
h2 = conv(s,f2,'same');

G0 = diag(a)
B = [b0 b1 b2];
G = B'*G0*B
a'
cert(1:7)'
a_cert = a' * cert(1:7)'
G0_cert = diag(a_cert)
B_cert = [b0 b1];
G_cert = B_cert'*G0_cert*B_cert
%G11 = conv(b1);
%G12 = conv();
%G22 = conv();

%%
c = inv(G)*[h0;h1;h2];
figure(8);
subplot(3,1,1);plot(c(1,:))
subplot(3,1,2);plot(c(2,:))
subplot(3,1,3);plot(c(3,:))

figure(7);
localsig=s(60-3:60+3);
reconsig=(B*c(:,60))';
diffsig=localsig-reconsig;
subplot(3,1,1);plot(localsig);
subplot(3,1,2);plot(reconsig);
subplot(3,1,3);plot(diffsig);


%% For preperation ex nr 3. (DONT RUN AGAIN) 
f = inv(G)*[f0'; f1'; f1']
h0_t = conv(s,f(1,:),'same');
h1_t = conv(s,f(2,:),'same');
h2_t = conv(s,f(3,:),'same');
figure(3)
subplot(3,1,1);plot(h0_t)
subplot(3,1,2);plot(h1_t)
subplot(3,1,3);plot(h2_t)
%% 

