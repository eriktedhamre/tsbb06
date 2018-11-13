%% Signal for part 1

k = 0:1:100;
s = sin(k/10);
figure(1);plot(s);
title('Signal (s) sin(k)')

%% Signal for part 4

cert = double(rand(1,101)>0.2);
scert = s.*cert;
figure(10);plot(scert);
title('Signal with uncertainty (scert) sin(k)')
s = scert;

%%
figure(6);plot(s);
title('Current signal')
x=(-3:3)';
b0=ones(7,1);
b1=x;
b2=x.^2;

% --- PLOTS
figure(2);
subplot(4,1,1);plot(b0,'-o');
title('Plot of b0')
subplot(4,1,2);plot(b1,'-o');
title('Plot b1')
subplot(4,1,3);plot(b2,'-o');
title('Plot b2')
% ----

a = exp(-x.^2/4);
figure(2);
subplot(4,1,4);plot(a,'-o');
title('Plot of a')
%%
f0 = b0.*a; f0 = f0(end:-1:1);
f1 = b1.*a; f1 = f1(end:-1:1);
f2 = b2.*a; f2 = f2(end:-1:1);
figure(3);
subplot(3,1,1);plot(f0,'-o');
title('f0')
subplot(3,1,2);plot(f1,'-o');
title('f1')
subplot(3,1,3);plot(f2,'-o');
title('f2')

h0 = conv(s,f0,'same');
h1 = conv(s,f1,'same');
h2 = conv(s,f2,'same');

G0 = diag(a)
B = [b0 b1 b2];
G = B'*G0*B
%% Considers uncertainty
%a_cert = a'.* cert(1:7)
%G0_cert = diag(a_cert)
%B_cert = [b0 b1];
%G_cert = B_cert'*G0_cert*B_cert


G11 = conv(cert, flip(b0.*a.*b0),'same')
G12 = conv(cert, flip(b0.*a.*b1),'same')
G22 = conv(cert, flip(b1.*a.*b1),'same')

detG = G11.*G22-G12.^2;
c0 = (G22.*h0-G12.*h1)./detG
c1 = (-G12.*h0+G11.*h1)./detG 
figure(21);
subplot(2,1,1);plot(c0)
title('C with uncertainty')
subplot(2,1,2);plot(c1)

%% Creates c without consideration to uncertainty
c = inv(G)*[h0;h1;h2];
figure(8);

% Coordinates
subplot(3,1,1);plot(c(1,:))
title('s(x)')
subplot(3,1,2);plot(c(2,:))
title('ds(x)/d(x)')
subplot(3,1,3);plot(c(3,:))
title(' ')

figure(7);
localsig=s(60-3:60+3);
reconsig=(B*c(:,60))';
diffsig=localsig-reconsig;

subplot(3,1,1);plot(localsig);
title('Real signal')
subplot(3,1,2);plot(reconsig);
title('Reconstructed signal')
subplot(3,1,3);plot(diffsig);
title('Diff')


%% Test For preperation ex nr 3. (DONT RUN AGAIN) 
f = inv(G)*[f0'; f1'; f1']
h0_t = conv(s,f(1,:),'same');
h1_t = conv(s,f(2,:),'same');
h2_t = conv(s,f(3,:),'same');
figure(3)
subplot(3,1,1);plot(h0_t)
subplot(3,1,2);plot(h1_t)
subplot(3,1,3);plot(h2_t)
%% 

im = double(imread('Scalespace0.png'));
figure(8);colormap(gray);imagesc(im);title('Original');

cert = double(rand(size(im)) > 0.60); imcert = im.*cert;
figure(9);colormap(gray);imagesc(imcert);title('Image * cert')

x = ones(7,1)*(-3:3)
y = x'
%a = exp(-(x.^2+y.^2)/4);
%a = exp((0.7./(log(x.^2 + y.^2))+0.6)/10);
a = exp((-(x.^2+y.^2))/4);
figure(10);mesh(a);

imlp = conv2(imcert, a, 'same');
figure(11);colormap(gray);imagesc(imlp); title('lowpass')

G = conv2(cert, flip(a), 'same');
c = imlp./G;
figure(14);colormap(gray);imagesc(c);title('Reconstructed');
figure(13);colormap(gray);imagesc(G);title('G');
