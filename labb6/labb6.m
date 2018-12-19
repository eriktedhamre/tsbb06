%% 3.1 Down-sampling and reconstruction of a discrete signal

s=randn(1024,1);
%%
S=fftshift(fft(s));
u=(-512:511)'*pi/512;
Rect8=(abs(u)<pi/8);
Sbl=S.*Rect8; %Create a pi/4-band-limited transform

sbl=ifft(ifftshift(Sbl));

figure(1);subplot(2,1,1);plot(0:1023,sbl);
title('band-limited signal');
subplot(2,1,2);plot(u,abs(Sbl));
title('Fourier transform of band-limited signal');

downsampl8=sbl(1:8:end);

upsampl8=zeros(1024,1);
upsampl8(1:8:end)=downsampl8;
figure(2);subplot(2,1,1);plot(0:1023,upsampl8);
title('The signal down-sampled with a factor 8');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl8))));
title('Fourier transform of down-sampled signal');

%%

rect8=ifft(fftshift(8*Rect8));
figure(3);
subplot(2,1,1);plot(-512:511,ifftshift(rect8));
title('Reconstructing filter rect8');
subplot(2,1,2);plot(u,abs(fftshift(fft(rect8))));

%%

for ix=0:127,
B8(:,ix+1)=circshift(rect8,8*ix);
end
figure(4);plot(0:1023,B8(:,[1 2 3]));
title('Reconstructing functions 3.1')

%%

sblrec8=B8*upsampl8(1:8:end);
figure(5);
subplot(2,1,1);plot(0:1023,sblrec8);
title('Reconstructed signal from factor 8 down-sampled signal');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec8))));
title('Fourier transform of reconstructed signal');
fprintf('Reconstruction from factor 8 down-sampled signal (no noise)\n');
fprintf('Reconstruction error: %e\n',norm(sbl-sblrec8));

%% 
sblrec8=B8*upsampl8(1:8:end);
figure(5);
subplot(2,1,1);plot(0:1023,sblrec8);
title('Reconstructed signal from factor 8 down-sampled signal');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec8))));
title('Fourier transform of reconstructed signal');
fprintf('Reconstruction from factor 8 down-sampled signal (no noise)\n');
fprintf('Reconstruction error: %e\n',norm(sbl-sblrec8));

err = sbl-sblrec8
mean_er = mean(err)
std_er = std(err)

%% 3.2  Sampling noise
sigma=0.1;
downsampl8n=sbl(1:8:end)+sigma*randn(128,1);
upsampl8n=zeros(1,1024);
upsampl8n(1:8:end)=downsampl8n;
figure(6);
subplot(2,1,1);plot(0:1023,upsampl8n);
title('The signal down-sampled with a factor 8 with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl8n))));
title('Fourier transform of down-sampled signal with noise');



sblrec8n=B8*downsampl8n;
figure(7);subplot(2,1,1);plot(0:1023,sblrec8n);
title('Reconstructed signal from factor 8 down-sampled signal with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec8n))));
title('Fourier transform of reconstructed signal with noise');
err=sbl-sblrec8n;
fprintf('Reconstruction from factor 8 down-sampled signal with noise\n');
fprintf('Reconstruction error: mean %f std %f\n',mean(err),std(err));

%Reconstruction from factor 8 down-sampled signal with noise
%Reconstruction error: mean 0.000658 std 0.097498
%Reconstruction from factor 8 down-sampled signal with noise
%Reconstruction error: mean 0.011691 std 0.101214
%Reconstruction from factor 8 down-sampled signal with noise
%Reconstruction error: mean -0.015419 std 0.094313
%Reconstruction from factor 8 down-sampled signal with noise
%Reconstruction error: mean -0.002645 std 0.091239

%% 3.3 Over-sampling, reconstruction with a basis
sigma=0.1;
downsampl4n=sbl(1:4:end)+sigma*randn(256,1);
upsampl4n=zeros(1,1024);
upsampl4n(1:4:end)=downsampl4n;
figure(8);
subplot(2,1,1);plot(0:1023,upsampl4n);
title('The signal down-sampled with a factor 4 with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl4n))));
title('Fourier transform of down-sampled signal with noise');



Rect4=(abs(u)<=pi/4);
rect4=ifft(fftshift(4*Rect4));
figure(9);
subplot(2,1,1);plot(-512:511,ifftshift(rect4));
title('Reconstructing filter rect4');
subplot(2,1,2);plot(u,abs(fftshift(fft(rect4))));

for ix=0:255,
    B4(:,ix+1)=circshift(rect4,4*ix);
end
figure(10);plot(0:1023,B4(:,[1 2 3 4])');
title('Reconstructing functions 3.3')

sblrec4b=B4*downsampl4n;
figure(11);subplot(2,1,1);plot(0:1023,sblrec4b);
title('Reconstructed signal from factor 4 down-sampled signal (basis)');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec4b))));
title('Fourier transform of reconstructed signal (basis)');
err=sbl-sblrec4b;
fprintf('Reconstruction from factor 4 down-sampled signal (basis)\n');
fprintf('Reconstruction error: mean %f std %f\n',mean(err),std(err));

% 3.4 Over-sampling, reconstruction with a frame

for ix=0:255,
F8(:,ix+1)=circshift(rect8,4*ix);
end
figure(13);plot(0:1023,F8(:,[1 2 3 4])');
title('Reconstructing functions 3.4')

%
sblrec4f=0.5*F8*downsampl4n;
figure(15);
subplot(2,1,1);plot(0:1023,sblrec4f);
title('Reconstructed signal from factor 4 down-sampled signal (frame)');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec4f))));
title('Fourier transform of reconstructed signal (frame)');
err=sbl-sblrec4f;
fprintf('Reconstruction from factor 4 down-sampled signal (frame)\n');
fprintf('Reconstruction error: mean %f std %f\n',mean(err),std(err));

%%
[U SB V] = svd(B4);
eigB =diag(SB)

[U SF V] = svd(F8);
eigF = diag(SF)
%% 4 DWT of 1D signals

l=256;
u=((-l/2):(l/2-1))*2*pi/l;
s0=rand(1,l);
S0=fftshift(fft(ifftshift(s0)));
S=S0.*(abs(u)<pi/4);
s=real(ifftshift(ifft(fftshift(S))));
%%
figure(1);plot(0:255,s);
title('Input signal');

[h0 g0 h1 g1]=wfilters('db3')

flength=length(h0);
u1=((-flength/2):(flength/2-1))*2*pi/flength;
figure(2);
subplot(4,1,1);plot(u1,abs(fftshift(fft(h0))));title('FT of h0');
subplot(4,1,2);plot(u1,abs(fftshift(fft(g0))));title('FT of g0');
subplot(4,1,3);plot(u1,abs(fftshift(fft(h1))));title('FT of h1');
subplot(4,1,4);plot(u1,abs(fftshift(fft(g1))));title('FT of g1');

%% TESTING 
figure(3)
subplot(4,1,1);plot(h0);title('Signal domain of h0');
subplot(4,1,2);plot(h1);title('Signal domain of h1');
subplot(4,1,3);plot(g0);title('Signal domain of g0');
subplot(4,1,4);plot(g1);title('Signal domain of g1');
%%
figure(3);
subplot(4,1,1);plot(u1,abs(fftshift(fft(h0))));title('FT of h0');
subplot(4,1,2);plot(u1,abs(fftshift(fft(g0))));title('FT of g0');
subplot(4,1,3);plot(u1,abs(fftshift(fft(h1))));title('FT of h1');
subplot(4,1,4);plot(u1,abs(fftshift(fft(g1))));title('FT of g1');

%%
dwtmode('per');  %Set periodic mode of filtering operations
[a d]=dwt(s,h0,g0);
figure(3);
subplot(2,1,1);plot(a);title('a');
subplot(2,1,2);plot(d);title('d');

srec1=idwt(a,d,h1,g1);
figure(4);plot(0:255,srec1);
title('Reconstructed signal');
%%
[a1 d1]=dwt(s,h0,g0);
[a2 d2]=dwt(a1,h0,g0);
figure(5);
subplot(3,1,1);plot(a2);title('a2');
subplot(3,1,2);plot(d2);title('d2');
subplot(3,1,3);plot(d1);title('d1');

%%
a1rec=idwt(a2,d2,h1,g1);
srec2=idwt(a1rec,d1,h1,g1);
figure(6);plot(srec2);title('Reconstructed signal');

%%
% Code snippet (A)
N=2;ad=s;p=length(s);figure(7);
p
for cnt=1:N,
    [a d]=dwt(ad(1:p),h0,g0);
    a
    d
    ad(1:p)=[a d];
    ad
    subplot(N+1,1,N+2-cnt);
    plot(d);title(sprintf('details level %d',cnt));
    p=p/2;
end
subplot(N+1,1,1);plot(a);
title(sprintf('approximation level %d',cnt));
figure(8);plot(ad);
title('Concatenated approximation and details');
%% 4.3 Simple signal compression

q = [16 3;16 2;16 1;16 0.1];  %Change this if needed
[ad bps] = quantisead(ad,q);  %Replace the channels with quantised values

%%
% Code snippet (B)
for cnt=1:N,
    ad(1:(2*p))=idwt(ad(1:p),ad((p+1):(2*p)),h1,g1);
    p=2*p;
end
figure(9);plot(ad);title('Reconstructed signal');

err=s-ad;

fprintf('Reconstruction error: mean %f std %f\n',mean(err),std(err));
