%tilt simulation
clear all

%% set parameters

energy=1000;
xlamds=1239.8/energy*1e-9;
prop_leng=0;

M=251;
N=251;
sigm_xy=3e-5;
sigm_z=0.5e-6;

leng_xy=20e-4;
leng_z=3e-5;

z=-(N-1)/2:1:(N-1)/2;
xy=-(M-1)/2:1:(M-1)/2;
dz=leng_z/N;
dxy=leng_xy/M;


k=2*pi/xlamds;
dk=2*pi/(leng_z);
K=k+dk*z;

P=1000*exp(-(dz*z).^2./sigm_z^2);

% !!!
gr_sigma=1/1123000;
Teta_g_i=0.017;
Teta_g_d=acos(cos(Teta_g_i)-xlamds/gr_sigma);
% !!!

%% create the field

XX=single(zeros(M,M,N));
for i=1:N
    XX(:,:,i)=fieldgaussian(M,leng_xy,1e-4,1e-4,0,0,xlamds,P(i));
end
%% show field

%[H{1}]=fieldplot(1,mean(XX,3),leng_xy,'horizontal projection',1);

figure(2)
subplot(2,1,1)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX).^2,2),M,N));
xlabel('z length [m]');
ylabel('y length [m]');
%axis('square')

subplot(2,1,2)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX).^2,1),M,N));
xlabel('z length [m]');
ylabel('x length [m]');


%% modify field

dTheta=dk*z*xlamds/(k*Teta_g_d*gr_sigma);

XXf=fftshift(fft(XX,[],3),3);
for i=1:N
    XXf(:,:,i)=abs(XXf(:,:,i)).*exp(1i*(angle(XXf(:,:,i))+dTheta(i)*k*meshgrid(xy)'*dxy));
    XXf(:,:,i)=prop_TF(XXf(:,:,i), leng_xy, xlamds,prop_leng);
end
spectrum=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);

% figure(5)
% plot(2*pi./K,spectrum)

% figure(2)
% imagesc(K,xy*dxy,reshape(mean(abs(XXf).^2,2),M,N));
% xlabel('z length [m]');
% ylabel('y length [m]');
% 
% figure(3)
% imagesc(K+k,xy*dxy,reshape(mean(abs(XXf).^2,1),M,N));
% xlabel('z length [m]');
% ylabel('x length [m]');

%% show filtered field

XX1=ifft(ifftshift(XXf,3),[],3);

figure(3)
subplot(2,1,1)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX1).^2,2),M,N));
xlabel('z length [m]');
ylabel('y length [m]');

subplot(2,1,2)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX1).^2,1),M,N));
xlabel('z length [m]');
ylabel('x length [m]');

%[H{4}]=fieldplot(4,XX1(:,:,151),leng_xy,'filtered field',1);