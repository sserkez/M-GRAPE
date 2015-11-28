%tilt simulation
clear all


%% set parameters

realfield=1;
prop_leng=0;

if realfield
    
    nm_f='1000_u1_tdp.out.dfl';
    nm_p='1000_u1_tdp.out';
    inread;
    [XX,N]=fieldimport_all(nm_f,M,1); %import field
    leng_xy=leng;
    dz=xlamds*zsep;
    dxy=dx;
    leng_z=dz*N;
    
    z=-(N-1)/2:1:(N-1)/2;
    xy=-(M-1)/2:1:(M-1)/2;
    
else
    
    energy=1000;
    xlamds=1239.8/energy*1e-9;
    
    M=251;
    N=351;
    sigm_xy=1e-5;
    sigm_z=0.5e-6;

    leng_xy=8e-4;
    leng_z=2e-5;

    z=-(N-1)/2:1:(N-1)/2;
    xy=-(M-1)/2:1:(M-1)/2;
    dz=leng_z/N;
    dxy=leng_xy/M;

    P=1000*exp(-(dz*z).^2./sigm_z^2);

    XX=single(zeros(M,M,N));

    for i=1:N
        XX(:,:,i)=fieldgaussian(M,leng_xy,1e-4,1e-4,0,0,xlamds,P(i)); %generate field
    end
    
%     XX(:,:,(N-1)/2:end)=0;
%     XX(1:(M-1)/2,:,:)=0;
    
end

k=2*pi/xlamds;
dk=2*pi/(leng_z);
K=k+dk*z;

% !!!
gr_sigma=1/1123000;
Teta_g_i=0.017;
Teta_g_d=acos(cos(Teta_g_i)-xlamds/gr_sigma);
% !!!

%% show field

%[H{1}]=fieldplot(1,mean(XX,3),leng_xy,'horizontal projection',1);

% XX=fftshift(fft(XX,[],3),3)./sqrt(N);
% XX=prop_TF(XX,leng,xlamds,-4.5);   
% Xt=ifft(ifftshift(XX,3),[],3).*sqrt(N);

figure(2)
subplot(2,1,1)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX).^2,2),M,N));
%xlim([-1.4e-5 -0.8e-5]);
xlabel('z length [m]');
ylabel('y length [m]');
title('x projection');
%axis('square')

subplot(2,1,2)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX).^2,1),M,N));
%xlim([-1.4e-5 -0.8e-5]);
xlabel('z length [m]');
ylabel('x length [m]');
title('y projection');


%% modify field

dTheta=dk*z*xlamds/(k*Teta_g_d*gr_sigma);

XX=fftshift(fft(XX,[],3),3);
XXf=XX;
clear XX

for i=1:N
    XXf(:,:,i)=abs(XXf(:,:,i)).*exp(1i*(angle(XXf(:,:,i))+dTheta(i)*k*meshgrid(xy)'*dxy));
    if prop_leng~=0
        XXf(:,:,i)=prop_TF(XXf(:,:,i), leng_xy, xlamds,prop_leng);
    end
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

XXf=ifft(ifftshift(XXf,3),[],3);
XX1=XXf;
clear XXf

figure(3)

subplot(2,1,1)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX1).^2,2),M,N));
%xlim([-1.4e-5 -0.8e-5]);
xlabel('z length [m]');
ylabel('y length [m]');
title('x projection');
%title(['z=',num2str(prop_leng),'m']);

subplot(2,1,2)
imagesc(z*dz,xy*dxy,reshape(mean(abs(XX1).^2,1),M,N));
%xlim([-1.4e-5 -0.8e-5]);
xlabel('z length [m]');
ylabel('x length [m]');
title('y projection');

%[H{4}]=fieldplot(4,XX1(:,:,151),leng_xy,'filtered field',1);