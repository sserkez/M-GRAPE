%split-step
clear all
%close all

calculation='CPU';

existsurface=1;


        N_a=fliplr([256,512,1024,2048]);
        dz_a=([0.25e-9,0.5e-9,1e-9,2e-9,3e-9,4e-9]);
        
N_a=2048*6;
dz_a=2e-9;

for m=1:numel(N_a)
    for j=1:numel(dz_a)

        N=N_a(m);
        dz=dz_a(j);
        
%N=1024;
L=4e-5;
xlamds=1.54e-9;
K0=2*pi/xlamds;
w0=20e-8;
%dz=10e-9;
Z=5e-5;
Theta=deg2rad(1);
x_offset=0;
surf_height=1e-7;
surf_height=L/2-5e-7;

dx=L/N;
k=round(Z/dz);
Z_sc=(0:k-1)*Z;

sc=-(N-1)/2:1:(N-1)/2;
x=sc*dx;
%  fx=linspace(-1/(2*dx)+1/L/2,1/(2*dx)-1/L/2,N);
fx=linspace(-2*pi/(2*dx)+1/L/2,2*pi/(2*dx)-1/L/2,N);


X=(exp(-(x+x_offset).^2./w0^2)).*exp(1i*Theta*K0*x);
disp(' ');
disp(['Transverse step ',num2str(dx),' m']);
disp(['Longitudinal step ',num2str(dz),' m']);
disp(['Transverse phase step ',num2str(dx*Theta/xlamds),' rad']);
% figure(1)
% plot(x,abs(X))
% 
%         H=exp(-1j*pi*xlamds*Z*fx.^2);
% 
%              u1=fftshift(fft(ifftshift(X)));
%              %fieldplot(301,u1,1,'field before',1)
%              %fieldplot(302,H,1,'propagator',1)
%              u1=H.*u1;                    
%              %fieldplot(303,u1,1,'field after',1) %!original
%              X2=fftshift(ifft(ifftshift(u1)));
% 
% figure(2)
% hold all
% plot(x,abs(X2))

%tic

%if calculation=='CPU'
%
epsilon=((1-0.0015)+1i*0.00024)-1;
%epsilon=-i;
%epsilon=0;
depsilon=(epsilon-1);
%depsilon=1;
%%Surface=zeros(N,k);
Ns=round(N/L*surf_height);
%Surface(N-Ns:N,:)=existsurface;

Surface=zeros(1,N);
Surface(N-Ns:N)=existsurface;

% figure
% imagesc(Surface)
%break
%Surface(:,1:50)=1;
%Epsilon=depsilon.*Surface;

    
    H=ifftshift(exp(-1j/(2*K0)*dz*fx.^2));
    %H=ifftshift(exp(0));
    M=exp(1i*K0*depsilon.*Surface*dz/2);
    %M=M(:,1);
    X_prop=zeros(N,k,'single');
    clear Surface
    
    u0=X;
    
    for i=1:k


%                  u0=fftshift(fft(ifftshift(u0)));
%                  u0=H.*u0;                    
%                  u0=fftshift(ifft(ifftshift(u0)));
                 
                 u0=ifft(H.*fft(u0));
                 
                 %u0=M(:,i)'.*u0;
                 
                 
                 u0=M.*u0;


%                  u0=fft(u0);
%                  %fieldplot(301,u1,1,'field before',1)
%                  %fieldplot(302,H,1,'propagator',1)
%                  u0=ifftshift(H).*u0;                    
%                  %fieldplot(303,u1,1,'field after',1) %!original
%                  u0=ifft(u0);

                 X_prop(:,i)=abs(u0).^2;

    end
    X3=u0;
    
    I=sum(X_prop,1);
    
    
    figure(5)
    imagesc(Z_sc,x,X_prop);
%    clear all
    figure(6)
    plot(Z_sc,I./max(I));
    
    clear X_prop
    
    disp(['Energy loss ',num2str((1-I(end)/I(1))*100),' %']);
    
    Eloss(m,j)=(1-I(end)/I(1))*100;
    
    end
end

% imagesc(N_a,dz_a,Eloss)
% figure
% imagesc(Eloss)

