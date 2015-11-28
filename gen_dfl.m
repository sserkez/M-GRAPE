xlamds=1.334e-9;%wavelength

M=201; %transverse mesh points MxM
N=1631; %longitudinal mesh points
zsep=16; %distance between slices in xlamds
mesh_width=8e-4; %transverse mesh size aka "leng"
mesh_length=zsep*xlamds*N; %longitudinal mesh size

pulse_length=7e-6; %pulse length in S
pulse_position=15e-6; %pulse position in S
pulse_energy=4e-9; %Joules
pulse_waist_width=15e-6;%waist rms width
pulse_waist_position=0;% waist position in Z



P0=10000; %dummy power
X=single(fieldgaussian(M,mesh_width,pulse_waist_width,pulse_waist_width,pulse_waist_position,pulse_waist_position,xlamds,P0)); %transverse intensity modulation + focusing
fieldplot(101,X,mesh_width,'flat field',1);

L=linspace(0,mesh_length,N);
X=repmat(X,[1,1,N]);
l_mod=single(exp(-(L-pulse_position).^2./pulse_length.^2));
l_mod=reshape(l_mod,[1,1,N]); %longitudinal intensity modulation
%break
l_mod=repmat(l_mod,[M,M,1]);

Xi=X.*l_mod;
clear X l_mod

E0=sum(sum(sum(abs(Xi).^2)))/N*mesh_length/3e8;
Xi=Xi.*sqrt(pulse_energy/E0);
% Xi(:,:,1:500)=0;
% Xi(:,:,1000:end)=0;
fieldplot3d(102,Xi,mesh_width,1,L,'bulk field',1);

Xf=fftshift(fft(Xi,[],3),3)./sqrt(N);

%calculate wavelength scale
sc=-(N-1)/2:1:(N-1)/2;
k=2*pi/xlamds;
dk=2*pi/(N*xlamds*zsep);
K=k+dk*sc;
Xlamds=2*pi./K;
%

fieldplot3d(103,Xf,mesh_width,0,Xlamds,'spectrum',1);
%Xi=ifft(ifftshift(Xf,3),[],3).*sqrt(N);

%fieldplot3d(104,Xi,mesh_width,1,L,'bulk field after cut',1);
%%
fieldexport(Xi,'C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_perfect_seed');