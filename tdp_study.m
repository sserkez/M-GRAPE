clear all

nm_p='1000_u1_tdp.out';
inread;

monoline=1.242e-9;

%M=151;
%Nslice=700;
%xlamds=1.770e-09 ;
%curlen=4.784e-06;
Npad=3;
%zsep=16;
%XX=fieldimport_tdp('700_tdp.out.dfl',M,Nslice,1);
[XX,Nslice]=fieldimport_all('1000_u1_tdp.out.dfl',M,1);
%XX(:,:,300:end)=0;
%XX=flipdim(XX,3);
XX=cat(3,zeros(M,M,Nslice*(Npad-1)),XX);
Nslice1=Nslice*Npad;

slice_scale=linspace(0,Nslice1*xlamds*zsep,Nslice1);
slice_power=reshape(sum(sum(abs(XX).^2)),1,[]);

figure(53)
plot(slice_scale, slice_power);
Int2_1=trapz(slice_scale, slice_power);


%XX(M,M,Nslice*(Npad))=0;%[XX zeros(1,Nslice*(Npad-1))];

sc=linspace(-Nslice1/2,Nslice1/2,Nslice1);

k0=2*pi/xlamds;
dk=1/(Nslice1*xlamds*zsep);
k=k0+2*pi*dk*sc;

XXf=fftshift(fft(XX,[],3),3);
clear XX
%spectrum=reshape(abs(XXf(ceil(M/2),ceil(M/2),:)).^2,1,[]);
spectrum=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);

figure(54);
plot(2*pi./k,spectrum);
Int1_1=trapz(k,spectrum);
%trapz(-2*pi./k,spectrum)
% XX1=ifft(ifftshift(XXf,3),[],3);
% slice_power1=reshape(sum(sum(abs(XX1).^2)),1,[]);
% slice_scale1=linspace(0,Nslice1*xlamds*zsep,Nslice1);
% figure(55)
% plot(slice_scale1, slice_power1);

%%
load('E:\WORK\New_matlab\Presentation\Resolution\1000.mat') % to get Power and resolshift
%Power=sqrt(Power);
Power=Power-min(Power)*0.9999;
% [k_n,filt_n,phs_n]=KKphase(9,xlamds,resolshift,Power/max(Power).*0.05);
lam = xlamds + resolshift'*xlamds;
k0= 2*pi./(lam(:)-xlamds+monoline);
Power=interp1([k(1) k0' k(end)],[Power(1) Power Power(end)],k);
filt=Power/max(Power)*0.05;
filt0=-log(filt)/2;
phs=-1*kkimbook2(k,filt0,0);
%phs=KKphase1(k,filt);

% filt=interp1([k(1) k_n k(end)],[min(filt_n) filt_n min(filt_n)],k);
% phs=interp1([k(1) k_n k(end)],[0 phs_n 0],k);
filt_cmplx=sqrt(filt).*exp(1i*phs);
%filt_cmplx=fliplr(filt_cmplx);

figure(54);
%plot(2*pi./k,spectrum/max(spectrum));
plot(2*pi./k,spectrum);
hold all
% plot(2*pi./k,filt/max(filt));
% plot(2*pi./k,phs/max(phs));
plot(2*pi./k,abs(filt_cmplx.^2)/max(abs(filt_cmplx.^2)));
plot(2*pi./k,angle(filt_cmplx)/max(angle(filt_cmplx)));
hold off

%XXff=XXf;

for i=1:size(XXf,1)
    for j=1:size(XXf,2)
        XXf(i,j,:)=reshape(XXf(i,j,:),1,[]).*filt_cmplx;
    end
end

%clear XXf

figure(54);
    spectrum=reshape(sum(sum(abs(XXf).^2)),1,[]);
    hold all
%    plot(2*pi./k,spectrum/max(spectrum));
    plot(2*pi./k,spectrum);
    Int1_2=trapz(k,spectrum);
    hold off
    Sigm_omega=findFWHM(3e8*k,spectrum);
    
figure(56);
     XXf=ifft(ifftshift(XXf,3),[],3);
     XX1=XXf;
     clear XXf
%    XX1=ifftshift(ifft(XXf,[],3),3);
    slice_power1=reshape(sum(sum(abs(XX1).^2)),1,[]);
    slice_scale1=linspace(0,Nslice1*xlamds*zsep,Nslice1);
    plot(slice_scale1, slice_power1);
    Int2_2=trapz(slice_scale1, slice_power1);
    Sigm_time=findFWHM(slice_scale1, slice_power1)*1e6*1e-15*3;
    
    disp(Sigm_omega*Sigm_time);
    %clear XX1 XXff
    %sum(slice_power)

% for i=1:numel(k)
%     [H{1}]=fieldplot(1,XX1(:,:,i),leng,'filtered field',1);
% end
%%
n=find(slice_scale1>4.7e-5,1, 'first');
fieldexport(XX1(:,:,n:end),'1000_tdp_filtered.dfl');
%%
%[XX1,Nslice]=fieldimport_all('1000_tdp_filtered.dfl',151,1);
%[H{1}]=fieldplot(1,mean(XX1,3),4E-4,'filtered field',1);