clear all

% nm_p='1000_tdp.out';
% nm_f='1000_tdp.out.dfl';

nm_p='c:\-D-\Work\LCLS\New_matlab\tdp\500_u1_tdp.out';
nm_f=[nm_p,'.dfl'];
load('C:\-D-\Work\LCLS\New_matlab\Steadystate\500_8_new.mat');
%ICF=0.8*0.058; % 500:0.4/0.8 700:0.6/1.3 1000:0.8/0.8    0.05 0.057 0.058
ICF=0.4*0.05; % 500:0.79/0.31 700:2/0.51 1000:0.61/0.52  (U8/U7)    0.05 0.057 0.058
% xlamds0=1.2432e-9;
% xlamds0=1.7652e-9;
 xlamds0=2.485e-9; %2.468
 xlamds0=2.489e-9;

linewigdth=1.5;
%lamdlim=[1.240 1.244];
lamdlim=[2.45 2.49];
%lamdlim=[1.758 1.774];
lenglim=[0.6e-4 1.3e-4];
%lenglim=[0.8e-4 1.4e-4];
figname='500_7';

inread;

%M=151;
%Nslice=700;
%xlamds=1.770e-09 ;
%curlen=4.784e-06;
Npad=2; %3
leng_u2=10e-4;%8e-4
M_u2=201;
%zsep=16;
%XX=fieldimport_tdp('700_tdp.out.dfl',M,Nslice,1);
[XX,Nslice]=fieldimport_all(nm_f,M,1);


%XX=flipdim(XX,3);
%XX(M,M,Nslice*(Npad))=0;%[XX zeros(1,Nslice*(Npad-1))];
slice_scale0=linspace(0,Nslice*xlamds*zsep,Nslice);

XX=cat(3,zeros(M,M,Nslice*(Npad-1)),XX);
Nslice1=Nslice*Npad;

slice_scale=linspace(0,Nslice1*xlamds*zsep,Nslice1);
slice_power=reshape(sum(sum(abs(XX).^2)),1,[]);

figure(53)
plot(slice_scale, slice_power,'linewidth',linewigdth);
xlabel('s[m]');
ylabel('P[W]');



%XX(M,M,Nslice*(Npad))=0;%[XX zeros(1,Nslice*(Npad-1))];

%         z=-(N-1)/2:1:(N-1)/2;
%         
%         k=2*pi/xlamds;
%         dk=2*pi/(leng_z);
%         K=k+dk*z;
%         Xlamds=2*pi./K;

sc=linspace(-Nslice1/2,Nslice1/2,Nslice1);%original
sc=-(Nslice1-1)/2:1:(Nslice1-1)/2;

k0=2*pi/xlamds;

dk=2*pi/(Nslice1*xlamds*zsep);
k=k0+dk*sc;

XXf=fftshift(fft(XX,[],3),3);
XXf=prop_TF(XXf,leng,xlamds,0); % one undulator backpropagate this case only!!!

XXplot=ifft(ifftshift(XXf,3),[],3);
XXplot=XXplot(:,:,Nslice1-Nslice+1:end);
[H{554}]=fieldplot3d(554,XXplot,leng,1,slice_scale0,'original field at end of undulator',1);

clear XX XXplot
%spectrum=reshape(abs(XXf(ceil(M/2),ceil(M/2),:)).^2,1,[]);
spectrum=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);

figure(54);
plot(2*pi./k*1e9,spectrum,'linewidth',linewigdth);
xlabel('\lambda[nm]');
ylabel('P(\lambda)[a.u.]');
xlim([min(2*pi./k*1e9) max(2*pi./k*1e9)])

% XX1=ifft(ifftshift(XXf,3),[],3);
% slice_power1=reshape(sum(sum(abs(XX1).^2)),1,[]);
% slice_scale1=linspace(0,Nslice1*xlamds*zsep,Nslice1);
% figure(55)
% plot(slice_scale1, slice_power1);

%%

%load('C:\-D-\WORK\LCLS\New_matlab\Resolution\1000.mat')

Power=inppowerout;
%Power=fliplr(inppowerout); %experimental!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
resolshift=idlpl_tilt;
Power=Power-min(Power)*0.999;
% [k_n,filt_n,phs_n]=KKphase(9,xlamds,resolshift,Power/max(Power).*0.05);

lam = xlamds0 + resolshift'*xlamds;
k0= 2*pi./lam(:);
Power(1)=Power(1)/100;
Power(end)=Power(end)/100;
Power=interp1([k(end) k0' k(1)],[Power(1) Power Power(end)],k);
filt=Power/max(Power)*sqrt(ICF);
phs=KKphase1(k,filt);

% filt=interp1([k(1) k_n k(end)],[min(filt_n) filt_n min(filt_n)],k);
% phs=interp1([k(1) k_n k(end)],[0 phs_n 0],k);
filt_cmplx=filt.*exp(1i*phs);
%filt_cmplx=fliplr(filt_cmplx);

figure(55);
plot(2*pi./k,spectrum/max(spectrum),'linewidth',1.5);
hold all
plot(2*pi./k,filt/max(filt),'linewidth',2);
plot(2*pi./k,phs/max(phs),'linewidth',2,'color','k');
hold off

%XXff=XXf;

for i=1:size(XXf,1)
    for j=1:size(XXf,2)
        XXf(i,j,:)=reshape(XXf(i,j,:),1,[]).*filt_cmplx;
    end
end

%clear XXf

spectrum2=reshape(sum(sum(abs(XXf).^2)),1,[]);
figure(55);
    hold all
    plot(2*pi./k,spectrum2/max(spectrum2),'linewidth',2,'color','r');
    hold off
    Sigm_omega=findFWHM(3e8*k,spectrum2);

%    XXf=prop_TF(XXf, leng, xlamds,5); %remove!!!
    
figure(58);
plot(2*pi./k*1e9,spectrum2,'linewidth',linewigdth);
xlabel('\lambda[nm]');
ylabel('P(\lambda)[a.u.]');
xlim([min(2*pi./k*1e9) max(2*pi./k*1e9)]);


figure(56);
    XX1=ifft(ifftshift(XXf,3),[],3);
    clear XXf
%    XX1=ifftshift(ifft(XXf,[],3),3);
    slice_power1=reshape(sum(sum(abs(XX1).^2)),1,[]);
    slice_scale1=linspace(0,Nslice1*xlamds*zsep,Nslice1);
    plot(slice_scale1, slice_power1,'linewidth',linewigdth);
    xlabel('s[m]');
    ylabel('P[W]');

    Sigm_time=findFWHM(slice_scale1, slice_power1)*1e6*1e-15*3;
    
    disp(Sigm_omega*Sigm_time);
    %clear XX1 XXff


% for i=1:numel(k)
%     [H{1}]=fieldplot(1,XX1(:,:,i),leng,'filtered field',1);
% end


Nslice2=Nslice;

%zsep=25;
fieldshiftsize=11e-6; %1000
%fieldshiftsize=6e-6; %
%fieldshiftsize=15e-6;  %500
%fieldshiftsize=0; %Remove
sliceshift=round(fieldshiftsize/xlamds/zsep);
%XX2=circshift(XX1,[0 0 sliceshift]);
XX2=XX1(:,:,Nslice1-Nslice-sliceshift+1:Nslice1-Nslice-sliceshift+Nslice2-1+1);
%slice_scale1=slice_scale(Nslice1-Nslice-sliceshift+1:Nslice1-Nslice-sliceshift+Nslice2-1+1);
slice_scale1=slice_scale(1:Nslice2);
%clear XX1

    XX3=zeros(M_u2,M_u2,size(XX2,3),'single');
    for i=1:size(XX2,3)
    XX3(:,:,i)=fieldinterpolate(XX2(:,:,i),leng,0,M_u2,leng_u2,'spline');
    end
%    clear XX2
%XX3=prop_TF(XX3, leng_u2, xlamds,10); %remove? 
[H{555}]=fieldplot3d(555,XX3,leng_u2,1,slice_scale1,'field at circshift',1);

break
%%

set(figure(53), 'Position', [100, 100, 450, 300]);
set(figure(54), 'Position', [100, 100, 450, 300]);
set(figure(56), 'Position', [100, 100, 450, 300]);
set(figure(58), 'Position', [100, 100, 450, 300]);
set(get(figure(53),'child'),'XLim',lenglim); 
set(get(figure(54),'child'),'XLim',lamdlim); 
set(get(figure(56),'child'),'XLim',lenglim); 
set(get(figure(58),'child'),'XLim',lamdlim); 
%saveas(figure(53),['fig\KK_',figmane,'_p_b'],'eps');

break

hgexport(figure(53), ['fig\KK_',figname,'_p_b'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(figure(56), ['fig\KK_',figname,'_p_a'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(figure(54), ['fig\KK_',figname,'_s_b'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(figure(58), ['fig\KK_',figname,'_s_a'], hgexport('factorystyle'), 'Format', 'eps');


%%
figure(58)
xlim([xlamds0*1e9-5e-4 xlamds0*1e9+5e-4]);
%xlim([1.2387 1.2393]);
set(figure(58), 'Position', [100, 100, 250, 200]);
set(get(figure(58),'child'),'XTick',lamdlim(1):0.0004:lamdlim(2)); 
break
hgexport(figure(58), ['fig\KK_',figname,'_s_a_enl'], hgexport('factorystyle'), 'Format', 'eps');

%%
fieldexport(XX3,'C:\-D-\Work\LCLS\tmp\3\damage\500_u1_1.out_fk.dfl');