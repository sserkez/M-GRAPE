clear all

nm_f='tdp\1000_u1_filtered_u8_DA.dfl';
nm_p='tdp\1000_u2_tdp_u8_DA.out';
nm_f=[nm_p,'.dfl'];
Npad=4;
inread;


[XX,Nslice]=fieldimport_all(nm_f,M,1);
[H{1}]=fieldplot3d(1,XX,leng,1,linspace(0,xlamds*zsep*Nslice,Nslice),'initial time domain',1);
XX=cat(3,zeros(M,M,Nslice*(Npad-1)),XX);
Nslice1=Nslice*Npad;

slice_scale=linspace(0,Nslice1*xlamds*zsep,Nslice1);
slice_power=reshape(sum(sum(abs(XX).^2)),1,[]);

% figure(53)
% plot(slice_scale, slice_power,'linewidth',2);
% xlabel('s[m]');
% ylabel('P[W]');

%sc=linspace(-Nslice1/2,Nslice1/2,Nslice1);%original
sc=-(Nslice1-1)/2:1:(Nslice1-1)/2;

k0=2*pi/xlamds;

dk=2*pi/(Nslice1*xlamds*zsep);
k=k0+dk*sc;

XXf=fftshift(fft(XX,[],3),3);
clear XX XXplot
%spectrum=reshape(abs(XXf(ceil(M/2),ceil(M/2),:)).^2,1,[]);
spectrum1=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);
clear XXf
%%
nm_f='tdp\1000_u1_filtered_u7_DA.dfl';
nm_p='tdp\1000_u2_tdp_u8_KK.out';
nm_f=[nm_p,'.dfl'];
inread;


[XX,Nslice]=fieldimport_all(nm_f,M,1);
[H{2}]=fieldplot3d(2,XX,leng,1,linspace(0,xlamds*zsep*Nslice,Nslice),'initial time domain',1);
XX=cat(3,zeros(M,M,Nslice*(Npad-1)),XX);
Nslice1=Nslice*Npad;

slice_scale=linspace(0,Nslice1*xlamds*zsep,Nslice1);
slice_power=reshape(sum(sum(abs(XX).^2)),1,[]);

% figure(53)
% plot(slice_scale, slice_power,'linewidth',2);
% xlabel('s[m]');
% ylabel('P[W]');

%sc=linspace(-Nslice1/2,Nslice1/2,Nslice1);%original
sc=-(Nslice1-1)/2:1:(Nslice1-1)/2;

k0=2*pi/xlamds;

dk=2*pi/(Nslice1*xlamds*zsep);
k=k0+dk*sc;

XXf=fftshift(fft(XX,[],3),3);
clear XX XXplot
%spectrum=reshape(abs(XXf(ceil(M/2),ceil(M/2),:)).^2,1,[]);
spectrum2=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);
clear XXf
%%

figure(3354);
plot(2*pi./k*1e9,spectrum1,'linewidth',2);
hold all
plot(2*pi./k*1e9,spectrum2,'linewidth',2,'color','r','linestyle','--');
hold off
xlabel('\lambda[nm]');
ylabel('P(\lambda)[a.u.]');
xlim([min(2*pi./k*1e9) max(2*pi./k*1e9)])

    legend('direct approach','phenomenological approach','Location','northeast')
    xlim([1.2386 1.2395]);
    %ylim([0 18e14]);
    %ylim([0 2.6e11]);
    %xlim([1.764 1.766]);
    ylabel('P(\lambda) [a.u.]');
    xlabel('\lambda [nm]');