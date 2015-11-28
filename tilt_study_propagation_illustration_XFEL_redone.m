%Mono for new XFEL fixed M2 position

%clear all;
close all;
fclose all;

nm_p='c:\-D-\Work\SASE3_SXRSS\1000_u1_tdp.out';
nm_f=[nm_p,'.dfl'];

scaledenergy=5e-6/7;

% resolution without slit calculation
t_start=tic;
showpictures=1;
realfield=1;
onetoone=0;
throughslit=0;
roughnessincluded=0;
aberrationsincluded=1;
dispersionon=1;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);

%% XFEL 700 run simulation parameters

interpM=1;%4 2 8 3 6
interpL=2;%3 2 3 4 2 %4,3,2 (8,5,1.3) for 930eV LCLS run simulation
Npad=2; %2 2 1 2 1

transmission=0.06; %recalculate

leng_u2=8e-4;%10e-4
M_u2=201;%201 %930eV %for LCLS run simulation

dfl_shift=-1.09;
u2_prop=0; %propagation within U2 (remove, for display only!!!)
fieldshiftsize=0e-6;

export1to1=1;

xlamds0=1.781e-9; %700eV 
dxlamds0=0.002e-9; %700eV

xlamds0=1.246e-9; %1000eV 
dxlamds0=0.002e-9; %1000eV


%%
d=outread(nm_p,1); 
leng=d.inp.leng;
M=d.inp.ncar;
xlamds=d.inp.xlamds;
zsep=d.inp.zsep;
dx=leng/M;
if (xlamds-xlamds0)/xlamds>0.01
    error ('wrong central wavelength defined')
end

P0=10000; %input power for model gaussian
%% Optical system parameters LCLS

opt_parms_1;


g_length=0.025;         %grating optical length [m]
m2_length=0.025;
energy=1239.8411/xlamds*1e-9;
gafoc_XFEL;

%resolshift=-1.5e-4:2e-5:1.5e-4; %scanned resolution shift via tilt on the grating (=dlambda/lambda)
resolshift=-0e-4;
slitwidth=5e-6;
iterations=size(resolshift,2);
Power=linspace(0,0,iterations); %input coupling factor

    %dTheta_g_d=xlamds*resolshift(iteration)/Theta_g_d/gr_sigma;
    t_iter=tic;
    %%
    leng_i=leng*interpL;

    if realfield
        
        [X,N]=fieldimport_all(nm_f,M,1); %import field
%X=flipdim(X,2);


%         X=cat(3,zeros(M,M,round(N*(Npad-1))),X);
%         N2=round(N*Npad);
        X=cat(3,zeros(M,M,N*(Npad-1)),X);
        N2=N*Npad;
        dz=xlamds*zsep;
        dxy=dx;
        leng_z=dz*N2;

        %scale the power!  notice! notice! notice! notice! notice! notice! %notice!
        E0=sum(sum(sum(abs(X).^2)))/N2*leng_z/3e8;
        X=X.*sqrt(scaledenergy/E0);
        %

        z=-(N2-1)/2:1:(N2-1)/2;
        xy=-(M-1)/2:1:(M-1)/2;
        Z=(z-min(z))*dz;
        
        k=2*pi/xlamds;
        dk=2*pi/(leng_z);
        K=k+dk*z;
        Xlamds=2*pi./K;
        dxlamds=Xlamds-xlamds;
        
        %try
sc=-(N2-1)/2:1:(N2-1)/2;

k=2*pi/xlamds;
dk=2*pi/(N2*xlamds*zsep);
K=k+dk*sc;
Xlamds=2*pi./K;
        %try
        

        %npfiltshift=0;
        if dxlamds0==0
            npfilt1=1;
            npfilt2=N2;
        else
            npfilt1=find(Xlamds>(xlamds0+dxlamds0/2),1,'last');%number preliminary filter
            npfilt2=find(Xlamds<(xlamds0-dxlamds0/2),1,'first');%number preliminary filter
        end
        %npfilt_mean=mean(npfilt1,npfilt2);
        Nf=npfilt2-npfilt1+1;
        Xlamdsf=Xlamds(npfilt1:npfilt2);
        Zf=Z(npfilt1:npfilt2);
        % end of preliminary filtration parameters

        Theta_g_d_0=acos(cos(Theta_g_i)-xlamds0*D0*1e3);
        
        if dxlamds0==0
        Nfshift=(Xlamdsf(round(size(Xlamdsf,2)/2))-xlamds0)/(Xlamdsf(1)-(Xlamdsf(2)));
        z0=linspace(-Nf/2-Nfshift,Nf/2-Nfshift,Nf);
        else
         z0=linspace(-Nf/2,Nf/2,Nf);   
        end
        
        if dispersionon
            dTheta_g_d=dk*z0*xlamds0/(k*Theta_g_d_0*gr_sigma);
        else
            dTheta_g_d=0*z0;
        end
        
        %dTheta_g_d=dTheta_g_d(npfilt1:npfilt2);
        
%         X(:,:,250:end)=0; %remove
%         X(:,:,1:200)=0; %remove

        [H{1}]=fieldplot3d(1,X,leng,1,Z,'initial time domain',showpictures);
        
        
            power_i_avg=reshape(sum(sum(abs(X).^2,1),2),1,[]);
            %power_i_avg=reshape(mean(mean(abs(X).^2,1),2),1,[]);
            %power_i_avg=reshape(max(max(abs(X).^2,[],1),[],2),1,[]);
            %power_i_c=reshape(abs(X((M-1)/2,(M-1)/2,:)).^2,1,[]);

            figure(331)
            plot(Z,power_i_avg,'linewidth',2);
            %hold all
            %plot(Z,power_i_c);
            %hold off
        
%            break
        
        X=fftshift(fft(X,[],3),3)./sqrt(N2);
        %X=abs(X).*exp(1i*(angle(X)+2*pi));
        
        
        X=prop_TF(X,leng,xlamds,dfl_shift); %-deltaU
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
        
        
        
        [H{222}]=fieldplot3d(2,X,leng,0,Xlamds,'u1 frequency domain',showpictures);
        Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
        [H{111}]=fieldplot3d(1,Xt,leng,1,linspace(0,max(Z),size(Xt,3)),'u1 time domain',showpictures);
 %%
% %         [Mx,My,N]=size(X);
% %         dx=leng/Mx;
% %         dy=leng/My;
% %         xscale=((Mx-1)/2+1-(1:Mx))*dx;
% %         yscale=((My-1)/2+1-(1:My))*dy;
% %         
% %         II=reshape(mean(abs(X).^2,1),My,N);
% %         figure
% %         imagesc(Xlamds,yscale,II,[0 max(max(II))]); % <- problem
% %         figure
% %         plot(Xlamds,II(75,:))
%         
% %         break
        



%         eps=1e-5;
%         ampl=linspace(eps,eps,N);
%         ampl(npfilt1:npfilt2)=1;
%         phs=-KKphase1(K,ampl);
        
%        [H{201}]=fieldplot(201,X,leng,'frequency domain KK',1);
        
        spectrum_i_avg=reshape(mean(mean(abs(X).^2,1),2),1,[]);
        spectrum_i_c=reshape(abs(X((M-1)/2,(M-1)/2,:)).^2,1,[]);
        
        figure(33)
        plot(Xlamds,spectrum_i_avg/max(spectrum_i_avg),'linewidth',1.5);
        ylabel('I [a.u.]');
        xlabel('wavelength [m]');
        set(figure(33),'name','Spectrum','numbertitle','off');
%         hold all
%         plot(Xlamds,spectrum_i_c/max(spectrum_i_c));
%         hold off
        
%%
%        X=X(:,:,npfilt1:npfilt2); %preliminary filtration
         X=X(:,:,npfilt1:npfilt2); %preliminary filtration
%          X(:,:,1:npfiltshift)=0;
%          X(:,:,end-npfiltshift:end)=0;
        
        
        %[H{1}]=fieldplot(1,mean(X,3),leng,'imported real FEL field',showpictures);
        
        if mod(round(interpM*M),2) == 0
            Mn=round(interpM*M)+1;
        else
            Mn=round(interpM*M);
        end
        Xi=zeros(Mn,Mn,Nf,'single');
        
%         X=prop_TF(X,leng,xlamds,-4.5);
%         [H{222}]=fieldplot3d(222,X,leng,0,Xlamds,'U1 frequency domain',showpictures);
            %Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
%             [H{111}]=fieldplot3d(111,Xt,leng,1,linspace(0,max(Z),size(Xt,3)),'U1 time domain',showpictures);
            
        for i=1:Nf
            Xi(:,:,i)=fieldinterpolate(X(:,:,i),leng,1,interpM,interpL,'linear');
        %[H{2}]=fieldplot(2,X,leng_i,'interpolated initial real field');
        end
        X=Xi;
        clear Xi
        
        
        X=prop_TF(X,leng_i,xlamds,U1_to_G+deltaU);
%         X=prop_TF(X,leng_i,xlamds,-4.5);
%         X=prop_TF(X,leng_i,xlamds,U1_to_G+deltaU);
        


        [H{3}]=fieldplot3d(3,X,leng_i,0,Xlamdsf,'field before grating freq',showpictures);
        
     
        
        Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
        [H{31}]=fieldplot3d(31,Xt,leng_i,1,linspace(0,max(Z),size(Xt,3)),'field before grating',showpictures);
        clear Xt
                
    else
        leng_i=leng_i;
        X=fieldgaussian(M*interpM,leng_i,s1,s1,-z1,-z1,xlamds,P0);
        %Xo=fieldgaussian(M,leng,s1,s1,-(z1-U1_to_G),-(z1-U1_to_G),xlamds,P0);
        [H{4}]=fieldplot(4,X,leng_i,'gaussian field just before grating',showpictures);
        Nf=1;
        N2=1;
        dTheta_g_d=0;
        Z=0;
        Xlamdsf=xlamds;
        Xlamds=xlamds;
    end


    
    for i=1:Nf
    X(:,:,i)=grating(X(:,:,i),leng_i,xlamds0,z1_prime,slitpos,R_g_tang,R_g_sag,Theta_g_i,Theta_g_d+dTheta_g_d(i),dTheta_g_d(i),D0,D1,D2,aberrationsincluded);
    
    X(:,:,i)=aperture(X(:,:,i),leng_i,g_length*Theta_g_d,'x'); %turn on
    end

    if roughnessincluded
        X=roughness_1(X, leng_i, Theta_g_i, Theta_g_d, xlamds, 'profile_G.mat');
    end
    
%       Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
%         [H{112}]=fieldplot3d(112,Xt,leng_i,1,linspace(0,max(Z),size(Xt,3)),'field after grating 0',showpictures);
%         clear Xt
         [H{312}]=fieldplot3d(312,X,leng_i,0,Xlamdsf,'field after grating 0 freq',showpictures);

%      X0=X;
     % break
      %%
%       X=prop_TF(X0, leng_i, xlamds,100*slitpos/100); % remove !!!!!!
%         Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
%         [H{113}]=fieldplot3d(113,Xt,leng_i,1,linspace(0,max(Z),size(Xt,3)),'field after grating',showpictures);
%         xlim([3e-5 8.4e-5])
%         clear Xt
%         [H{313}]=fieldplot3d(313,X,leng_i,0,Xlamdsf,'field after grating freq',showpictures);
%         %xlim([1.235e-9 1.245e-9])
    %%   
%       X=prop_TF(X0, leng_i, xlamds,125*slitpos/100); % remove !!!!!!
%         Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
%         [H{114}]=fieldplot3d(114,Xt,leng_i,1,linspace(0,max(Z),size(Xt,3)),'field after grating',showpictures);
%         clear Xt
%         [H{314}]=fieldplot3d(3,X,leng_i,0,Xlamdsf,'field after grating 125 freq',showpictures);
%       X=prop_TF(X0, leng_i, xlamds,200*slitpos/100); % remove !!!!!!
%         Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
%         [H{115}]=fieldplot3d(115,Xt,leng_i,1,linspace(0,max(Z),size(Xt,3)),'field after grating',showpictures);
%         clear Xt
%         [H{315}]=fieldplot3d(315,X,leng_i,0,Xlamdsf,'field after grating 200 freq',showpictures);
        
% break


%     %slit
%      X=prop_TF(X, leng_i, xlamds,slitpos);
%       [H{6}]=fieldplot3d(6,X,leng_i,0,Xlamdsf,'slit frequency domain',showpictures);
% %     [H{6}]=fieldplot3d(6,X1,leng_i,0,Xlamdsf,'slit frequency domain',1);
%      [H{601}]=fieldplot3d(601,ifft(ifftshift(X,3),[],3),leng_i,0,Zf,'slit time domain',1);
% %     clear X1
%     %endslit

        
%         X=prop_TF(X, leng_i, xlamds,slitpos);
%         
%         %%
%     slitinterp=1/5;
%     [H{8}]=fieldplot(8,X,leng_i,'field bi',showpictures);
%     Xi=fieldinterpolate(X,leng_i,1,slitinterp,slitinterp,'spline');
%     Xi=aperture(Xi,leng_i,slitwidth,'x');
%     [H{9}]=fieldplot(9,Xi,leng_i*slitinterp,'field ai',showpictures);
%     [z_opt,fwhm_min,fwhm_0]=waistscan(Xi,leng_i*slitinterp,xlamds,-0.02:0.001:0,2);
%     fprintf('optimal waist is %.0f%% smaller\n',(fwhm_0-fwhm_min)/fwhm_0*100);
%         %%
% %         X=aperture(X,leng_i,slitwidth,'x');
% %         [H{7}]=fieldplot(7,X,leng_i,'field at slit');
% %         X=prop_TF(X, leng_i, xlamds,z3_t_offset-slitpos);
% %         X=mirror(X,leng_i,K,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Theta_m_i,'y',aberrationsincluded);
        
    if roughnessincluded
        X=prop_TF(X, leng_i, xlamds,0.04);
        [H{399}]=fieldplot3d(399,X,leng_i,0,Xlamdsf,'spectrum at M1',showpictures);
        Theta_m1=(Theta_g_i+Theta_g_d)/2;
        X=roughness_1(X, leng_i, Theta_m1, Theta_m1, xlamds, 'profile_M1.mat');
        X=prop_TF(X, leng_i, xlamds,z3_t_offset-0.04);
    else
        X=prop_TF(X, leng_i, xlamds,z3_t_offset);
    end
        

  
    %[H{6}]=fieldplot(6,X,leng_i,'field at mirror 2',1);
    for i=1:Nf
    X(:,:,i)=mirror(X(:,:,i),leng_i,xlamds,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Theta_m_i,'x',aberrationsincluded);
    X(:,:,i)=aperture(X(:,:,i),leng_i,m2_length*Theta_m_i,'x');
    end
    
    [H{400}]=fieldplot3d(400,X,leng_i,0,Xlamdsf,'spectrum at M2',showpictures);
    
    %[H{10}]=fieldplot(10,X,leng_i,'field at M2',showpictures);
    
    if roughnessincluded
        X=roughness_1(X, leng_i, Theta_m_i, Theta_m_i, xlamds, 'profile_M2.mat');
        X=prop_TF(X, leng_i, xlamds,0.13);
        X=roughness_1(X, leng_i, Theta_m_i, Theta_m_i, xlamds, 'profile_M3.mat');       
        X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset-0.13);
    else
        X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset);
    end
    
    %X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset); %remove 0.5!!!
 
    [H{11}]=fieldplot3d(11,X,leng_i,0,Xlamdsf,'spectrum at U2 (not interpolated)',showpictures);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  leng_u2=leng;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %X=prop_TF(X, leng_i, xlamds,u2_prop); %remove?  
if export1to1==1
    M_u2=size(X,1);
    leng_u2=leng_i;
else
    if (size(X,1)~=M_u2 || leng_i~=leng_u2)
        Xi=zeros(M_u2,M_u2,Nf);
        for i=1:Nf
            Xi(:,:,i)=fieldinterpolate(X(:,:,i),leng_i,0,M_u2,leng_u2,'spline');
        end
        X=Xi; clear Xi
    end
end
     X=prop_TF(X, leng_u2, xlamds,u2_prop); %remove?    

    
    [H{12}]=fieldplot3d(12,X,leng_u2,0,Xlamdsf,'spectrum at U2 (interpolated)',showpictures);


    
    if realfield~=1
        break
    end
    
%     for i=1:size(X,1)
%         for j=1:size(X,2)
%             X(i,j,:)=reshape(X(i,j,:),1,[]).*exp(1i*phs(npfilt1:npfilt2));
%         end
%     end

    Xpad=single(zeros(M_u2,M_u2,N2));
    if realfield 
       Xpad(:,:,npfilt1:npfilt2)=X;
%         X=cat(3,zeros(M_u2,M_u2,npfilt1-1),X);
%         X=cat(3,X,zeros(M_u2,M_u2,N2-npfilt2));
    end
    disp(N2);
    disp(size(X,3));
    X=Xpad;
    clear Xpad
    
    [H{121}]=fieldplot3d(121,X,leng_u2,0,Xlamds,'spectrum at U2',showpictures);
        
        spectrum_i_avg_f=reshape(mean(mean(abs(X).^2,1),2),1,[]);
        figure(33)
        hold all
        plot(Xlamds,spectrum_i_avg_f/max(spectrum_i_avg_f),'linewidth',1.5);
        hold off
    
    Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
    %Xf=fftshift(fft(Xt,[],3),3);
    clear X
    
    Xt=Xt.*sqrt(transmission);
    [H{13}]=fieldplot3d(13,Xt,leng_u2,1,Z,'field at U2 absorbed',showpictures);
    %xlim([3e-5 8e-5]);

    
    nshift=double(round(fieldshiftsize/dz));
    Xt=circshift(Xt,[0 0 nshift]);
%     Xt=Xt(:,:,(Npad-1)*N+1:end);
     Xt=Xt(:,:,end-N+1:end);

    [H{132}]=fieldplot3d(132,Xt,leng_u2,1,Z(1:N),'field at circshift',1);
    
    npad=0;
    Ls=0;
    Lf=20;
    ns=find(Z>=Ls*1e-6,1,'first');
    nf=find(Z<=Lf*1e-6,1,'last');

    break
    
%% Save Figures
% folder='U7_1000_horn_3m';
%     h = get(0,'children');
%     %mkdir(['Fig\',nm_f,'\']);
%      mkdir(['Fig\',folder,'\']);
%     for i=1:length(h)
%         if h(i)~=33 & h(i)~=331
%             saveas(h(i), ['Fig\',folder,'\',get(h(i),'Name')], 'fig');
%             %saveas(h(i), ['Fig\',folder,'\',get(h(i),'Name')], 'eps');
%             saveas(h(i), ['Fig\',folder,'\',get(h(i),'Name')], 'png');
%         end
%     end
%% Export

    fieldexport(cat(3,Xt,single(zeros(M_u2,M_u2,npad)).*(1+1i)),[nm_p,'_f7.dfl']);
    toc(t_start)
   % clear variables
    disp('done');
