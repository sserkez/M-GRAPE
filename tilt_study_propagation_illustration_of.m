%Mono

clear all;
%close all;
fclose all;
delete('prop.out.dfl','prop.out');


% resolution without slit calculation
t_start=tic;
showpictures=1;
realfield=1;
onetoone=0;
throughslit=0;
roughnessincluded=1;
aberrationsincluded=1;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);

interpM=3;%2 7
interpL=3;%2 3
Npad=2; %2 1

%export parameters
leng_u2=10e-4;%8e-4
M_u2=201;%151
u2_prop=0; %propagation within U2 (remove, for display only!!!)

deltaU=3.87;
%deltaU=0;

%preliminary filtration
%xlamds0=2.47e-9; %500eV
%dxlamds0=0.0015e-9; %500eV 0.0004e-9

% xlamds0=1.7650e-9; %700eV
% dxlamds0=0.0011e-9; %700eV 0.0004e-9

 xlamds0=1.240e-9; %1000eV
 dxlamds0=0.0008e-9; %1000eV 0.0004e-9
 
 xlamds0=1.3330e-9; %1000eV
 dxlamds0=0.0008e-9; %1000eV 0.0004e-9

disp(' ');
disp(['interpN   =',num2str(interpM)]);
disp(['interpleng=',num2str(interpL)]);

nm_f='500_u1_tdp.out.dfl';
nm_p='500_u1_tdp.out';

%if realfield
    inread;

%    K=2*pi/xlamds;
% else
%     energy=700;
%     M=301;
%     leng=5e-4;
% 
%     xlamds=1239.8/energy/1e9;
%     K=2*pi/xlamds;
% end

P0=10000; %input power for model gaussian
%% Optical system parameters

D0=1123; %[l/mm]
D1=1.6; %[l/mm2]
D2=0.002; %[l/mm3]

% gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)        
% gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) 
gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)     
gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k)

g_length=0.025;         %grating optical length [m]
%R_g_tang=195;           %tangential radius of grating   
R_g_tang=185;           %tangential radius of grating new according to PRD 
R_g_sag=0.18;           %sagittal radius of grating   
%Teta_g_i=1*0.01745;     %incidence angle      
Teta_g_i=1.04*2*pi/360;     %incidence angle new according to PRD   [rad]

m2_length=0.025;
R_m_tang=23.2;          %tangential radius of mirror M2 [m]
Teta_m_i=0.015;         %M2 incidence angle original [rad]
%Teta_m_i=0.02;         %M2 incidence angle original

z3_t_offset=1.53;       %Grating to refocusing mirror distance original
%z3_t_offset=1.59;       %Grating to refocusing mirror distance 

slitpos=1.350;          %grating to slit distance

U1_to_G=1.265;      %First undulator to grating [m]

G_to_U2=3.353;       %Grating to second undulator [m]

gafoc;

% %resolshift=-1.5e-4:2e-5:1.5e-4; %scanned resolution shift via tilt on the grating (=dlambda/lambda)
% resolshift=-0e-4;
% slitwidth=5e-6;
% iterations=size(resolshift,2);
% Power=linspace(0,0,iterations); %input coupling factor

    %dTeta_g_d=xlamds*resolshift(iteration)/Teta_g_d/gr_sigma;
    t_iter=tic;
    %%
    leng_i=leng*interpL;

    if realfield
        
        [X,N]=fieldimport_all(nm_f,M,1); %import field

        X=cat(3,zeros(M,M,N*(Npad-1)),X);
        N2=N*Npad;
        %leng_xy=leng;
        dz=xlamds*zsep;
        dxy=dx;
        leng_z=dz*N2;
        %clear leng

        z=-(N2-1)/2:1:(N2-1)/2;
        xy=-(M-1)/2:1:(M-1)/2;
        Z=(z-min(z))*dz;
        
        k=2*pi/xlamds;
        dk=2*pi/(leng_z);
        K=k+dk*z;
        Xlamds=2*pi./K;
        dxlamds=Xlamds-xlamds;
        
        %preliminary filtration parameters

%         xlamds0=2.481e-9; %500eV
%         dxlamds0=0.001e-9; %500eV
%         xlamds0=1.773e-9; %700eV
%         dxlamds0=0.001e-9; %700eV
        %xlamds0=1.2389e-9;
        %xlamds0=xlamds;

        npfiltshift=0;
        npfilt1=find(Xlamds>(xlamds0+dxlamds0/2),1,'last');%number preliminary filter
        npfilt2=find(Xlamds<(xlamds0-dxlamds0/2),1,'first');%number preliminary filter
        npfilt_mean=mean(npfilt1,npfilt2);
        Nf=npfilt2-npfilt1+1;
        Xlamdsf=Xlamds(npfilt1:npfilt2);
        Zf=Z(npfilt1:npfilt2);
        % end of preliminary filtration parameters

        Teta_g_d_0=acos(cos(Teta_g_i)-xlamds0*D0*1e3);
        z0=linspace(-Nf/2,Nf/2,Nf);
        dTeta_g_d=dk*z0*xlamds0/(k*Teta_g_d_0*gr_sigma);
        %dTeta_g_d=dTeta_g_d(npfilt1:npfilt2);
        
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
        
        [H{2}]=fieldplot3d(2,X,leng,0,Xlamds,'initial frequency domain',showpictures);
        
 
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
        plot(Xlamds,spectrum_i_avg/max(spectrum_i_avg),'linewidth',2);
        hold all
        plot(Xlamds,spectrum_i_c/max(spectrum_i_c));
        hold off
        

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
        
        
        for i=1:Nf
            Xi(:,:,i)=fieldinterpolate(X(:,:,i),leng,1,interpM,interpL,'linear');
        %[H{2}]=fieldplot(2,X,leng_i,'interpolated initial real field');
        end
        X=Xi;
        clear Xi
        
        X=prop_TF(X,leng_i,xlamds,-4.5);
        X=prop_TF(X,leng_i,xlamds,U1_to_G+deltaU);
        


        [H{3}]=fieldplot3d(3,X,leng_i,0,Xlamdsf,'field before grating',showpictures);
        
                
    else
        leng_i=leng_i;
        X=fieldgaussian(M*interpM,leng_i,s1,s1,-z1,-z1,xlamds,P0);
        Xo=fieldgaussian(M,leng,s1,s1,-(z1-U1_to_G),-(z1-U1_to_G),xlamds,P0);
        [H{4}]=fieldplot(4,X,leng_i,'gaussian field just before grating',showpictures);
    end

    for i=1:Nf
    X(:,:,i)=grating(X(:,:,i),leng_i,xlamds0,z1_prime,slitpos,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d+dTeta_g_d(i),dTeta_g_d(i),D0,D1,D2,aberrationsincluded);
    %X(:,:,i)=grating(X(:,:,i),leng_i,xlamds,z1_prime,slitpos,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d,dTeta_g_d(i),D0,D1,D2,aberrationsincluded);
    X(:,:,i)=aperture(X(:,:,i),leng_i,g_length*Teta_g_d,'y');
    end
    
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
%     Xi=aperture(Xi,leng_i,slitwidth,'y');
%     [H{9}]=fieldplot(9,Xi,leng_i*slitinterp,'field ai',showpictures);
%     [z_opt,fwhm_min,fwhm_0]=waistscan(Xi,leng_i*slitinterp,xlamds,-0.02:0.001:0,2);
%     fprintf('optimal waist is %.0f%% smaller\n',(fwhm_0-fwhm_min)/fwhm_0*100);
%         %%
% %         X=aperture(X,leng_i,slitwidth,'y');
% %         [H{7}]=fieldplot(7,X,leng_i,'field at slit');
% %         X=prop_TF(X, leng_i, xlamds,z3_t_offset-slitpos);
% %         X=mirror(X,leng_i,K,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Teta_m_i,'y',aberrationsincluded);
        

        
        
    X=prop_TF(X, leng_i, xlamds,z3_t_offset);
    %[H{6}]=fieldplot(6,X,leng_i,'field at mirror 2',1);
    for i=1:Nf
    X(:,:,i)=mirror(X(:,:,i),leng_i,xlamds,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Teta_m_i,'x',aberrationsincluded);
    X(:,:,i)=aperture(X(:,:,i),leng_i,m2_length*Teta_m_i,'x');
    end
    
    
    %[H{10}]=fieldplot(10,X,leng_i,'field at M2',showpictures);
    X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset); %remove 0.5!!!
 
    [H{11}]=fieldplot3d(11,X,leng_i,0,Xlamdsf,'spectrum at U2 (not interpolated)',showpictures);

    %X=prop_TF(X, leng_i, xlamds,u2_prop); %remove?  

    if size(X,1)~=M_u2 || leng_i~=leng_u2
        Xi=zeros(M_u2,M_u2,Nf);
        for i=1:Nf
        Xi(:,:,i)=fieldinterpolate(X(:,:,i),leng_i,0,M_u2,leng_u2,'spline');
        end
        X=Xi; clear Xi
    end
    
%     X0=X;
%     break
% 
%     %%
%     for u2_prop=0;
    X=prop_TF(X, leng_u2, xlamds,u2_prop); %remove?    
    [H{12}]=fieldplot3d(12,X,leng_u2,0,Xlamdsf,'spectrum at U2 (interpolated)',showpictures);

%     for i=1:size(X,1)
%         for j=1:size(X,2)
%             X(i,j,:)=reshape(X(i,j,:),1,[]).*exp(1i*phs(npfilt1:npfilt2));
%         end
%     end

    Xpad=single(zeros(M_u2,M_u2,N2));
    Xpad(:,:,npfilt1:npfilt2)=X;
    X=Xpad;
    clear Xpad
    
%    [H{121}]=fieldplot3d(121,X,leng_u2,0,Xlamds,'spectrum at U2',showpictures);
        
        spectrum_i_avg_f=reshape(mean(mean(abs(X).^2,1),2),1,[]);
        figure(33)
        hold all
        plot(Xlamds,spectrum_i_avg_f/max(spectrum_i_avg_f),'linewidth',2);
        hold off
    
    Xt=ifft(ifftshift(X,3),[],3).*sqrt(N2);
    %Xf=fftshift(fft(Xt,[],3),3);
    clear X
    
    Xt=Xt.*sqrt(0.05);
    [H{13}]=fieldplot3d(13,Xt,leng_u2,1,Z,'field at U2*5%',showpictures);
%     end
    %%
    fieldshiftsize=7e-6;
    nshift=round(fieldshiftsize/dz);
    Xt=circshift(Xt,[0 0 nshift]);
    Xt=Xt(:,:,(Npad-1)*N+1:end);
    [H{132}]=fieldplot3d(132,Xt,leng_u2,1,Z(1:N),'field at circshift',1);
    
    npad=0;
    Ls=0;
    Lf=20;
    ns=find(Z>=Ls*1e-6,1,'first');
    nf=find(Z<=Lf*1e-6,1,'last');

    break
%% Export
    fieldexport(cat(3,Xt,single(zeros(M_u2,M_u2,npad)).*(1+1i)),'700_u1_filtered_u7_10_201.dfl');
    clear all