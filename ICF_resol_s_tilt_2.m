clear all;
%close all;
fclose all;

% resolution without slit calculation
t_start=tic;
showpictures=1;
realfield=0;
roughnessincluded=0;
aberrationsincluded=1;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);
%showpictures=1;
 idlpl_tilt=(-10:2:10).*1e-5;
 %idlpl_tilt=0e-5;

 inppowerout=linspace(0,0,size(idlpl_tilt,2));
 slitsize=2e-6;

 interpN=8;%12
 interpleng=2;%5
 
 disp(' ');
 disp(['interpN   =',num2str(interpN)]);
 disp(['interpleng=',num2str(interpleng)]);


nm_f='1000.out.dfl';
nm_p='1000.out';

nm='1000_u2.out';
% energy=900;
 P0=1500;
% 
%% Optical system parameters

energy=1000;
xlamds=1239.8/energy*1e-9;

D0=1123; %[l/mm]
D1=1.6; %[l/mm2]
D2=0.002; %[l/mm3]

% gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)         !!!!!!
% gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) !!
gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)         !!!!!!
gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k) !!

g_length=0.025;         %grating optical length [m]
%R_g_tang=195;           %tangential radius of grating     !!!!!!
R_g_tang=185;           %tangential radius of grating new according to PRD !!!!!!
R_g_sag=0.18;           %sagital radius of grating     !!!!!
%Teta_g_i=1*0.01745;     %incidence angle         !!!!!!!
Teta_g_i=1.04*2*pi/360;     %incidence angle new according to PRD        !!!!!!!

m2_length=0.025;
R_m_tang=23.2;          %tangential radius of mirror M2 !!!
Teta_m_i=0.015;         %M2 incidence angle original!!!
%Teta_m_i=0.016;         %M2 incidence angle original!!!
z3_t_offset=1.53;       %Grating to refocusing mirror distance original !!!
%z3_t_offset=1.52;       %Grating to refocusing mirror distance !!!
slitpos=1.350;          %slit position !!!

g_to_und=3.353;
G_to_U2=3.353; %use one of these

deltaU=3.87;
%deltaU=0;

gafoc;

%%

% inread;
%disp(['accuracy=',num2str(interpN/interpleng)]);
leng=8e-4;
M=151;


if mod(interpN*M,2) == 0
    M2=interpN*M+1;
else
    M2=interpN*M;
end
M1=M;
M=M2;
leng2=leng*interpleng;
leng1=leng;
leng=leng2;

dx=leng/M;
dx_prime=dx;
dy=dx;
K=2*pi/xlamds;

for index=1:size(idlpl_tilt,2);
    disp(['iteration # ',num2str(index),' of ',num2str(size(idlpl_tilt,2))]);
    t_iter=tic;
    
    P0=1500;
    
    leng=leng2;
    M=M2;
    dx=leng/M;
    xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;    
    
    dlpl_tilt=idlpl_tilt(index);

dxlamds=xlamds*dlpl_tilt;

dTeta_d=dxlamds/Teta_g_d/gr_sigma;

%                                 wx=w4_s;
%                                 wy=w2_t;
%                                 
%                                 posx=-(-z4_s+z2_t);
%                                 posy=-0;
%                                 
% 
%                             %yoffset=0e-6;
% 
%                                 qx=(1j*pi*(wx^2)/xlamds+posx);
%                                 qy=(1j*pi*(wy^2)/xlamds+posy);
%                                 [xx,yy]=meshgrid((M-1)/2+1-(1:M));
%                                 X=exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dx+yoffset).^2./qy)./2);
%             Xp=X.*conj(X);
%             X=sqrt(Xp*P0/sum(sum(abs(Xp)))).*exp(1i*angle(X));

%                     qx=(1j*pi*(w1_prime^2)/xlamds-z1_prime);
%                     qy=(1j*pi*(w1^2)/xlamds-z1);
%                     [xx,yy]=meshgrid((M-1)/2+1-(1:M));
%                     X=exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dx).^2./qy)./2);
%                     
%                     Xp=X.*conj(X);
%                     X=sqrt(Xp*P0/sum(sum(abs(Xp)))).*exp(1i*angle(X));

if realfield==1

% fd=fopen([nm,'wof.out.dfl'],'r'); %%% !!! changed
fd=fopen(nm_f,'r');

if fd==-1
         error('no such file');
end


dxo=leng1*2/M1;
xscaleo=((M1-1)/2+1-(1:M1))*dxo;
x=fread(fd,2*M1*M1,'double');
b=x(1:2:end)+1i*x(2:2:end);
X=single(reshape(b,M1,M1));
clear x b
minX=min(min(abs(X)));

[xxo,yyo]=meshgrid((M1-1)/2+1-(1:M1));
[xx,yy]=meshgrid((M-1)/2+1-(1:M));
X=interp2(xxo*dxo,yyo*dxo,X,xx*dx,yy*dx,'linear');
X(isnan(X))=minX;
clear xx0 yy0 

X=prop_TF(X,leng,xlamds,2.77);

[sigma,~,~]=mygaussfit(xscale*1e6,abs(X(:,floor(M/2))));
disp(['N=',num2str(sigma*1e-6*sqrt(2)/Teta_g_i/gr_sigma)]);

X=interp2(xx,yy,X,xx,yy.*(Teta_g_i/Teta_g_d),'linear');

leng=leng2;
    M=M2;
    dx=leng/M;
    xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;

elseif realfield==0
    % [xxo,yyo]=meshgrid((Mo-1)/2+1-(1:Mo));
    % [xx,yy]=meshgrid((M-1)/2+1-(1:M));
    % X=interp2(xxo*dxo,yyo*dxo,X,xx*dx,yy.*(Teta_g_i/Teta_g_d)*dx,'linear');
    % X(isnan(X))=0;
    % clear xx0 yy0 

    %   second variant of analytical sourse
% qx=(1j*pi*(w1^2)/xlamds-z1);
% qy=(1j*pi*(w1^2)/xlamds-z1);
% [xx,yy]=meshgrid((M-1)/2+1-(1:M));
% X=exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dx).^2./qy)./2);
% X=interp2(xx,yy,X,xx,yy.*(Teta_g_i/Teta_g_d),'linear');
% 
%  first variant of analytical sourse

    leng=leng2;
    M=M2;
    dx=leng/M;
    xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;

                    qy=(1j*pi*(w1_prime^2)/xlamds-z1_prime);
                    qx=(1j*pi*(w1^2)/xlamds-z1);
                    [xx,yy]=meshgrid((M-1)/2+1-(1:M));
X=single(exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dx).^2./qy)./2));

[sigma,~,~]=mygaussfit(xscale*1e6,abs(X(floor(M/2),:)));
disp(['N=',num2str(sigma*1e-6*sqrt(2)/Teta_g_i/gr_sigma)]);
    
else
    error('invalid realfield value');
end
    

  if dTeta_d~=0
X=abs(X).*exp(1i*(angle(X)+dTeta_d*K*xx.*dx));
  end

 
  
        if showpictures==1 || index==1
                fig3010=figure(3010);
                clf;
                set(fig3010,'name','Intensity_on_gr','numbertitle','off');
                subplot(1,2,1);
                imagesc(xscale,yscale,abs(X));                                   %int-ty
                axis equal tight;
                colorbar('SouthOutside');
                 subplot(1,2,2);
                 imagesc(xscale,yscale,angle(X));                                   %int-ty
                 axis equal tight;
                 colorbar('SouthOutside');
        end


 [sigma,~,~]=mygaussfit(xscale*1e6,abs(X(floor(M/2),:)).^2);
%[sigma,~,~]=mygaussfit(xscale*1e6,abs(X(:,floor(M/2))).^2);
            disp(['horiz waist after grating (2 sigma intensity)=',num2str(sigma*2)]);
           
[sigma,~,~]=mygaussfit(xscale*1e6,abs(X(floor(M/2),:)));
            disp(['horiz waist after grating (w amplitude)=',num2str(sigma*sqrt(2))]);


            
if roughnessincluded==1
[X,~,~]=roughness(X, leng, Teta_g_i, Teta_g_d, K, 125);

    if showpictures==1 || index==1
                    fig3011=figure(3011);
                    clf;
                    set(fig3011,'name','Intensity_on_gr_rough','numbertitle','off');
                    subplot(1,2,1);
                    imagesc(xscale,yscale,abs(X));                                   %int-ty
                    axis equal tight;
                    colorbar('SouthOutside');
                     subplot(1,2,2);
                     imagesc(xscale,yscale,angle(X));                                   %int-ty
                     axis equal tight;
                     colorbar('SouthOutside');
    end
end

% X(:,:,i)=grating(X(:,:,i),leng_i,xlamds0,z1_prime,slitpos,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d+dTeta_g_d(i),dTeta_g_d(i),D0,D1,D2,aberrationsincluded);
X=grating(X,leng,xlamds,z1_prime,slitpos,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d,0,D0,D1,D2,aberrationsincluded);
X=aperture(X,leng,g_length*Teta_g_d,'y');

% % % X=prop_TF(X.', leng, xlamds,0.04);
% % % 
% % % if roughnessincluded==1
% % % [X,XX,hh]=roughness(X, leng, (Teta_g_i+Teta_g_d)/2, (Teta_g_i+Teta_g_d)/2, K, 15);
% % % 
% % % 
% % % if showpictures==1
% % % figure(5012)
% % % plot(XX,hh,'k')
% % % hold all
% % % plot(xscale./tan(Teta_g_d/2),abs(X(:,floor(M/2))).^2./max(abs(X(:,floor(M/2))).^2).*max(hh));
% % % hold off
% % % title(['FFT RMS ' num2str(std(hh)) '   Mean ' num2str(mean(hh))]);
% % % grid on
% % % end
% % % end
% % % 
% % % 
% % %         if showpictures==1 || index==1
% % %                 fig3012=figure(3012);
% % %                 clf;
% % %                 set(fig3012,'name','Intensity_on_mirr','numbertitle','off');
% % %                 subplot(1,2,1);
% % %                 imagesc(xscale,yscale,abs(X));                                   %int-ty
% % %                 axis equal tight;
% % %                 colorbar('SouthOutside');
% % %                  subplot(1,2,2);
% % %                  imagesc(xscale,yscale,angle(X));                                   %int-ty
% % %                  axis equal tight;
% % %                  colorbar('SouthOutside');
% % %         end
% % %         
% % %         
% % % X=prop_TF(X, leng, xlamds,slitpos-0.04);
% % % 
% % % 
% % %                       xp=X.*conj(X);
% % % X=sqrt(xp*P0/sum(sum(abs(xp)))).*exp(1i*angle(X));
% % %                       clear xp
% % % 
% % %     xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;
% % % 
% % % 
% % %             if showpictures==1
% % %             fig302=figure(302);
% % %             clf;
% % %             set(fig302,'name','Intensity_slit','numbertitle','off');
% % %             subplot(1,2,1);
% % %             imagesc(xscale,yscale,abs(X));                                   %int-ty
% % %             axis equal tight;
% % %             colorbar('SouthOutside');
% % %              subplot(1,2,2);
% % %              imagesc(xscale,yscale,angle(X));                                   %int-ty
% % %              axis equal tight;
% % %              colorbar('SouthOutside');
% % %             end
% % % 
% % %  %SLIT%%%%%%%%%%%%%%%%%%%%%%%%
% % %             [sigma,~,~]=mygaussfit(xscale*1e6,abs(X(:,floor(M/2))));
% % %             %sigma
% % %             disp(['waist on slit (w, amplpitude)=',num2str(sigma*sqrt(2))]);
% % %             
% % %             [sigma,mu,A]=mygaussfit(xscale*1e6,abs(X(:,floor(M/2))).^2);
% % %             %sigma
% % %             sigma=sigma/1e6;
% % %             mu=mu/1e6;
% % %             disp(['waist on slit (2 sigma intensity)=',num2str(sigma*2)]);
% % %             disp(['resolution (sigma->FWHM)=',num2str(1/(sigma*2*sqrt(2*log(2))/(slitpos.*xlamds./Teta_g_d./gr_sigma)))]);
% % %             
% % %             
% % % slitint=abs(X(:,floor(M/2))).^2;
% % % slitintscale=xscale;
% % %             
% % %             slitFWHM=findFWHM(xscale*1e6,slitint);
% % %             slitFWHM=abs(slitFWHM/1e6);
% % %             disp(['resolution (FWHM)=',num2str(1/(slitFWHM/(slitpos.*xlamds./Teta_g_d./gr_sigma)))]);
% % % 
% % % 
% % %             
% % %             if slitsize~=0
% % %             distr_ft=double(abs(slitintscale)<slitsize/2);
% % %             distr_ft=distr_ft/sum(distr_ft);
% % %             slitint_conv=conv(slitint,distr_ft);
% % %             slitint_conv=slitint_conv((numel(slitint_conv)-numel(slitint))/2+1:(numel(slitint_conv)+numel(slitint))/2);
% % % 
% % %             
% % %             slitFWHM=findFWHM(xscale*1e6,slitint_conv);
% % %             slitFWHM=abs(slitFWHM/1e6);
% % %             disp(['resolution (FWHM_convolved with ',num2str(slitsize*1e6),' um)=',num2str(1/(slitFWHM/(slitpos.*xlamds./Teta_g_d./gr_sigma)))]);
% % %             end
% % %             
% % %             if showpictures==1
% % %                 if realfield
% % %                     figure(9911)
% % %                 else
% % %                     figure(991)
% % %                 end
% % %             plot(slitintscale,slitint);
% % %             hold all
% % %                 if slitsize~=0
% % %                 plot(slitintscale,slitint_conv);
% % %                 plot(slitintscale,distr_ft);
% % %                 end
% % %             %plot(slitintscale,A * exp( -(slitintscale-mu).^2 ./ (2*sigma^2) ));
% % %             hold off
% % %             xlim([-2e-5 2e-5]);
% % %             end
% % %    
% % % %ENDSLIT%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % 
% % %             %sum(sum(abs(X)))
% % % 
% % %                         xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;
% % % 
% % %             %islitsize(index)=dx*floor((slitsize-dx)/dx)+dx;
% % %             %slitsize=islitsize(index);
% % %             %X=slit(X,leng,slitsize,'y');
% % %             %islitsize(index)=dx*floor(slitsize/dx);
% % % X=prop_TF(X, leng, xlamds, z3_t_offset-slitpos); %z3_t

X=prop_TF(X, leng, xlamds,z3_t_offset);

% % % X=lens_tang(X,leng,K,z3_t,7.2-2.77-z3_t_offset,f_m_t,R_m_tang,R_m_tang,Teta_m_i,Teta_m_i,'y',aberrationincluded); %f_m_t
X=mirror(X,leng,xlamds,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Teta_m_i,'x',aberrationsincluded);
X=aperture(X,leng,m2_length*Teta_m_i,'x');
    
    
            if index==1 && showpictures==1
            fig1141=figure(1141);
            clf;
            set(fig1141,'name','Intensity_M2','numbertitle','off');
            imagesc(xscale,yscale,abs(X));                                   %int-ty
            axis equal tight;
            %title(['Intensity, z= ',num2str(Z),' m']);
            colorbar('SouthOutside');
            end

% % % X=prop_TF(X, leng, xlamds,7.2-2.77-z3_t_offset); %z4_t
X=prop_TF(X, leng, xlamds,G_to_U2-z3_t_offset); %remove 0.5!!!
%break

    if showpictures==1
    fig1141=figure(1141);
    clf;
    set(fig1141,'name','Intensity_out','numbertitle','off');
    imagesc(xscale,yscale,abs(X));                                   %int-ty
    axis equal tight;
    %title(['Intensity, z= ',num2str(Z),' m']);
    colorbar('SouthOutside');
    end

        % 
        % % fig115=figure(115);
        % % clf;
        % % set(fig115,'name','Phase_11','numbertitle','off');
        % % imagesc(xscale,yscale,angle(X));                               %phase
        % % axis equal tight;
        % % %    title(['Phase, z= ',num2str(Z),' m']);
        % % colorbar('SouthOutside');

    if interpleng~=1 || interpN~=1
    P0=sum(sum(abs(X(round(M2/2*(1-1/interpleng)):round(M2-M2/2*(1-1/interpleng)),round(M2/2*(1-1/interpleng)):round(M2-M2/2*(1-1/interpleng)))).^2));
     M=M1;
     [XI,YI]=meshgrid((M-1)/2+1-(1:M));
     dx1=leng1/M;
     Xi=interp2(xx*dx,yy*dx,X,XI.*dx1,YI.*dx1);
                Xp=Xi.*conj(Xi);
                X=sqrt(Xp.*P0./sum(sum(abs(Xp)))).*exp(1i*angle(Xi));
     clear Xi Xp
    end
    
        [xx,yy]=meshgrid((M-1)/2+1-(1:M));

        %sum(sum(abs(X).^2))

        M=M1; leng=leng1;
        dx=leng/M;
        xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;

        if showpictures==1
        fig1142=figure(1142);
                clf;
                set(fig1142,'name','Intensity_out_i','numbertitle','off');
                subplot(1,2,1);
                imagesc(xscale,yscale,abs(X));
                axis equal tight;
                %title(['Intensity, z= ',num2str(Z),' m']);
                colorbar('SouthOutside');

                subplot(1,2,2);
                imagesc(xscale,yscale,angle(X));
                axis equal tight;
                %title(['Intensity, z= ',num2str(Z),' m']);
                colorbar('SouthOutside');
        end

% fig115=figure(115);
% clf;
% set(fig115,'name','Phase_11','numbertitle','off');
% imagesc(xscale,yscale,angle(X));
% axis equal tight;
% %    title(['Phase, z= ',num2str(Z),' m']);
% colorbar('SouthOutside');
% break

x3=reshape(X,1,M*M);
x4=zeros(1,2*M*M);
x4(1:2:end)=real(x3);
x4(2:2:end)=imag(x3);
fd=fopen('prop_field.dfl','wb');
fwrite(fd,x4,'double');
fclose(fd);
clear x3 x4


dos('genesis301');

fd=fopen([nm,'.dfl'],'r');

nm_p=nm;
inread;

%      rbeam=sqrt(rxbeam^2+rybeam^2);
%      ray=sqrt(zrayl*xlamds/pi);
%      if dgrid==0
%         leng=rmax0*(rbeam+ray); %%%! *2
%      else
%         leng=dgrid*2;
%      end
%      dx=leng/M;     dy=dx;
%      K=2*pi/xlamds;

     if fd==-1
         error('no such file');
     end

     x=fread(fd,2*M*M,'double');
     b=x(1:2:end)+1i*x(2:2:end);
     x1=reshape(b,M,M);
     clear x b
     
     inppower=sum(sum(x1.*conj(x1)));
     inppowerout(index)=inppower;

I1=abs(x1);

if showpictures==1
fig2=figure(107);
clf;
set(fig2,'name','Intensity_GENESIS','numbertitle','off');
imagesc(xscale,yscale,I1);
axis equal tight
colorbar('SouthOutside');
%title(['Intensity, z= 0m']);
V=axis;

% 29.04.2013
%text(double(V(1)),double(V(3)),sprintf('  \\lambda = %.3e\n  P_{sum} = %.3e', xlamds, inppower),...
%       'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10,'color','white');

end
   
figresol=figure(109);
clf
set(figresol,'name','resolution without slit','numbertitle','off');
%plot(iyoffset,inppowerout/max(inppowerout));
%plot(iyoffset*gr_sigma/slitpos/xlamds,inppowerout/max(inppowerout));
plot(idlpl_tilt(1:index),inppowerout(1:index));
%findFWHM(iyoffset./xlamds./slitpos.*Teta_g_d.*gr_sigma,inppowerout);
time=toc(t_iter);
disp(['time per iteration  ',num2str(time),' sec.']);
disp(['estimated finish    ',num2str(time*(size(idlpl_tilt,2)-index)/60),' min.']);
end


clear X x1 xx yy

% figresol=figure(109);
% clf
% %set(figresol,'name','resolution without slit','numbertitle','off');
% %plot(iyoffset,inppowerout/max(inppowerout));
% %plot(iyoffset*gr_sigma/slitpos/xlamds,inppowerout/max(inppowerout));
%  plot(islitsize,inppowerout./max(inppowerout));

resolution=1/findFWHM(idlpl_tilt,inppowerout/max(inppowerout));
% 
% zzexp_arg=(idlpl_tilt)'.*1e5;
%    zzexp_func=(inppowerout)';
% toc(t_start)

 zzexp_arg=slitintscale.';
 zzexp_func=slitint;

disp(' ');
disp(['done, resolution=',num2str(resolution)]);