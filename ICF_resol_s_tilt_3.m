%% Propagation parameters

%Theta_g_d=0;
clear all;
close all;
fclose all;
delete('prop.out.dfl','prop.out');

t_start=tic;
showpictures=0;
realfield=1;
onetoone=0;
throughslit=0;
roughnessincluded=1;
aberrationsincluded=1;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);

nm_p='500_u1.out';
nm_f=[nm_p,'.dfl'];
nm_inp_u2='geninp_500_u2';
nm_p_u2='500_u2.out';
nm_f_u2=[nm_p_u2,'.dfl'];


 idlpl_tilt=(-25:2:25).*1e-5;
 %idlpl_tilt=(-8:2:8).*1e-5;
 %idlpl_tilt=0e-5;
 inppowerout=linspace(0,0,size(idlpl_tilt,2));
 slitsize=2e-6;

interpM=4;%2 8 3 6
interpL=3;%2 3 4 2

%export parameters
leng_u2=3e-4;%10e-4
M_u2=201;%201


u2_prop=0; %propagation within U2 (remove, for display only!!!)

deltaU=3.87;
%deltaU=0;

disp(' ');
disp(['interpN   =',num2str(interpM)]);
disp(['interpL=',num2str(interpL)]);


dos(['COPY ',nm_inp_u2,' geninp']);

inread;

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
%Theta_g_i=1*0.01745;     %incidence angle      
Theta_g_i=1.04*2*pi/360;     %incidence angle new according to PRD   [rad]

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
% break
    
if mod(interpM*M,2) == 0
    M2=interpM*M+1;
else
    M2=interpM*M;
end

leng_i=leng*interpL;

dxlamds=xlamds*idlpl_tilt;
dTheta_g_d=dxlamds/Theta_g_d/gr_sigma;

for index=1:size(idlpl_tilt,2);
    disp(['iteration # ',num2str(index),' of ',num2str(size(idlpl_tilt,2))]);
    t_iter=tic;
    
       

        [X,N]=fieldimport_all(nm_f,M,1);
        X=prop_TF(X,leng,xlamds,-4.8);
        [H{1}]=fieldplot3d(1,X,leng,0,xlamds,'U1',showpictures);
        X=prop_TF(X,leng,xlamds,U1_to_G+deltaU);
        [H{2}]=fieldplot3d(2,X,leng,0,xlamds,'before G',showpictures);
        break
        X=fieldinterpolate(X,leng,1,interpM,interpL,'cubic');
        X=grating(X,leng_i,xlamds+dxlamds(index),z1_prime,slitpos,R_g_tang,R_g_sag,Theta_g_i,Theta_g_d+dTheta_g_d(index),dTheta_g_d(index),D0,D1,D2,aberrationsincluded);
        X=aperture(X,leng_i,g_length*Theta_g_d+dTheta_g_d(index),'x');
        if roughnessincluded
            X=roughness_1(X, leng_i, Theta_g_i, Theta_g_d+dTheta_g_d(index), xlamds, 'profile_G.mat');
        end
        
        [H{3}]=fieldplot3d(3,X,leng_i,0,xlamds,'after G',showpictures);
        
        if roughnessincluded
            X=prop_TF(X, leng_i, xlamds,0.04);
            [H{399}]=fieldplot3d(399,X,leng_i,0,xlamds,'spectrum at M1',showpictures);
            Theta_m1=(Theta_g_i+Theta_g_d+dTheta_g_d(index))/2;
            X=roughness_1(X, leng_i, Theta_m1, Theta_m1, xlamds, 'profile_M1.mat');
            X=prop_TF(X, leng_i, xlamds,z3_t_offset-0.04);
        else
            X=prop_TF(X, leng_i, xlamds,z3_t_offset);
        end
        
        X=mirror(X,leng_i,xlamds,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Teta_m_i,'x',aberrationsincluded);
        X=aperture(X,leng_i,m2_length*Teta_m_i,'x');
        [H{4}]=fieldplot3d(4,X,leng_i,0,xlamds,'after M2',showpictures);

        if roughnessincluded
            X=roughness_1(X, leng_i, Teta_m_i, Teta_m_i, xlamds, 'profile_M2.mat');
            X=prop_TF(X, leng_i, xlamds,0.13);
            X=roughness_1(X, leng_i, Teta_m_i, Teta_m_i, xlamds, 'profile_M3.mat');       
            X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset-0.13);
        else
            X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset);
        end
        
        [H{5}]=fieldplot3d(5,X,leng_i,0,xlamds,'U2',showpictures);
        
        if (size(X,1)~=M_u2 || leng_i~=leng_u2)            
                X=fieldinterpolate(X,leng_i,0,M_u2,leng_u2,'cubic');
        end
        
        X=X*sqrt(0.05)/100;
        
        [H{12}]=fieldplot3d(12,X,leng_u2,0,xlamds,'U2_i',1);
        
        fieldexport(X,'field_prop.dfl');

        dos(['genesis301']);
        
        [X,N]=fieldimport_all(nm_f_u2,M_u2,1);
        [H{15}]=fieldplot3d(15,X,1,0,xlamds,'after ampl',1);
        inppowerout(index)=sum(sum(abs(X).^2));
        
        
        figresol=figure(109);
        clf
        set(figresol,'name','resolution without slit','numbertitle','off');
        plot(idlpl_tilt(1:index),inppowerout(1:index));
        time=toc(t_iter);
        disp(['time per iteration  ',num2str(time),' sec.']);
        disp(['estimated finish    ',num2str(time*(size(idlpl_tilt,2)-index)/60),' min.']);
        
end
%%
[X,N]=fieldimport_all(nm_f,M,1);
X=prop_TF(X,leng,xlamds,-4.8);
X=fieldinterpolate(X,leng,0,M_u2,leng_u2,'cubic');
X=X*sqrt(0.05)/100;
[H{16}]=fieldplot3d(16,X,leng_u2,0,xlamds,'1to1',1);
fieldexport(X,'field_prop.dfl');
dos(['genesis301']);
[X,N]=fieldimport_all(nm_f_u2,M_u2,1);
[H{17}]=fieldplot3d(17,X,leng_u2,0,xlamds,'1to1 after ampl',1);
inppowerout_1to1=sum(sum(abs(X).^2));
ICF=max(inppowerout)/inppowerout_1to1;
clear X x1 xx yy

% figresol=figure(109);
% clf
% %set(figresol,'name','resolution without slit','numbertitle','off');
% %plot(iyoffset,inppowerout/max(inppowerout));
% %plot(iyoffset*gr_sigma/slitpos/xlamds,inppowerout/max(inppowerout));
%  plot(islitsize,inppowerout./max(inppowerout));

resolution=1/findFWHM(idlpl_tilt,inppowerout/max(inppowerout));
% 
zzexp_arg=(idlpl_tilt)'.*1e5;
   zzexp_func=(inppowerout)';
toc(t_start)

figure(109);
text(1,1,sprintf('Resolution = %.0f\nICF = %.3f', resolution, ICF),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');

%  zzexp_arg=slitintscale.';
%  zzexp_func=slitint;

disp(' ');
disp(['done, resolution=',num2str(resolution)]);


%%
load('C:\-D-\Work\LCLS\New_matlab\Steadystate\500_8_new.mat'); % delete!!!!!!!!!
k0=2*pi./(xlamds*(1+idlpl_tilt));
k0l=2*pi./(xlamds*(1+20*idlpl_tilt));
k=linspace(k0l(1),k0l(end),20*numel(k0));
I_function=inppowerout;

% I_function=I_function-max(I_function)*0.02;
% I_function(I_function<0)=0;
% I_function=I_function+max(I_function)*0.02;

I_function=I_function/max(I_function)*0.05*ICF;
I_function=interp1([k0l(1) k0 k0l(end)],[I_function(1)/5 I_function I_function(end)/5],k,'cubic');
[Phase,haxes]=KKphase1(k,I_function);
xlim(haxes(1),[idlpl_tilt(1) idlpl_tilt(end)]);
xlim(haxes(2),[idlpl_tilt(1) idlpl_tilt(end)]);
ylim(haxes(1),[0 0.12]);
set(haxes(2),'YTick',-3:1:3)
set(haxes(1),'YTick',0:0.01:0.12)
set(figure(1001), 'Position', [100, 100, 550, 450]);
%%
% k0=2*pi./(xlamds*(1+idlpl_tilt));
% %k0l=2*pi./(xlamds*(1+2*idlpl_tilt));
% k=linspace(k0(1),k0(end),numel(k0));
% I_function=inppowerout;
% 
% % I_function=I_function-max(I_function)*0.02;
% % I_function(I_function<0)=0;
% % I_function=I_function+max(I_function)*0.02;
% 
% 
% I_function=I_function/max(I_function)*0.058;
% I_function=interp1(k0,I_function,k,'cubic');
% Phase=KKphase1(k,I_function);
% %xlim[idlpl_tilt(1) idlpl_tilt(end)]
% % figure(576)
% % plotyy(2*pi./k,I_function,2*pi./k,Phase,'linewidth',2);
 break
%% save
save('500_7_new.mat', 'idlpl_tilt', 'inppowerout', 'ICF', 'Phase');