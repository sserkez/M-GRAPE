%% Propagation parameters
cd('c:\-D-\Work\SASE3_SXRSS\ICF\');
%%
%Theta_g_d=0;
clear all;
close all;
fclose all;
delete('prop.out.dfl','prop.out');

runrun=0;
scan_und=0;
zz=[-1:0.2:10];

t_start=tic;
showpictures=1;
realfield=1;
onetoone=0;
throughslit=1;
roughnessincluded=0;
r_files_fldr='c:\-D-\Work\LCLS\New_matlab\';
aberrationsincluded=1;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);


transmission=0.111;
%energy=300;
nm_icf_mat='c:\-D-\Work\SASE3_SXRSS\ICF\ICF_700.mat';
nm_p=['c:\-D-\Work\SASE3_SXRSS\700_u1.out'];
nm_f=[nm_p,'.dfl'];

nm_inp_u2='c:\-D-\Work\SASE3_SXRSS\ICF\SASE3_700_u2.in';
nm_p_u2='c:\-D-\Work\SASE3_SXRSS\ICF\700_u2.out';
nm_f_u2=[nm_p_u2,'.dfl'];
dos(['COPY ',nm_inp_u2,' c:\-D-\Work\SASE3_SXRSS\ICF\geninp']);

%xlamds0=1.781e-9; %700eV 

 idlpl_tilt=(-20:1:20).*1e-5;

 %idlpl_tilt=(-8:2:8).*1e-5;
 idlpl_tilt=200e-5;
  idlpl_tilt=-9e-5;
 inppowerout=linspace(0,0,size(idlpl_tilt,2));
 slitsize=2e-6;

interpM=8;%2 8 3 6
interpL=6;%2 3 4 2

%export parameters
leng_u2=8e-4;%10e-4
M_u2=201;%201
leng_u2=3e-4;%10e-4
M_u2=201;%201

 u2_prop=0; %propagation within U2 (remove, for display only!!!)
% u2_prop=z4_t;


disp(' ');
disp(['interpN   =',num2str(interpM)]);
disp(['interpL=',num2str(interpL)]);


%



%% Optical system parameters

d=outread(nm_p,1,1,2); 
leng=d.inp.leng;
M=d.inp.ncar;
xlamds=d.inp.xlamds;
zsep=d.inp.zsep;
dx=leng/M;
% if (xlamds-xlamds0)/xlamds>0.01
%     error ('wrong central wavelength defined')
% end

opt_parms_5;

dfl_shift=-1.09;

g_length=0.025;         %grating optical length [m]
m2_length=0.025;
%remove remove  remove!!!
% g_length=0.05;         %grating optical length [m]
% m2_length=0.05;

energy=1239.8411/xlamds*1e-9;

gafoc_XFEL;
%C_20_m=-Theta_m_i/R_m_tang/(z3_t_offset-slitpos)+1/(z3_t_offset-slitpos)^2-Theta_m_i/R_m_tang/z4_t+1/z4_t^2% break
%u2_prop=z4_t; %REMOVEREMOVEREMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
        X=prop_TF(X,leng,xlamds,dfl_shift);
        [H{1}]=fieldplot3d(1,X,leng,0,xlamds,'U1',showpictures);
        X=fieldinterpolate(X,leng,1,interpM,interpL,'spline');
        X=prop_TF(X,leng_i,xlamds,U1_to_G+deltaU);
        [H{2}]=fieldplot3d(2,X,leng_i,0,xlamds,'before G',showpictures);
        
%         figure(6575)
%         plot(linspace(-leng_i/2,leng_i/2,size(X,1))./Theta_g_i,sum(abs(X).^2,1));
%        break
        
        X=grating(X,leng_i,xlamds+dxlamds(index),z1_prime,slitpos,R_g_tang,R_g_sag,Theta_g_i,Theta_g_d+dTheta_g_d(index),dTheta_g_d(index),D0,D1,D2,aberrationsincluded);
        X=aperture(X,leng_i,g_length*Theta_g_d+dTheta_g_d(index),'x');
        if roughnessincluded && 1
            X=roughness_1(X, leng_i, Theta_g_i, Theta_g_d+dTheta_g_d(index), xlamds, [r_files_fldr,'profile_G.mat']);
        end
        
        [H{3}]=fieldplot3d(3,X,leng_i,0,xlamds,'after G',showpictures);
        
%         figure(6576)
%         plot(linspace(-leng_i/2,leng_i/2,size(X,1))./Theta_g_d,sum(abs(X).^2,1));
%         break
        
%         X1=fieldinterpolate_a(X,leng_i,leng_i,1,11,0.09,1,1,'spline');
%         X1=prop_TF_a(X1, leng_i,leng_i, xlamds,slitpos);
%         %waistscan_1(68,X1,leng_i,xlamds,[-0.001:0.001:0.001]);
%         %waistscan_1_a(68,X1,leng_i,leng_i,xlamds,[-0.02:0.002:0.02])
%         %break
%         [H{5213}]=fieldplot3d_a(5213,X1,leng_i,leng_i,0,xlamds,'spectrum slit',showpictures);
% 
%        figure(6576)
%         plot(linspace(-leng_i/2,leng_i/2,size(X1,2))/slitpos*Theta_g_d*gr_sigma./xlamds,sum(abs(X1).^2,1));
%         1/findFWHM(linspace(-leng_i/2,leng_i/2,size(X1,2))/slitpos*Theta_g_d*gr_sigma./xlamds,sum(abs(X1).^2,1))
%         break
        
        
        if roughnessincluded && 1
            X=prop_TF(X, leng_i, xlamds,0.06);
            [H{399}]=fieldplot3d(399,X,leng_i,0,xlamds,'spectrum at M1',showpictures);
            Theta_m1=(Theta_g_i+Theta_g_d+dTheta_g_d(index))/2;
            X=roughness_1(X, leng_i, Theta_m1, Theta_m1, xlamds, [r_files_fldr,'profile_M1.mat']);
            X=prop_TF(X, leng_i, xlamds,z3_t_offset-0.06);
        else
            X=prop_TF(X, leng_i, xlamds,z3_t_offset);
        end
        
        
        
        
        if showpictures
%             interpleng_x_s=0.07;
%             interpleng_y_s=0.3;
%             interpM_x_s=4;
%             interpM_y_s=0.05;
%             leng_x_s=interpleng_x_s*leng_i;
%             leng_y_s=interpleng_y_s*leng_i;
%             Xt=fieldinterpolate_a(X,leng,leng,1,interpM_x_s,interpM_y_s,interpleng_x_s,interpleng_y_s,'spline');
%             [H{5}]=fieldplot3d_a(40,Xt,leng_x_s,leng_y_s,0,xlamds,'before M2',showpictures);
%             %%
%             Xt1=prop_TF_a(Xt, leng_x_s,leng_y_s, xlamds,-z3_t_offset+slitpos);
%             [H{4}]=fieldplot3d_a(30,Xt1,leng_x_s,leng_y_s,0,xlamds,'at slit',showpictures);
%             waistscan_1_a(346,Xt1,leng_x_s,leng_y_s,xlamds,[-0.05:0.005:0.05]);
%             clear Xt1
%             
        end
        
        

        
        
        
        [H{5}]=fieldplot3d(40,X,leng_i,0,xlamds,'before M2',showpictures);
        X=mirror(X,leng_i,xlamds,z3_t_offset-slitpos,z4_t,R_m_tang,R_m_tang,Theta_m_i,'x',aberrationsincluded);
        X=aperture(X,leng_i,m2_length*Theta_m_i,'x');
        [H{6}]=fieldplot3d(4,X,leng_i,0,xlamds,'after M2',showpictures);

%         X=prop_TF(X, leng_i, xlamds,0.3528);
        
        if roughnessincluded && 1
            X=roughness_1(X, leng_i, Theta_m_i, Theta_m_i, xlamds, [r_files_fldr,'profile_M2.mat']);
            X=prop_TF(X, leng_i, xlamds,0.13);
            X=roughness_1(X, leng_i, Theta_m_i, Theta_m_i, xlamds, [r_files_fldr,'profile_M3.mat']);       
            X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset+u2_prop-0.13);
        else
            X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset+u2_prop);
        end
        
        [H{5}]=fieldplot3d(5,X,leng_i,0,xlamds,'U2+u2_prop',showpictures);

        
% calculation of resolving power

% [My,Mx,N]=size(X);
% Xx=X(round((Mx+1)/2),:,:)';
% Xy=X(:,round((My+1)/2),:);
% Xx=abs(Xx).^2;
% Xy=abs(Xy).^2;
% Xx=Xx./sum(sum(Xx));
% Xy=Xy./sum(sum(Xy));
% %     Ix=sum(sum(abs(X).^2,3),1)'./sum(sum(sum(abs(X).^2))); %projected
% %     Iy=sum(sum(abs(X).^2,3),2)./sum(sum(sum(abs(X).^2)));
% Ix=Xx;
% Iy=Xy;
% figure(4675)
% % disp1=xlamds/Theta_g_d/gr_sigma;
% % disp2=disp1*(z3_t_offset-f_m_t)/f_m_t;
% % disp3=disp2*(G_to_U2-(f_m_t*z3_t_offset/(z3_t_offset-f_m_t))+u2_prop);
% x=z3_t_offset*f_m_t/(z3_t_offset-f_m_t);
% multiplier=f_m_t*Theta_g_d*gr_sigma/(G_to_U2-z3_t_offset-x+u2_prop)/(z3_t_offset-f_m_t);
% 
% % theta_2=dTheta_g_d*(z3_t_offset-f_m_t)/f_m_t;
% % theta_2*(G_to_U2-(f_m_t*z3_t_offset/(z3_t_offset-f_m_t)));
% plot(linspace(-leng_i/2,leng_i/2,size(X,2)).*multiplier./xlamds,Ix);
% 1/findFWHM(linspace(-leng_i/2,leng_i/2,size(X,2)).*multiplier./xlamds,Ix)
% break
        
        if (size(X,1)~=M_u2 || leng_i~=leng_u2)            
                X=fieldinterpolate(X,leng_i,0,M_u2,leng_u2,'spline');
        end
       
        X=X*sqrt(0.05)/100;


    
        
        [H{12}]=fieldplot3d(12,X,leng_u2,0,xlamds,'U2_i',1);
if runrun==0
break
end
        fieldexport(X,'c:\-D-\Work\SASE3_SXRSS\ICF\field_prop.dfl');

        dos(['c:\-D-\Work\SASE3_SXRSS\ICF\genesis301.exe']);
        
        [X1,N]=fieldimport_all(nm_f_u2,M_u2,1);
        [H{15}]=fieldplot3d(15,X1,1,0,xlamds,'after ampl',1);
        inppowerout(index)=sum(sum(abs(X1).^2));
        
        
        figresol=figure(109);
        clf
        set(figresol,'name','resolution without slit','numbertitle','off');
        plot(idlpl_tilt(1:index),inppowerout(1:index));
        time=toc(t_iter);
        disp(['time per iteration  ',num2str(time),' sec.']);
        disp(['estimated finish    ',num2str(time*(size(idlpl_tilt,2)-index)/60),' min.']);
        if index~=1
            xlim([idlpl_tilt(1) idlpl_tilt(index)]);
        end
end
%%
break
if scan_und==1
    %X=prop_TF(X,d(Di).inp.leng,d(Di).inp.xlamds,-(12-11.28));
    waistscan_1(68,X,leng_u2,xlamds,zz);
    h1=figure(68);
    for i=1:3
    subplot(3,1,i)
    xlim([zz(1) zz(end)]);
    end
end


if runrun==0
break
end
[X1,N]=fieldimport_all(nm_f,M,1);
X1=prop_TF(X1,leng,xlamds,dfl_shift);
X1=fieldinterpolate(X1,leng,0,M_u2,leng_u2,'cubic');
X1=X1*sqrt(0.05)/100;
[H{16}]=fieldplot3d(16,X1,leng_u2,0,xlamds,'1to1',1);
fieldexport(X1,'field_prop.dfl');
dos(['genesis301']);
[X1,N]=fieldimport_all(nm_f_u2,M_u2,1);
[H{17}]=fieldplot3d(17,X1,leng_u2,0,xlamds,'1to1 after ampl',1);
inppowerout_1to1=sum(sum(abs(X1).^2));

clear X1 x1 xx yy

[X1,N]=fieldimport_all(nm_f,M,1);
X1=prop_TF(X1,leng,xlamds,dfl_shift-2);
X1=fieldinterpolate(X1,leng,0,M_u2,leng_u2,'spline');
X1=X1*sqrt(0.05)/100;
[H{16}]=fieldplot3d(16,X1,leng_u2,0,xlamds,'1to1',1);
fieldexport(X1,'field_prop.dfl');
dos(['genesis301']);
[X1,N]=fieldimport_all(nm_f_u2,M_u2,1);
[H{17}]=fieldplot3d(17,X1,leng_u2,0,xlamds,'1to1 after ampl',1);
inppowerout_best=sum(sum(abs(X1).^2));

clear X1 x1 xx yy



% inppowerout_1to1 =1.5044e+05; %700
% inppowerout_best=2.7608e+05; %700

%ICF=0;
ICF=max(inppowerout)/inppowerout_1to1;
ICF_best=max(inppowerout)/inppowerout_best;

resolution=1/findFWHM(idlpl_tilt,inppowerout/max(inppowerout));
% 
zzexp_arg=(idlpl_tilt)'.*1e5;
   zzexp_func=(inppowerout)';
toc(t_start)

%idlpl_tilt=idlpl_tilt+4e-5;
%%

figresol=figure(109);
hold all
        clf
        set(figresol,'name','resolution without slit','numbertitle','off');
        plot(idlpl_tilt,inppowerout(1:index)./inppowerout_1to1,'linewidth',2);
        %plot(idlpl_tilt,ans,'linewidth',2);
        xlabel('\Delta\lambda/\lambda');
        ylabel('ICF');
        try
            xlim([idlpl_tilt(1) idlpl_tilt(end)])
            xlim([-1.5e-4 1.5e-4])
        end
% text(1,1,sprintf('Resolution = %.0f \nICF_{1to1} = %.3f \nICF_{best} = %.3f ', resolution, ICF, ICF_best),...
%                     'HorizontalAlignment','right','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized');
text(1,1,sprintf('Resolving power = %.0f ', resolution),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
text(0,1,sprintf(' E=%.0f eV', round(energy/10)*10),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');

%  zzexp_arg=slitintscale.';
%  zzexp_func=slitint;
set(figure(109), 'Position', [100, 100, 550, 450]);
disp(' ');
disp(['done, resolution=',num2str(resolution)]);

% break
%%
% save([nm_icf_mat,'rough.mat'],'xlamds', 'idlpl_tilt','inppowerout','inppowerout_1to1', 'transmission')
% 
% %%
% 
% load(nm_icf_mat);
%%


coeff=3;
%load('C:\-D-\Work\LCLS\New_matlab\Steadystate\500_8_new.mat'); % delete!!!!!!!!!
k0=2*pi./(xlamds*(1+idlpl_tilt));
k0l=2*pi./(xlamds*(1+coeff*idlpl_tilt));
k=linspace(k0l(1),k0l(end),coeff*numel(k0));
I_function=inppowerout(1:index)./inppowerout_1to1.*transmission;

% % % % % I_function=I_function-max(I_function)*0.02;
% % % % % I_function(I_function<0)=0;
% % % % % I_function=I_function+max(I_function)*0.02;
% % % % 
% % % % %I_function=I_function/max(I_function)*0.05*ICF;
 I_function=interp1([k0l(1) k0 k0l(end)],[I_function(1)/5 I_function I_function(end)/5],k,'cubic');

 [Phase,haxes]=KKphase1(k,I_function);

%[Phase,haxes]=KKphase1(k0,I_function);
xlim(haxes(1),[idlpl_tilt(1) idlpl_tilt(end)]);
xlim(haxes(2),[idlpl_tilt(1) idlpl_tilt(end)]);
ylim(haxes(1),[0 0.17]);
set(haxes(2),'YTick',-3:1:3)
set(haxes(1),'YTick',0:0.02:0.17)
set(figure(1001), 'Position', [100, 100, 550, 450]);

cd('c:\Users\Svitozar\GoogleDrive\Documents\Dissertation\Figures\XFEL');
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