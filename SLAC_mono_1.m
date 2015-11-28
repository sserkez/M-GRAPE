%Mono

clear all;
%close all;
fclose all;

% resolution without slit calculation
t_start=tic;
showpictures=1;
realfield=0;
roughnessincluded=0;
aberrationsincluded=0;
disp('---------------------------------------------------');
disp(['realfield =',num2str(realfield)]);
disp(['roughness =',num2str(roughnessincluded)]);
disp(['abberation=',num2str(aberrationsincluded)]);

interpM=3;%12
interpL=2;%5



disp(' ');
disp(['interpN   =',num2str(interpM)]);
disp(['interpleng=',num2str(interpL)]);

nm_f='700.out.dfl';
nm_p='700.out';

if realfield
    inread;
    K=2*pi/xlamds;
else
    energy=700;
    M=301;
    leng=10e-4;

    xlamds=1239.8/energy/1e9;
    K=2*pi/xlamds;
end
   
P0=1000; %inout power for model gaussian
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
Teta_g_i=1.04*0.01745;     %incidence angle new according to PRD   [rad]


R_m_tang=23.2;          %tangential radius of mirror M2 [m]
Teta_m_i=0.015;         %M2 incidence angle original [rad]
%Teta_m_i=0.016;         %M2 incidence angle original
z3_t_offset=1.53;       %Grating to refocusing mirror distance original
%z3_t_offset=1.52;       %Grating to refocusing mirror distance 
slitpos=1.350;          %grating to slit distance

U1_to_G=1.265;      %First undulator to grating [m]
G_to_U2=3.353;       %Grating to second undulator [m]

gafoc;

resolshift=1e-3; %scanned resolution shift via tilt on the grating (=dlambda/lambda)
dTeta_g_d=xlamds*resolshift/Teta_g_d/gr_sigma;
%%
leng_i=leng*interpL;

if realfield
    X=fieldimport(nm_f,M,1);
    [H{1}]=fieldplot(1,X,leng,'initial real field');
    X=fieldinterpolate(X,leng,1,interpM,interpL,'linear');
    [H{2}]=fieldplot(2,X,leng_i,'interpolated initial real field');
    X=prop_TF(X,leng_i,xlamds,U1_to_G);
    [H{3}]=fieldplot(3,X,leng_i,'real field at grating');
else
    leng_i=leng_i;
    X=fieldgaussian(M*interpM,leng_i,s1,s1,-z1,-z1,xlamds,P0);
    [H{4}]=fieldplot(4,X,leng_i,'gaussian field before grating');
end

%X=aperture(X,leng_i,g_length*Teta_g_i,'y');
X=grating(X,leng_i,K,z1_prime,slitpos,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d,dTeta_g_d,D0,D1,D2,aberrationsincluded);
X=aperture(X,leng_i,g_length*Teta_g_d,'y');
%X=lens(X,leng_i,K,z1_prime,slitpos,f_g_s,f_g_t,R_g_tang,R_g_sag,Teta_g_i,Teta_g_d,D2*1e9,abberationincluded);
%[H{5}]=fieldplot(5,X,leng_i,'field after grating');
[H{6}]=fieldplot(6,X,leng_i,'field after grating');

X_noslit=prop_TF(X, leng_i, xlamds,z3_t_offset);
%X=prop_TF(X, leng_i, xlamds,0.04);
X=prop_TF(X, leng_i, xlamds,slitpos);
[H{7}]=fieldplot(7,X,leng_i,'field at slit');

%%
% slitinterp=1/5;
% [H{8}]=fieldplot(8,X,leng_i,'field bi');
% Xi=fieldinterpolate(X,leng_i,1,slitinterp,slitinterp,'spline');
% [H{9}]=fieldplot(9,Xi,leng_i*slitinterp,'field ai');
% [z_opt,fwhm_min,fwhm_0]=waistscan(Xi,leng_i*slitinterp,xlamds,-0.02:0.001:0,2);
% fprintf('optimal waist is %.0f%% smaller\n',(fwhm_0-fwhm_min)/fwhm_0*100);
%%
%X=prop_TF(X, leng_i, xlamds,z3_t_offset-slitpos);
%X=mirror(X,leng_i,K,z3_t_offset-z2_t,z4_t,R_m_tang,R_m_tang,Teta_m_i,'y',aberrationincluded);
X=mirror(X_noslit,leng_i,K,z3_t_offset-z2_t,z4_t,R_m_tang,R_m_tang,Teta_m_i,'y',aberrationsincluded);
[H{10}]=fieldplot(10,X,leng_i,'field at M2');
X=prop_TF(X, leng_i, xlamds,G_to_U2-z3_t_offset);
[H{11}]=fieldplot(11,X,leng_i,'field at U2');

% [z_opt,fwhm_min,fwhm_0]=waistscan(X,leng_i,xlamds,-1:0.2:3,2);
% fprintf('optimal waist is %.0f%% smaller\n',(fwhm_0-fwhm_min)/fwhm_0*100);
%%
X=fieldinterpolate(X,leng_i,0,301,leng,'spline');
X=X*sqrt(20/sum(sum(abs(X).^2)));
[H{12}]=fieldplot(12,X,leng,'interpolated at U2');

     Pi=sum(sum(abs(X).^2));
     X=X*sqrt(20/sum(sum(abs(X).^2)));

fieldexport(X,'prop_field.dfl')

dos('genesis301');
    leng_u2=4e-4;
    Xn=fieldimport('700_prop.out.dfl',M,1);
    [H{13}]=fieldplot(13,Xn,leng_u2,'amplified field');
