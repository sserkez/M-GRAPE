function[X_out]=grating_a(X_in,leng_x,leng_y,xlamds,z1,z2,R_t,R_s,Theta_i,Theta_d,dTheta_d,D0,D1,D2,aberration)
% X_in is the input complex field
% X_out is the output complex field
% leng is transverse size of the mesh (dgrid*2)
% xlamds is wavelength
% z1 and z2 are distances to source and image in order to calculate 
%   aberration effect (we use semi-analytical approach to propagate, 
%   so we need to provide this information)
% R_t is tangential radius of curvature
% R_s is sagittal radius of curvature
% Teta_i and Teta_d are incidence and reflection angles (obviously one can 
%   calculate it now internaly within function. I keep dummy input variable
%   to save time on revriting all versiong of script code.
% dTheta_d is angle deviation of the field from principle ray direction
%   after grating. It is used to introduce different angles of reflected
%   field from grating at different  wavelengths if time-dependent (3D) field
%   is propagated, resulting in pusle-front tilt.
% D0, D1, D2 are grating coefficients (in out case 1123[l/mm] 1.6[l/mm2] 0.002[l/mm3])
% aberration is logical value, that defines whether include aberration effects, or not

Theta_d=acos(cos(Theta_i)-xlamds*D0*1e3);
K=2*pi/xlamds;
gr_asym=Theta_i/Theta_d;  %grating asymmetry factor
gr_sigma=1/(D0*1e3);         %sigma parameter of grating (1/k)
gr_alpha=D1*1e6/(D0*1e3)^2;  %alpha parameter of grating (old=1,8) (n1/k/k) !!
D2=D2*1e9;

if D2==0
    n=0;
    D2=1;
else
    n=1;
end

% sagital focal distance calculation

    f_s=R_s/(sin(Theta_i)+sin(Theta_d));
    %disp(['f_g_s=',num2str(f_s),'m']);

% tangential focal distance calculation
    % vls contribution
    f_t_vls=gr_sigma^2*sin(Theta_d)^2/xlamds/gr_alpha; %!!! original!!!!
    %f_g_t_vls=gr_sigma^2/xlamds/gr_alpha;
    % curvature contribution
    f_t_c=R_t/(1/Theta_d+Theta_i/Theta_d^2); %tangential focal distance due to curvature

% final tangential focal distance
    f_t=1/(1/f_t_vls+1/f_t_c);
    %f_g_t=1.07;

%w1_prime=w1*(Teta_d/Teta_i);
%disp(['f_g_t=',num2str(f_t),'m,  f_g_t_c=',num2str(f_t_c),'m,  f_g_t_vls=',num2str(f_t_vls),'m']);
%disp(['z1''=',num2str(z1_prime),'m,  s1''=',num2str(w1_prime/2*1e6),'mikr,   zr''=',num2str(zr_prime),'m']);

[My,Mx]=size(X_in);
dx=leng_x/Mx;
dy=leng_y/My;
[xx,yy]=meshgrid(((Mx-1)/2+1-(1:Mx)).*dx,((My-1)/2+1-(1:My)).*dy);
%C_20=(Teta_i^2/z1+Teta_i/R_t+Teta_d^2/z2+Teta_d/R_t)/2;
if aberration
C_30=(n*pi*2/K/3*D2)+(sin(Theta_i)^2/z1-sin(Theta_i)/R_t)*cos(Theta_i)/2/z1-(sin(Theta_d)^2/z2-sin(Theta_d)/R_t)*cos(Theta_d)/2/z2;
C_12=(-Theta_i/R_s/z1+1/z1^2+Theta_d/R_s/z2-1/z2^2)/2;
else
C_30=0;
C_12=0;
end
% (n*pi*2/K/3*n2)
% (sin(Teta_i)^2/z1-sin(Teta_i)/R_t)*cos(Teta_i)/2/z1-(sin(Teta_d)^2/z2-sin(Teta_d)/R_t)*cos(Teta_d)/2/z2

% enegry conservation calculations
dMx=round(Mx*(1-gr_asym)/2);
P_in=sum(sum(abs(X_in(:,dMx:Mx-dMx)).^2));
%

X=X_in;

if size(X,1)~=1
    X=interp2(xx,yy,X,xx.*(gr_asym),yy,'linear');
    P_out=sum(sum(abs(X).^2));
else
    X=interp1(xx,X,xx.*(gr_asym),'linear');
    P_out=sum(sum(abs(X).^2));
end

if P_out~=0
    X=X*sqrt(P_in/P_out);
end

    X=X.*exp(-1i*K/f_t/2.*(xx.^2));          %focusing sag
    X=X.*exp(-1i*K/f_s/2.*(yy.^2));             %focusing tang
    X=X.*exp(1i*K*C_30*(xx.^3)/Theta_d^3);       %aberration 30 
    X=X.*exp(1i*K*C_12*yy.^2.*xx/Theta_d);    %aberration 12

    %X=abs(X).*exp(1i*(angle(X)+dTheta_d*K*xx)); %phase tilt
    X=X.*exp(1i*dTheta_d*K*xx);
    
    X_out=X;
    
%X_out=interp2(xx,yy,X_out,xx,yy.*(gr_asym),'linear');
    
%     x=x.*exp(1i*K*C_30*(yy.^3)/Teta_i^3);       %abberation 30
%     xout=x.*exp(1i*K*C_12*xx.^2.*yy/Teta_i);    %abberation 12

% disp(['1/f=',num2str(1/f/2)]);
% disp(['C_20=',num2str(C_20/Teta_i^2)]);
