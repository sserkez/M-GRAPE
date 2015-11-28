function[X_out]=mirror(X_in,leng,xlamds,z1,z2,R_t,R_s,Teta_i,dim,abberationsincluded)
K=2*pi/xlamds;
Teta_d=Teta_i;
f=R_t/(1/Teta_d+Teta_i/Teta_d^2);

[M,~]=size(X_in);
dx=leng/M;
[xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);

if abberationsincluded
%C_20=(Teta_i^2/z1+Teta_i/R_t+Teta_d^2/z2+Teta_d/R_t)/2;
C_30=(sin(Teta_i)^2/z1-sin(Teta_i)/R_t)*cos(Teta_i)/2/z1-(sin(Teta_d)^2/z2-sin(Teta_d)/R_t)*cos(Teta_d)/2/z2;
%C_30=(Teta_i^2/z1-Teta_i/R_t)/2/z1-(Teta_d^2/z2-Teta_d/R_t)/2/z2;
C_12=(-Teta_i/R_s/z1+1/z1^2-Teta_d/R_s/z2+1/z2^2)/2;
else
 C_30=0;
 C_12=0;
end
 
if dim=='x'
    X=X_in.*exp(-1i*K/f/2.*(xx.^2));             %focusing
    X=X.*exp(1i*K*C_30*(xx.^3)/Teta_i^3);       %abberation 30
    X_out=X.*exp(1i*K*C_12*yy.^2.*xx/Teta_i);    %abberation 12
elseif dim=='y'
    X=X_in.*exp(-1i*K/f/2.*(yy.^2));             %focusing
    X=X.*exp(1i*K*C_30*(yy.^3)/Teta_i^3);       %abberation 30
    X_out=X.*exp(1i*K*C_12*xx.^2.*yy/Teta_i);    %abberation 12
else
    error('argument dim must be a string value of ''x'' or ''y''');
end

% disp(['1/f=',num2str(1/f/2)]);
% disp(['C_20=',num2str(C_20/Teta_i^2)]);
