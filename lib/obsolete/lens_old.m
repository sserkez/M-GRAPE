%gr_asym - grating asymmetry parameter = Teta_g_i/Teta_g_d

function[xout]=grating(xin,leng,K,z1,z2,f_t,f_s,R_t,R_s,Teta_i,Teta_d,n2,aberration)

if n2==0
    n=0;
    n2=1;
else
    n=1;
end

[M,~]=size(xin);
dx=leng/M;
[xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);
%C_20=(Teta_i^2/z1+Teta_i/R_t+Teta_d^2/z2+Teta_d/R_t)/2;
if aberration
C_30=(n*pi*2/K/3*n2)+(sin(Teta_i)^2/z1-sin(Teta_i)/R_t)*cos(Teta_i)/2/z1-(sin(Teta_d)^2/z2-sin(Teta_d)/R_t)*cos(Teta_d)/2/z2;
C_12=(-Teta_i/R_s/z1+1/z1^2+Teta_d/R_s/z2-1/z2^2)/2;
else
C_30=0;
C_12=0;
end

% (n*pi*2/K/3*n2)
% (sin(Teta_i)^2/z1-sin(Teta_i)/R_t)*cos(Teta_i)/2/z1-(sin(Teta_d)^2/z2-sin(Teta_d)/R_t)*cos(Teta_d)/2/z2

x=xin;
    x=x.*exp(-1i*K/f_s/2.*(yy.^2));          %focusing sag
    x=x.*exp(-1i*K/f_t/2.*(xx.^2));             %focusing tang
    x=x.*exp(1i*K*C_30*(xx.^3)/Teta_d^3);       %abberation 30 
    xout=x.*exp(1i*K*C_12*yy.^2.*xx/Teta_d);    %abberation 12
 
    
%     x=x.*exp(1i*K*C_30*(yy.^3)/Teta_i^3);       %abberation 30
%     xout=x.*exp(1i*K*C_12*xx.^2.*yy/Teta_i);    %abberation 12



% disp(['1/f=',num2str(1/f/2)]);
% disp(['C_20=',num2str(C_20/Teta_i^2)]);
