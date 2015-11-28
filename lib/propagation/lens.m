%gr_asym - grating asymmetry parameter = Teta_g_i/Teta_g_d

function[xout]=lens(xin,leng,xlamds,f)
K=2*pi/xlamds;
if numel(size(xin))==1
[M,~]=size(xin);
dx=leng/M;
[xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);
%C_20=(Teta_i^2/z1+Teta_i/R_t+Teta_d^2/z2+Teta_d/R_t)/2;

 xout=xin.*exp(-1i*K/f/2.*(yy.^2+xx.^2));          %focusing
else
    
[M,~,N]=size(xin);
dx=leng/M;
[xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);
%C_20=(Teta_i^2/z1+Teta_i/R_t+Teta_d^2/z2+Teta_d/R_t)/2;
xout=xin;
for i=1:N
 xout(:,:,i)=xin(:,:,i).*exp(-1i*K/f/2.*(yy.^2+xx.^2));          %focusing
end
end
 
 
    
%     x=x.*exp(1i*K*C_30*(yy.^3)/Teta_i^3);       %abberation 30
%     xout=x.*exp(1i*K*C_12*xx.^2.*yy/Teta_i);    %abberation 12



% disp(['1/f=',num2str(1/f/2)]);
% disp(['C_20=',num2str(C_20/Teta_i^2)]);
