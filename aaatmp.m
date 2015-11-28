leng=50e-6;
M=301;
fwhmint=10e-6;  % FWHM Int

X=fieldgaussian(M,leng,fwhmint/2.35,fwhmint/2.35,0,0,xlamds);
[H{7}]=fieldplot(7,X,leng,leng,'tryfield');
% sy=fwhmint/2.35;  % sigmaint
% sx=sy; zx=0; zy=0;
% dx=leng/M; dy=dx;
% wx=sx*sqrt(2);  % W ampl.
% wy=sy*sqrt(2);
% qy=(2j*pi*(wy^2)/xlamds)+zy; %?
% qx=(2j*pi*(wx^2)/xlamds)+zx;
% [xx,yy]=meshgrid((M-1)/2+1-(1:M));
% K=2*pi/xlamds;
% X=exp(-1j.*K./2.*((xx.*dx).^2./qx+(yy*dx).^2./qy)); %?
% %X=exp(-1.*((xx.*dx).^2+(yy*dx).^2)./wx^2);
% handle=fieldplot(5,X,leng,leng,'try');
% findFWHM(xx(1,:)*dx,abs(X(50,:)).^2)

%%
leng=100e-6;
M=100;
dx=leng/M;
[xx,yy]=meshgrid((M-1)/2+1-(1:M));
X=exp(-(xx*dx).^2/(50e-6)^2);
