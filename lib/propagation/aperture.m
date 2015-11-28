function[xout]=aperture(xin,leng,slitsize,dim)

[M,~]=size(xin);
dx=leng/M;
[xx,~]=meshgrid(((M-1)/2+1-(1:M)).*dx);
%[xx,~]=meshgrid(-leng/2+dx/2:dx:leng/2-dx/2,-leng/2+dx/2:dx:leng/2-dx/2);
xx=abs(xx)<=slitsize/2;
if dim=='x'
    xout=xin.*xx;
elseif dim=='y'
    xout=xin.*xx';
else
    error('argument dim must be a string value of ''x'' or ''y''');
end
