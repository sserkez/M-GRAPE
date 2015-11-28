function[xout]=aperture_a(xin,leng_x,leng_y,slitsize,dim)

[My,Mx]=size(xin);
dx=leng_x/Mx;
dy=leng_y/My;
[xx,yy]=meshgrid(((Mx-1)/2+1-(1:Mx)).*dx,((My-1)/2+1-(1:My)).*dy);
%[xx,~]=meshgrid(-leng/2+dx/2:dx:leng/2-dx/2,-leng/2+dx/2:dx:leng/2-dx/2);
xx=abs(xx)<=slitsize/2;
yy=abs(yy)<=slitsize/2;
if dim=='x'
    xout=xin.*xx;
elseif dim=='y'
    xout=xin.*yy';
else
    error('argument dim must be a string value of ''x'' or ''y''');
end
