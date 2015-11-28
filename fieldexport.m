%imports field and pafameters from file.out.dfl and file.out
%correspondingly

function fieldexport(X_in,filename)
    
%     fd=fopen(filename,'r');
%     
%     if fd==-1
%              error('no such file');
%     end
%     
%     x=fread(fd,2*M*M,'double');
%     b=x(1:2:end)+1i*x(2:2:end);
%     if singlepromt
%         X=single(reshape(b,M,M));
%     else
%         X=reshape(b,M,M);
%     end
%     
%     fclose(fd);
X_in=permute(X_in,[2 1 3]);
[Mx,My,Mz]=size(X_in);
if Mx~=My
    error('export field has different number of points along transverse dimentions');
end

X_row=reshape(X_in,1,Mx*Mx*Mz);
clear X_in
X_out=zeros(1,2*Mx*Mx*Mz,'single');
X_out(1:2:end)=real(X_row);
X_out(2:2:end)=imag(X_row);
clear X_row
fd=fopen(filename,'wb');
fwrite(fd,X_out,'double');
fclose(fd);