function array=readfield(nm, M)

% nm='LCLS_1.out.dfl'; %filename
% M=301; %number of transverse points
fd=fopen(nm,'r');

if fd==-1
         error('no such file');
end

size=dir(nm);
nslice=size.bytes/M/M/16;
array=single(zeros(M,M,nslice)+1i*zeros(M,M,nslice));
b=zeros(M,M)+1i*zeros(M,M);

for iii=1:nslice             
        x=fread(fd,2*M*M,'double');
        b=x(1:2:end)+1i*x(2:2:end);
        x0=reshape(b,M,M);
        array(:,:,iii)=single(x0);
        clear x b
end

fclose(fd);