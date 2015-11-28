%imports field and pafameters from file.out.dfl and file.out
%correspondingly

function X=fieldimport_tdp(filename,M,N,singlepromt)
    
    fd=fopen(filename,'r');
    
%      FileInfo = dir([pwd,'\',filename]);
%      Nslice=FileInfo.bytes/(2*M*M*8);
    
    if fd==-1
             error('no such file');
    end
    
    if singlepromt
        x=single(fread(fd,M*M*N*2*8,'double'));
    else
        x=fread(fd,M*M*N*2*8,'double');
    end
    b=x(1:2:end)+1i*x(2:2:end);
    b=b(1:M*M*N,:);
    size(b)
    X=reshape(b,M,M,N);
    %X=reshape(b,M,M);
        
    fclose(fd);