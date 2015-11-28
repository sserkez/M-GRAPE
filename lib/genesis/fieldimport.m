%imports field and pafameters from file.out.dfl and file.out
%correspondingly

function X=fieldimport(filename,M,singlepromt)
    
    fd=fopen(filename,'r');
    
%      FileInfo = dir([pwd,'\',filename]);
%      Nslice=FileInfo.bytes/(2*M*M*8);
    
    if fd==-1
             error('no such file');
    end
    
    if singlepromt
        x=single(fread(fd,M*M*2,'double'));
    else
        x=fread(fd,M*M*2,'double');
    end
    
    b=x(1:2:end)+1i*x(2:2:end);
    
    %X=reshape(b,M,M,Nslice);
    X=reshape(b,M,M);
        
    fclose(fd);