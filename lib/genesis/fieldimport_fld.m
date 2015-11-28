%imports field and pafameters from file.out.dfl and file.out
%correspondingly

%filename=[nm_p{Di},'.dfl']; M=d(Di).inp.ncar; singlepromt=1;
function [X,Nslice]=fieldimport_all(filename,M,singlepromt)
    
    fd=fopen(filename,'r');
    
    if fd<0
             error('no such file');
    end
    
      FileInfo = dir([pwd,'\',filename]);
      FileInfo = dir(filename);
      %FileInfo.bytes
      Nslice=int32(FileInfo.bytes/(2*M*M*8));
      M=int32(M);
    
    
    if singlepromt
        x=single(fread(fd,M*M*2*Nslice,'double'));
    else
        x=fread(fd,M*M*2*Nslice,'double');
    end
    
    b=x(1:2:end)+1i*x(2:2:end);
    X=reshape(b,M,M,Nslice);
    X=permute(X,[2 1 3]);
    
    fclose(fd);
    Nslice=single(Nslice);