%read data from *.in file
%fdin=fopen([nm,'.in'],'r');
%function fieldparmsimport(filename)
fdin=fopen(nm_p,'r');

for lines=1:200
        tline = fgetl(fdin);
        if ~ischar(tline), break, end
        
        if findstr(tline,'npart')
            k=findstr(tline,'='); 
            tline(k)=' ';
            k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            npart=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'ncar')
            k=findstr(tline,'='); 
            tline(k)=' ';
            k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            M=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'rxbeam')
            k=findstr(tline,'='); 
            tline(k)=' ';
            k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            rxbeam=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'zrayl')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            zrayl=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'rybeam')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            rybeam=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'xlamds')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            xlamds=sscanf(tline,'%*s%f');
        end
        
        
        if findstr(tline,'rmax0 ')
            k=findstr(tline,'=');
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            rmax0=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'dgrid')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            dgrid=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'nslice')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            nslice=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'zsep')
            k=findstr(tline,'='); 
            tline(k)=' ';
                k=findstr('D',tline);
                if ~isempty(k)
                        tline(k)='E';
                end
            zsep=sscanf(tline,'%*s%f');
        end
        
        if findstr(tline,'$end')
            break
        end
end

     
     if dgrid==0
        rbeam=sqrt(rxbeam^2+rybeam^2);
        ray=sqrt(zrayl*xlamds/pi);
        leng=rmax0*(rbeam+ray); %%%! *2
     else
         leng=dgrid*2;
     end
     
     clear rxbeam rybeam rbeam ray rmax0
     
     dx=leng/M;     dy=dx;
     K=2*pi/xlamds;

energy=1239.8/xlamds/1e9;
fclose(fdin);
clear fdin k tline lines