function dist=distread(filename)
% reads Genesis distribution files

%filename = 'C:\-D-\WORK\LCLS\tmp\3\ebeams\150pC1p5kA3p5GeV.dist';

delimiter = ' ';

formatSpec = '%f%f%f%f%f%f%*s%*s%[^\n\r]';

    disp(' ');
    disp(' -Reading distribution file from');
    disp(['  -',filename]);
    
    fd = fopen(filename,'r');
    if fd==-1
        error('file not found')
    end

for i=1:20
    M=fgetl(fd);
    if M(1)=='#'||M(1)=='?'
        disp(M)
        k=strfind(M,'=');
        if ~isempty(k)
            var=genvarname(lower(sscanf(M(2:k-1),'%c')));
            val=single(sscanf(M(k+1:end),'%f'));
            if isempty(val)
                val=sscanf(M(k+1:end),'%s');
            end
            eval(['dist.' var '=val;']);
        end
    else
        break
    end
end
clear M var val k

%break
frewind(fd)
dataArray = textscan(fd, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,i-1, 'ReturnOnError', false);

fclose(fd);


dist.X = dataArray{:, 1};
dist.PX = dataArray{:, 2};
dist.Y = dataArray{:, 3};
dist.PY = dataArray{:, 4};
dist.T = dataArray{:, 5};
dist.G = dataArray{:, 6};

if numel(dist.X)~=dist.size
    disp(' ');
    disp('  !!! Number of elements !!!');
    disp('  not consistent with header');
    disp('  ');
end

gavg=mean(dist.G);

dist.emitx=sqrt(mean(dist.X.^2).*mean(dist.PX.^2)-mean(dist.X.*dist.PX).^2).*gavg;
dist.emity=sqrt(mean(dist.Y.^2).*mean(dist.PY.^2)-mean(dist.Y.*dist.PY).^2).*gavg;
dist.betax=mean(dist.X.*dist.X).*gavg./dist.emitx;
dist.betay=mean(dist.Y.*dist.Y).*gavg./dist.emity;
dist.alphax=-mean(dist.X.*dist.PX).*gavg./dist.emitx;
dist.alphay=-mean(dist.Y.*dist.PY).*gavg./dist.emity;


    disp(['Beam energy = ',num2str(mean(dist.G)*3e8*3e8*9.1e-31/1.6e-19/1e9),' GeV'])
    disp('  -done');

clearvars filename delimiter startRow formatSpec fileID dataArray ans fd