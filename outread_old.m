% reads .out file
clear all
fclose all
nm_p='tdp\700_u2_tdp_u7_DA.out';
fdin=fopen(nm_p,'r');
%% verifying beginning of the file
for i=1:8
    tline = fgetl(fdin);
    if strcmp(sscanf(tline,'%s'),'$newrun') || strcmp(sscanf(tline,'%s'),'$NEWRUN') 
        disp(' -starting import of input parameters');
        break
    else
        tline = fgetl(fdin);
        %error('start string "$newrun" not found');
    end
end
%% reading input parameters
for i=1:200
    tline = fgetl(fdin);
    if strcmp(sscanf(tline,'%s'),'$end') || strcmp(sscanf(tline,'%s'),'$END')
        disp('  +done');
        clear val var i k
        break
    end

    done=0;
    if done~=1
        k=strfind(tline,'D');
            if ~isempty(k)
                    tline(k)='E';
                    done=1;
            end
    end
    clear done k
    k=strfind(tline,'='); 
    var=genvarname(sscanf(tline(1:k-1),'%c'));
    val=single(sscanf(tline(k+1:end),'%f'));
    if isempty(val)
        val=sscanf(tline(k+1:end),'%s');
    end
    eval(['inp.' var '=val;']);
end
%% reading var names:  z[m] aw qfld
for i=1:20
    %tline = fgetl(fdin);
    if isempty(strfind(tline,'z[m]'));
        tline = fgetl(fdin);
    else
        [~,countvars]=sscanf(tline,'%s');
        vars=cell(1,countvars);
        vars{1}=sscanf(tline,'%s',1);
        for m=2:countvars
            k=strfind(tline,vars{m-1});
            vars{m}=sscanf(tline(k+size(vars{m-1},2):end),'%s',1);
        end
        disp(' -starting import of magnetic parameters');
        clear k i m
        vars{1}='z';
        break
    end
    
    if i==20; error('Import error'); end
end

%% reading var values:  z[m] aw qfld

tline = fgetl(fdin);
for i=1: inp.zstop/inp.xlamd+3
    if isempty(sscanf(tline,'%f'));
        %disp(i);
        disp('  +done'); 
        clear i
        break 
    end
    vals(i,:)=single(sscanf(tline,'%f'));
    
    tline = fgetl(fdin);
end    
    
for i=1:countvars
eval(['outp.' vars{i} '=vals(:,' num2str(i) ')'';']);
end
outp.Zn=size(vals,1);
clear vals i countvars
    
%% reading rest var names

for i=1:6 %find initial line
    if isempty(strfind(tline,'slice'));
        tline = fgetl(fdin);
    else
        break
        %error('start string "$newrun" not found');
    end
end

slicen(1)=single(sscanf(tline(strfind(tline,'slice')+7:end),'%f'));
fgetl(fdin); tline = fgetl(fdin);
current(1)=single(sscanf(tline,'%f',1));
fgetl(fdin); fgetl(fdin); tline = fgetl(fdin);

k=strfind(tline,'-'); tline(k)='_';
k=strfind(tline,'<'); tline(k)='';
k=strfind(tline,'>'); tline(k)='';

C=textscan(tline,'%s');
vars=C{1}';
clear C

tline = fgetl(fdin);

%% reading output values for first slice 

vals=single(zeros(outp.Zn,size(vars,2)));  
for m=1:outp.Zn
vals(m,:)=single(sscanf(tline,'%f'));
tline = fgetl(fdin);
end

for i=1:size(vars,2);
eval(['outp.' vars{i} '=vals(:,' num2str(i) ')'';']);
end

%% reading output values for 2:end slices 

for I=2:200%inp.nslice
    
    tline = fgetl(fdin);
    slicen(I)=single(sscanf(tline(strfind(tline,'slice')+7:end),'%f'));
    fgetl(fdin); tline = fgetl(fdin);
    current(I)=single(sscanf(tline,'%f',1));
    fgetl(fdin); fgetl(fdin); fgetl(fdin); tline = fgetl(fdin);
    
    %vals=single(zeros(outp.Zn,size(vars,2)));  
    for m=1:outp.Zn
    vals(m,:)=single(sscanf(tline,'%f'));
    tline = fgetl(fdin);
    end

    for i=1:size(vars,2);
    eval(['outp.' vars{i} '(' num2str(I) ',:)=vals(:,' num2str(i) ')'';']);
    end
    
end

%     if    isempty(strfind(tline,'z[m]'));
%         tline = fgetl(fdin);
%     else
%         [~,countvars]=sscanf(tline,'%s');
%         vars=cell(1,countvars);
%         vars{1}=sscanf(tline,'%s',1);
%         for m=2:countvars
%             k=strfind(tline,vars{m-1});
%             vars{m}=sscanf(tline(k+size(vars{m-1},2):end),'%s',1);
%         end
%         disp(' -starting import of magnetic parameters');
%         break
%     end
%     
%     if i==20; error('Import error'); end



break
%%
if ~ischar(tline), break, end
        
if findstr(tline,'$newrun')
    k=findstr(tline,'='); 
    tline(k)=' ';
    k=findstr('D',tline);
        if ~isempty(k)
                tline(k)='E';
        end
    npart=sscanf(tline,'%*s%f');
end

