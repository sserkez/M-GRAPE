%    nm_p='C:\-D-\Work\LCLS\tmp\930_u2_tdp_40.out';
%    readoutput=1;
 

 
 
function [d]=outread(nm_p,readoutput)
scriptversion='1.3';

    fprintf('\r');
    disp(['---',nm_p,'---']);

if ~exist([nm_p,'.mat'],'file') && ~exist(nm_p,'file')
    error('Both .out and .out.mat files not found');
end



    
    if ~exist(nm_p,'file')
        date_out=0;
    else
        file_out = dir(nm_p);
        date_out=datenum(file_out.date);
    end
    
    if ~exist([nm_p,'.mat'],'file')
        date_mat=0;
    else
        file_mat = dir([nm_p,'.mat']);
        date_mat=datenum(file_mat.date);
    end
    

if (exist([nm_p,'.mat'],'file') && date_mat>date_out && 1)
    file_mat = dir([nm_p,'.mat']);
    file_out = dir([nm_p,'.out']);
    %fprintf('\r');
    if date_out==0
        disp(' - .out file not found');
    end
    disp(' -Loading from present .mat file:');
    disp(['  -',nm_p,'.mat']);
    load ([nm_p,'.mat']);
    
    warntext='    !Warning, .mat file was recorded with an older version of outread script; consider reassembly';
    if any(strcmp('scriptversion',fieldnames(d.parm)))
        if ~strcmp(d.parm.scriptversion,scriptversion);
            disp(warntext);
            %disp(warntext);
        end
    else
        disp(warntext);
    end
    
    disp('  +done');
    disp(' ');
else

    intN=2;
    c=299792458;

    disp(' -reading .out file');
    tic;
    fdin=fopen(nm_p,'r');
    if fdin==-1
        error('file not found');
    end
    A = textscan(fdin,'%s','delimiter','\n');
    A=A{1,1};
    t=toc;
    %disp('   +done');
    fprintf('   +done in %.1f seconds \r',t);
    fclose all; clear ans fdin;
    %% Input parameters import
    disp(' -import of input parameters');
    is=find(strcmp(A(1:10),'$newrun'))+1;
    ie=find(strcmp(A(1:200),'$end'))-1;

    for i=is:ie
    k=strfind(A{i},'D');
    A{i}(k)='E';

    k=strfind(A{i},'='); 
    var=genvarname(lower(sscanf(A{i}(1:k-1),'%c')));
    val=single(sscanf(A{i}(k+1:end),'%f'));
    if isempty(val)
        val=sscanf(A{i}(k+1:end),'%s');
    end
    eval(['inp.' var '=val;']);
    end
    outp.Zn=round(inp.zstop/inp.xlamd/inp.iphsty+1);
    outp.Sn=round((numel(A)-(181+outp.Zn+2+7))/(outp.Zn+8))/inp.ishsty;
    outp.Zn=single(sscanf(A{ie+4},'%f'));
    outp.Sn=single(sscanf(A{ie+5},'%f'));

    if inp.dgrid==0
        rbeam=sqrt(inp.rxbeam^2+inp.rybeam^2);
        ray=sqrt(inp.zrayl*inp.xlamds/pi*(1+(inp.zwaist/inp.zrayl)^2)); %not cross-checked
        inp.leng=inp.rmax0*(rbeam+ray);
    else
        inp.leng=inp.dgrid*2;
    end
    
    disp('   +done');
    clear i is ie k val var

    %% Magnetic parameters import
    disp(' -import of magnetic parameters');
    vars={'z','aw','qfld'};
    vals=zeros(outp.Zn,3,'single'); %reserving memory

    is=find(strcmp(A(1:200),'z[m]          aw            qfld      '));


    for i=1:outp.Zn%reformatting text to numbers
    vals(i,:)=single(sscanf(A{is+i},'%f'));
    end
    for i=1:size(vars,2);%assigning to variables
    eval(['outp.' vars{i} '=vals(:,' num2str(i) ')'';']);
    end
    outp.Zscale=outp.z;
    clear outp.z
    disp('   +done');

    % for i=1:5
    %     if strcmp(sscanf(A{182+outp.Zn+i},'%17c',1),'=================');
    %         break
    %     end
    % end
    % 
    % 181+outp.Zn+2
    % 181+outp.Zn+2+7 початок +(K-1)*(outp.Zn+8)
    % 192+2*outp.Zn+(K-1)*(outp.Zn+8) кінець

    if readoutput
        %% Output parameters input
        disp(' -import of output parameters');
        pause(0.01);
        is1=is+outp.Zn;
        text=A{is1+7};
        k=strfind(text,'-'); text(k)='_';
        k=strfind(text,'<'); text(k)='';
        k=strfind(text,'>'); text(k)='';
        C=textscan(text,'%s');
        vars=C{1}';
        inp.Vn=round(numel(vars));
        clear C k text is

        vals=zeros(outp.Zn,outp.Sn,inp.Vn,'single');
        tic 
        fprintf('   -converting ACSII: %5.1f%% ',(0)*100);
        %convert all
        if readoutput==1
            for I=1:outp.Sn %reformatting text to numbers
                 AA=A((1:outp.Zn)+is1+7+(I-1)*(outp.Zn+7));

                 for i=1:outp.Zn
                     %vals(i,I,:)=sscanf(A{i+is1+7+(I-1)*(outp.Zn+7)},'%f');
                     vals(i,I,:)=sscanf(AA{i},'%f');
                 end

                 fprintf(repmat('\b',1,7));
                 fprintf('%5.1f%% ',(I/outp.Sn)*100);         
            end
            fprintf('\r');
             %disp('     +done');
            t=toc;
            fprintf('   +done in %.1f seconds \r',t);
            for i=1:size(vars,2);%assigning to variables
            eval(['outp.' vars{i} '.v=vals(:,:,' num2str(i) ')'';']);
            eval(['outp.' vars{i} '.mean_S=mean(outp.' vars{i} '.v,1);']);
            eval(['outp.' vars{i} '.mean_Z=mean(outp.' vars{i} '.v,2);']);
            eval(['outp.' vars{i} '.max_S=max(outp.' vars{i} '.v,[],1);']);
            eval(['outp.' vars{i} '.max_Z=max(outp.' vars{i} '.v,[],2);']);
            end

        else
            for I=1:outp.Sn %reformatting text to numbers
                 AA=A((1:outp.Zn)+is1+7+(I-1)*(outp.Zn+7));

                 i=outp.Zn;
                     %vals(i,I,:)=sscanf(A{i+is1+7+(I-1)*(outp.Zn+7)},'%f');
                     vals(i,I,:)=sscanf(AA{i},'%f');

                 fprintf(repmat('\b',1,7));
                 fprintf('%5.1f%% ',(I/outp.Sn)*100);         
            end
            fprintf('\r');
            t=toc;
            fprintf('   +done in %.1f seconds \r',t);
            for i=1:size(vars,2);%assigning to variables
            eval(['outp.' vars{i} '.last=vals(:,:,' num2str(i) ')'';']);
    %         eval(['outp.' vars{i} '.mean_S=mean(outp.' vars{i} '.v,1);']);
    %         eval(['outp.' vars{i} '.mean_Z=mean(outp.' vars{i} '.v,2);']);
    %         eval(['outp.' vars{i} '.max_S=max(outp.' vars{i} '.v,[],1);']);
    %         eval(['outp.' vars{i} '.max_Z=max(outp.' vars{i} '.v,[],2);']);
            end
        end


        for I=1:outp.Sn
            %outp.current(I,:)=sscanf(A{185+outp.Zn+(I-1)*(outp.Zn+7)},'%f',1);
            outp.current(I,:)=sscanf(A{is1+4+(I-1)*(outp.Zn+7)},'%f',1);
        end
        inp.charge=sum(outp.current*inp.zsep*inp.xlamds/c);
        outp.Sscale=linspace(0,inp.xlamds*inp.zsep*outp.Sn,outp.Sn)';


        sp1=sum(outp.power.v,1);
        %sp2=sum(outp.power.v,2);
        outp.power.peakpos=sum(repmat(outp.Sscale,1,outp.Zn).*outp.power.v,1)./sp1;
        %outp.power.std=sqrt(sum(((repmat(outp.Sscale,1,outp.Zn)-repmat(outp.power.peakpos,outp.Sn,1)).^2.*outp.power.v./repmat(sp1,outp.Sn,1)),1));
        outp.power.std=sqrt(sum(((repmat(outp.Sscale,1,outp.Zn)-repmat(outp.power.peakpos,outp.Sn,1)).^2.*outp.power.v),1)./sp1);
        %outp.power.peakpos=sum(repmat(outp.Sscale,1,size(outp.power.v,2)).*outp.power.v./repmat(sp1,outp.Sn,1));


        %outp.power.mean_S_2s=sp1./(max(outp.Sscale)-min(outp.Sscale)).*outp.power.std.*2;
        outp.power.E=outp.power.mean_S*inp.xlamds*inp.zsep*inp.nslice/c;
        
        outp.power.weight=outp.power.v./repmat(sum(outp.power.v,1),outp.Sn,1);
        outp.r_size.mean_S_norm=sum(outp.r_size.v.*outp.power.weight,1);
        clear outp.z vals vars i I A AA is1 sp1
        %% Spectrum calculation
        if inp.itdp==1
            tic
            disp(' -analysis of spectrum');
            sc=linspace(-outp.Sn/2,outp.Sn/2,outp.Sn*intN);%original
            %sc=-(outp.Sn-1)/2:1:(outp.Sn-1)/2;

            k0=2*pi/inp.xlamds;
            dk=2*pi/(outp.Sn*inp.xlamds*inp.zsep);
            outp.Kscale=k0+dk*sc';
            outp.Lamdscale=2*pi./outp.Kscale;
            outp.Ln=outp.Sn*intN;
            
            fprintf('   -bandwidth fit: %5.1f%% ',(0)*100);
            outp.spectrum_mid.v=abs(fftshift(fft(sqrt(outp.p_mid.v).*exp(1i*outp.phi_mid.v),intN*outp.Sn,1),1)).^2/sqrt(outp.Sn)./(2*inp.leng/inp.ncar)^2./1e10; %.*exp(-i*outp.phi_mid.v)
            outp.spectrum_mid.mean_S=mean(outp.spectrum_mid.v,1);
            outp.spectrum_mid.mean_Z=mean(outp.spectrum_mid.v,2);
            outp.spectrum_mid.max_S=max(outp.spectrum_mid.v,[],1);
            outp.spectrum_mid.max_Z=max(outp.spectrum_mid.v,[],2);
            %outp.spectrum_mid.lamdpos_max=outp.Lamdscale(outp.spectrum_mid.v)
            outp.spectrum_mid.lamdpos=sum(repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2)).*outp.spectrum_mid.v,1)./sum(outp.spectrum_mid.v,1);
            outp.spectrum_mid.lamdpos2=sum(repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2)).*outp.spectrum_mid.v.^2,1)./sum(outp.spectrum_mid.v.^2,1);
            %outp.spectrum_mid.std_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1))./sum(outp.spectrum_mid.v,1);
            outp.spectrum_mid.std_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1)./sum(outp.spectrum_mid.v,1));
            outp.spectrum_mid.std_lamd1=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos2,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1)./sum(outp.spectrum_mid.v,1));
            outp.spectrum_mid.std_lamd2=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v.^2,2))-repmat(outp.spectrum_mid.lamdpos2,size(outp.spectrum_mid.v.^2,1),1)).^2.*outp.spectrum_mid.v.^2),1)./sum(outp.spectrum_mid.v.^2,1));

%             for i=1:outp.Zn
%     %           [outp.spectrum_mid.s_gfit(i),~]=gaussfit(outp.Lamdscale*1e9,outp.spectrum_mid.v(:,i),outp.spectrum_mid.std_lamd2(i)*1e9,outp.spectrum_mid.lamdpos(i)*1e9);
%                 try
%                     [outp.spectrum_mid.s_gfit(i),~]=gaussfit(outp.Lamdscale*1e9,outp.spectrum_mid.v(:,i)./max(outp.spectrum_mid.v(:,i)));
%                 catch
%                 end
%                 fprintf(repmat('\b',1,7));
%                 fprintf('%5.1f%% ',(i/outp.Zn)*100); 
%             end
%             outp.spectrum_mid.s_gfit=outp.spectrum_mid.s_gfit.*1e-9;
%              fprintf('\r');

            
            clear sc k0 dk
            t=toc;
            fprintf('   +done in %.1f seconds \r',t);
        end
        %% 

    %     nm_f=[nm_p,'.dfl'];
    %     fd=fopen(nm_f);
    % 
    %     if fd~=-1
    %         disp(' -output field found');
    %         disp(['   -importing "',nm_f,'"']);
    %         [XX,Nslice]=fieldimport_all(nm_f,inp.ncar,1);
    % 
    %         if Nslice~=1 && inp.itdp~=1
    %             error(['numerous slices in ',nm_f,' while calculation was done in steadystate, accordind to flag itdp~=1']);
    %         end
    % 
    %         if Nslice~=1        
    %             XX=XX(:,:,Nslice-outp.Sn+1:end);
    %             XXf=cat(3,zeros(inp.ncar,inp.ncar,outp.Sn*(intN-1)),XX);
    %             XXf=fftshift(fft(XXf,[],3),3);
    %             outp.spectrum_whole.v=reshape(sum(sum(abs(XXf).^2,1),2),1,[]);
    %             outp.field=XX;
    %         else
    %             outp.field=XX;
    %         end
    %         disp('   +done');
    %     end


    end
    %% Data cataloging

    d.inp=inp;
    d.outp=outp;
    d.nm_p=nm_p;
    d.parm.readoutput=readoutput;
    d.parm.scriptversion=scriptversion;
    %clear inp outp
    
    fprintf(' -saving data to .mat file\r');
    save([nm_p,'.mat'],'d');
    fprintf('   +done\r');
    
    fprintf(' +end of import');
    disp(' ');
    
end

return
