%    nm_p='C:\-D-\Work\LCLS\tmp\930_u2_tdp_40.out';
%    readoutput=1;
 
% readm=0 - do not read
% readm=1 - force read 
% readm=2 - auto

 
function [d]=outread(nm_p,readoutput,readdfl,readm)
scriptversion='2.7';

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
        disp('  - .mat file found');
        file_mat = dir([nm_p,'.mat']);
        date_mat=datenum(file_mat.date);
    end
    
    if ~exist([nm_p,'.dfl'],'file')
        date_dfl=0;
    else
        disp('  - .dfl file found');
        file_dfl = dir([nm_p,'.dfl']);
        date_dfl=datenum(file_dfl.date);
    end
   
if (exist([nm_p,'.mat'],'file'))
    try
        matfile_parm=load ([nm_p,'.mat'],'parm');
        matfile_vervion=matfile_parm.parm.scriptversion;
        matfile_dflpresent=matfile_parm.parm.dflpresent;
    catch
        matfile_vervion='0'; %no file
        matfile_dflpresent='0'; %no file
    end
else
    matfile_vervion='0'; %no file
end

if (readm~=0  && 1 && exist([nm_p,'.mat'],'file') && date_mat>date_out &&  strcmp(scriptversion,matfile_vervion) && ~(readdfl&&~matfile_dflpresent)) || readm==1
%     file_mat = dir([nm_p,'.mat']);
%     file_out = dir([nm_p,'.out']);
    %fprintf('\r');
    if date_out==0
        disp(' - .out file not found');
        if ~strcmp(scriptversion,matfile_vervion)
            warntext='    ! WARNING\n    .mat file was recorded with a different version of outread script.\n    No file for reanalysis. Problems may appear.\n';
            fprintf(warntext);
        end
    end
    disp(' -Loading from present .mat file:');
    disp(['  -',nm_p,'.mat']);
    d=load([nm_p,'.mat'],'d');
    d=d.d;
    

%     if any(strcmp('scriptversion',fieldnames(d.parm)))
%         if ~strcmp(d.parm.scriptversion,scriptversion);
%             disp(warntext);
%             %disp(warntext);
%         end
%     else
%         disp(warntext);
%     end
    
    disp('  +done');
%    disp(' ');
else

    if ~strcmp(scriptversion,matfile_vervion)
            warntext='    ! .mat file was recorded with a different version of outread script.\n    Recalculating.\n';
            fprintf(warntext);
    end
    
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
        if inp.nslice~=0
            outp.power.E=outp.power.mean_S*inp.xlamds*inp.zsep*inp.nslice/c;
        else
            outp.power.E=outp.power.mean_S*inp.xlamds*inp.zsep*outp.Sn/c;
        end
        
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
            fprintf('\n   +done in %.1f seconds \r',t);
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
    %%
    
%     if readdfl && date_dfl>date_out
%         disp(' -import of .dfl field file,');
%         disp(['  generated ',num2str(round((date_dfl-date_out)*24*60*60)),' sec after .out']);
%         tic
%         [Xt,dfl.nslice]=fieldimport_all([nm_p,'.dfl'],inp.ncar,1);
%         t=toc;
% %       disp('  +done');
%         dfl.Sscale=linspace(0,inp.xlamds*inp.zsep*dfl.nslice,dfl.nslice);
%         dfl.E=sum(sum(sum(abs(Xt).^2)))/dfl.nslice*abs(dfl.Sscale(end)-dfl.Sscale(1))/3e8;
%         dfl.power=reshape(sum(sum(abs(Xt).^2,1),2),1,[]);
%         fprintf('   +done in %.1f seconds \r',t);
% 
%         disp(' -calculation of .dfl field transverse distribution');        
%         
%         Xfz=fftshift(fft2(ifftshift(Xt)));
%         dfl.Xfz=sqrt(sum(abs(Xfz).^2,3)).*exp(1i.*angle(mean(Xfz,3)));
%         clear Xfz
%         dfl.X=sqrt(sum(abs(Xt).^2,3)).*exp(1i.*angle(mean(Xt,3))); %2D field
%         disp('   +done');
%         
%         disp(' -calculation of .dfl field file spectrum');
%         tic
%         [Xf,dfl.Lamdscale]=dfl_time2freq(Xt,dfl.Sscale,inp.xlamds);
%         t=toc;
%         fprintf('   +done in %.1f seconds \r',t);
%         clear Xt
%         dfl.spectrum=reshape(sum(sum(abs(Xf).^2,1),2),1,[]);
%         clear Xf
%         
%         dflpresent=1;
%     else
%         dflpresent=0;
%     end
    
if readdfl %&& date_dfl>date_out
        disp(' -import of .dfl field file,');
        disp(['  generated ',num2str(round((date_dfl-date_out)*24*60*60)),' sec after .out']);
        tic
        [Xt,dfl.nslice]=fieldimport_all([nm_p,'.dfl'],inp.ncar,1);
        t=toc;
%       disp('  +done');
        dfl.Sscale=linspace(0,inp.xlamds*inp.zsep*dfl.nslice,dfl.nslice);
        dfl.E=sum(sum(sum(abs(Xt).^2)))/dfl.nslice*abs(dfl.Sscale(end)-dfl.Sscale(1))/3e8;
        dfl.power=reshape(sum(sum(abs(Xt).^2,1),2),1,[]);
        [Xf,dfl.Lamdscale]=dfl_time2freq(Xt,dfl.Sscale,inp.xlamds);
        clear Xt
        dfl.spectrum=reshape(sum(sum(abs(Xf).^2,1),2),1,[]);
                disp(' -backpropagation of .dfl'); 
        Xf=prop_TF(Xf,inp.leng,inp.xlamds,-1.1);

        disp(' -calculation of .dfl field transverse distribution');        
        
        Xfz=fftshift(fft2(ifftshift(Xf)));

        dfl.X=sqrt(sum(abs(Xf).^2,3)).*exp(1i.*angle(mean(Xf,3))); %2D field
        clear Xf
        dfl.Xfz=sqrt(sum(abs(Xfz).^2,3)).*exp(1i.*angle(mean(Xfz,3)));
        clear Xfz       
        disp('   +done');
        dflpresent=1;
else
        dflpresent=0;
end
    
    %% Data cataloging

    d.inp=inp;
    d.outp=outp;
    d.nm_p=nm_p;
    parm.readoutput=readoutput;
    parm.scriptversion=scriptversion;
    parm.dflpresent=dflpresent;
    d.parm=parm;
    if readdfl && date_dfl>date_out
        d.dfl=dfl;
    end
    %clear inp outp
    
    fprintf(' -saving data to .mat file\r');
    save([nm_p,'.mat'],'d');
%     save([nm_p,'.mat'],'scriptversion','-append');
    save([nm_p,'.mat'],'parm','-append');
    fprintf('   +done\r');
    
    fprintf(' +end of import');
  %  disp(' ');

end

return
