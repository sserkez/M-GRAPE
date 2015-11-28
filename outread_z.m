clear all
nm_p='C:\-D-\Work\SASE3_chicane\run4\U1.1.out';
 readoutput=1;
 Z=100;
 %Z=[];
 
 
%function [d]=outread(nm_p,readoutput)
intN=2;
c=299792458;

fprintf('\r');
disp(['---',nm_p,'---']);
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
var=genvarname(sscanf(A{i}(1:k-1),'%c'));
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

    if numel(Z)==0

    vals=zeros(outp.Zn,outp.Sn,inp.Vn,'single');
    tic 
    fprintf('   -converting ACSII: %5.1f%% ',(0)*100);
    %convert all

        for I=1:outp.Sn %reformatting text to numbers
             AA=A((1:outp.Zn)+is1+7+(I-1)*(outp.Zn+7));

             parfor i=1:outp.Zn
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
%         
    elseif numel(Z)==1
        vals=zeros(1,outp.Sn,inp.Vn,'single');
            if Z>max(outp.Zscale)
                disp('Z parameter exceeds data limits');
                Z=max(outp.Zscale);
            elseif Z<min(d(Di).outp.Zscale)
                Z=min(outp.Zscale);
                disp('Z parameter exceeds data limits');
            end
            outp.Zi=find(outp.Zscale>=Z,1,'first');
        tic 
        fprintf('   -converting ACSII: %5.1f%% ',(0)*100);
        for I=1:outp.Sn %reformatting text to numbers
             AA=A((1:outp.Zn)+is1+7+(I-1)*(outp.Zn+7));

             i=round(outp.Zi);
                 %vals(i,I,:)=sscanf(A{i+is1+7+(I-1)*(outp.Zn+7)},'%f');
                 vals(1,I,:)=sscanf(AA{i},'%f');

             fprintf(repmat('\b',1,7));
             fprintf('%5.1f%% ',(I/outp.Sn)*100);         
        end
        fprintf('\r');
        t=toc;
        fprintf('   +done in %.1f seconds \r',t);
        for i=1:size(vars,2);%assigning to variables
        eval(['outp.' vars{i} '.v=vals(:,:,' num2str(i) ')'';']);
%         eval(['outp.' vars{i} '.mean_S=mean(outp.' vars{i} '.v,1);']);
%         eval(['outp.' vars{i} '.mean_Z=mean(outp.' vars{i} '.v,2);']);
%         eval(['outp.' vars{i} '.max_S=max(outp.' vars{i} '.v,[],1);']);
%         eval(['outp.' vars{i} '.max_Z=max(outp.' vars{i} '.v,[],2);']);
        end
        
     else error('numel(Z)>1');
     end
    

    for I=1:outp.Sn
        %outp.current(I,:)=sscanf(A{185+outp.Zn+(I-1)*(outp.Zn+7)},'%f',1);
        outp.current(I,:)=sscanf(A{is1+4+(I-1)*(outp.Zn+7)},'%f',1);
    end
    inp.charge=sum(outp.current*inp.zsep*inp.xlamds/c);
    outp.Sscale=linspace(0,inp.xlamds*inp.zsep*outp.Sn,outp.Sn)';
    
    if numel(Z)==0
        sp1=sum(outp.power.v,1);
        outp.power.peakpos=sum(repmat(outp.Sscale,1,outp.Zn).*outp.power.v,1)./sp1;
        outp.power.std=sqrt(sum(((repmat(outp.Sscale,1,outp.Zn)-repmat(outp.power.peakpos,outp.Sn,1)).^2.*outp.power.v),1)./sp1);
        outp.power.E=outp.power.mean_S*inp.xlamds*inp.zsep*inp.nslice/c;
    else
        sp1=sum(outp.power.v,1);    
        outp.power.peakpos=sum(repmat(outp.Sscale,1,size(outp.power.v,2)).*outp.power.v,1)./sp1;    
        outp.power.std=sqrt(sum(((repmat(outp.Sscale,1,size(outp.power.v,2))-repmat(outp.power.peakpos,outp.Sn,1)).^2.*outp.power.v),1)./sp1);
        outp.power.E=mean(outp.power.v)*inp.xlamds*inp.zsep*inp.nslice/c;
    %outp.power.mean_S_2s=sp1./(max(outp.Sscale)-min(outp.Sscale)).*outp.power.std.*2;
    end
    clear outp.z vals vars i I A AA is1 sp1
    %% Spectrum calculation

        disp(' -calculation of spectrum');

        sc=linspace(-outp.Sn/2,outp.Sn/2,outp.Sn*intN);%original
        %sc=-(outp.Sn-1)/2:1:(outp.Sn-1)/2;

        k0=2*pi/inp.xlamds;
        dk=2*pi/(outp.Sn*inp.xlamds*inp.zsep);
        outp.Kscale=k0+dk*sc';
        outp.Lamdscale=2*pi./outp.Kscale;
        outp.Ln=outp.Sn*intN;
if numel(Z)==0
        outp.spectrum_mid.v=abs(fftshift(fft(sqrt(outp.p_mid.v).*exp(1i*outp.phi_mid.v),intN*outp.Sn,1),1)).^2/sqrt(outp.Sn); %.*exp(-i*outp.phi_mid.v)
        outp.spectrum_mid.mean_S=mean(outp.spectrum_mid.v,1);
        outp.spectrum_mid.mean_Z=mean(outp.spectrum_mid.v,2);
        outp.spectrum_mid.max_S=max(outp.spectrum_mid.v,[],1);
        outp.spectrum_mid.max_Z=max(outp.spectrum_mid.v,[],2);
        outp.spectrum_mid.lamdpos=sum(repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2)).*outp.spectrum_mid.v,1)./sum(outp.spectrum_mid.v,1);
%         outp.spectrum_mid.std_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1))./sum(outp.spectrum_mid.v,1);
        outp.spectrum_mid.std_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1)./sum(outp.spectrum_mid.v,1));
elseif numel(Z)==1
        outp.spectrum_mid.v=abs(fftshift(fft(sqrt(outp.p_mid.v).*exp(1i*outp.phi_mid.v),intN*outp.Sn,1),1)).^2/sqrt(outp.Sn); %.*exp(-i*outp.phi_mid.v)
        outp.spectrum_mid.lamdpos=sum(repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2)).*outp.spectrum_mid.v,1)./sum(outp.spectrum_mid.v,1);
        outp.spectrum_mid.std_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum_mid.v,2))-repmat(outp.spectrum_mid.lamdpos,size(outp.spectrum_mid.v,1),1)).^2.*outp.spectrum_mid.v),1)./sum(outp.spectrum_mid.v,1));
else error('numel(Z)>1');
end
        clear sc k0 dk
        disp('   +done');
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
    
   

%% Data cataloging

d.inp=inp;
d.outp=outp;
d.nm_p=nm_p;
%clear inp outp
fprintf(' +end of import\r');

return
