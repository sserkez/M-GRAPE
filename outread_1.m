% reads .out file
fclose all;
clear all
DiN=1;
nm_p{1}='tdp\1000_u2_tdp_u8_DA.out'; %Data
%nm_p{2}='tdp\1000_u2_tdp_u7_KK.out';

for Di=1:DiN %Data index
    
    fdin=fopen(nm_p{Di},'r');
    A = textscan(fdin,'%s','delimiter','\n');
    A=A{1,1};
    fclose all; clear ans fdin;
    %% Input parameters import
    disp(' -starting import of input parameters');
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
    inp.Zn=round(inp.zstop/inp.xlamd+1);
    inp.Sn=round((numel(A)-(181+inp.Zn+2+7))/(inp.Zn+8));
    disp('   +done');
    clear i is ie k val var
    %% Magnetic parameters import
    disp(' -starting import of magnetic parameters');
    vars={'z','aw','qfld'};
    vals=zeros(inp.Zn,3,'single'); %reserving memory

    for i=1:inp.Zn%reformatting text to numbers
    vals(i,:)=single(sscanf(A{181+i},'%f'));
    end
    for i=1:size(vars,2);%assigning to variables
    eval(['outp.' vars{i} '=vals(:,' num2str(i) ')'';']);
    end
    outp.Zscale=outp.z;
    clear outp.z
    disp('   +done');

    % for i=1:5
    %     if strcmp(sscanf(A{182+inp.Zn+i},'%17c',1),'=================');
    %         break
    %     end
    % end
    % 
    % 181+inp.Zn+2
    % 181+inp.Zn+2+7 початок +(K-1)*(inp.Zn+8)
    % 192+2*inp.Zn+(K-1)*(inp.Zn+8) кінець
    %% Output parameters input
    disp(' -starting import of output parameters');
    text=A{188+inp.Zn};
    k=strfind(text,'-'); text(k)='_';
    k=strfind(text,'<'); text(k)='';
    k=strfind(text,'>'); text(k)='';
    C=textscan(text,'%s');
    vars=C{1}';
    inp.Vn=round(numel(vars));
    clear C k text

    vals=zeros(inp.Zn,inp.Sn,inp.Vn,'single');

    tic 
    for I=1:inp.Sn %reformatting text to numbers
         %for i=181+inp.Zn+2+7+(I-1)*(inp.Zn+8):190+2*inp.Zn+(I-1)*(inp.Zn+8)
         for i=1:188+2*inp.Zn+(I-1)*(inp.Zn+7)-(180+inp.Zn+2+7+(I-1)*(inp.Zn+7))+1    
             %tic
             %val(i,I,:)=single(sscanf(A{i+181+inp.Zn+2+7+(I-1)*(inp.Zn+8)},'%f')); %
             vals(i,I,:)=       sscanf(A{i+180+inp.Zn+2+7+(I-1)*(inp.Zn+7)-1},'%f'); %
             %toc
         %disp(((I-1)/inp.Sn+(i-1)/inp.Zn/inp.Sn)*100);
         end
         fprintf('   -converting ACSII: %.1f%% \r',(I/inp.Sn)*100);
         %fprintf('blablabla %d \r', intvar);
     end
    toc

    for i=1:size(vars,2);%assigning to variables
    eval(['outp.' vars{i} '.v=vals(:,:,' num2str(i) ')'';']);
    eval(['outp.' vars{i} '.mean_S=mean(outp.' vars{i} '.v,1);']);
    eval(['outp.' vars{i} '.mean_Z=mean(outp.' vars{i} '.v,2);']);
    eval(['outp.' vars{i} '.max_S=max(outp.' vars{i} '.v,[],1);']);
    eval(['outp.' vars{i} '.max_Z=max(outp.' vars{i} '.v,[],2);']);
    end

    for I=1:inp.Sn
        outp.current(I,:)=sscanf(A{185+inp.Zn+(I-1)*(inp.Zn+7)},'%f',1);
    end

    outp.Sscale=linspace(0,inp.xlamds*inp.zsep*inp.Sn,inp.Sn)';
    clear outp.z vals vars i I A
    disp('   +done');


    %% Spectrum calculation

    intN=3;
    sc=linspace(-inp.Sn/2,inp.Sn/2,inp.Sn*intN);%original
    %sc=-(inp.Sn-1)/2:1:(inp.Sn-1)/2;

    k0=2*pi/inp.xlamds;
    dk=2*pi/(inp.Sn*inp.xlamds*inp.zsep);
    outp.Kscale=k0+dk*sc';
    outp.Lamdscale=2*pi./outp.Kscale*1e9;

    outp.spectrum.v=abs(fftshift(fft(sqrt(outp.p_mid.v).*exp(1i*outp.phi_mid.v),intN*inp.Sn,1),1)).^2/sqrt(inp.Sn); %.*exp(-i*outp.phi_mid.v)
    outp.spectrum.mean_S=mean(outp.spectrum.v,1);
    outp.spectrum.mean_Z=mean(outp.spectrum.v,2);
    outp.spectrum.max_S=max(outp.spectrum.v,[],1);
    outp.spectrum.max_Z=max(outp.spectrum.v,[],2);
    outp.spectrum.lamdpos=sum(repmat(outp.Lamdscale,1,size(outp.spectrum.v,2)).*outp.spectrum.v,1)./sum(outp.spectrum.v,1);
    outp.spectrum.rms_lamd=sqrt(sum(((repmat(outp.Lamdscale,1,size(outp.spectrum.v,2))-repmat(outp.spectrum.lamdpos,size(outp.spectrum.v,1),1)).*outp.spectrum.v).^2,1))./sum(outp.spectrum.v,1);
    clear sc k0 dk 
    %%
    
    d(Di).inp=inp;
    d(Di).outp=outp;
    d(Di).nm_p=nm_p(Di);
    %clear inp outp
end

%% plots
Di=1;

fig1=figure(1);
set(fig1,'name',['Electrons ',nm_p{Di}],'numbertitle','off');
H(Di).h1=subplot(4,1,1);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.aw,d(Di).outp.Zscale,d(Di).outp.qfld);
set(haxes,{'ycolor'},{'b';'r'})
ylabel(haxes(1),'Undulator\newlinea_W');
ylabel(haxes(2),'Quadrupole\newlinedB/dx [T/m]');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2,'color','r');
for i=1:2
    ax(i) = get(hline(i),'Parent');
    set(ax(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
clear ax
%---------
H(Di).h2=subplot(4,1,2);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.xrms.mean_S,d(Di).outp.Zscale,d(Di).outp.x.mean_S);
hline(3)=line(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S, 'Parent', haxes(1), 'linestyle','--','color','b','linewidth',2.5);
hline(4)=line(d(Di).outp.Zscale,d(Di).outp.y.mean_S, 'Parent', haxes(2), 'linestyle','--','color',[0 .5 0],'linewidth',1.5);
ylabel(haxes(1),'Beam size\newline\sigma [m]');
ylabel(haxes(2),'Beam centroid\newline-x - -y [m]');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',1.5,'color','b','linewidth',2.5);
set(hline(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',1.5);
for i=1:2
    ax(i) = get(hline(i),'Parent');
    set(ax(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
clear ax
%---------
H(Di).h3=subplot(4,1,3);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.energy.mean_S,d(Di).outp.Zscale,d(Di).outp.e_spread.mean_S);
ylabel(haxes(1),'energy (\Delta\gamma)');
ylabel(haxes(2),'spread (\sigma_\gamma)');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2);
for i=1:2
    ax(i) = get(hline(i),'Parent');
    set(ax(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
clear ax
%---------
H(Di).h3=subplot(4,1,4);
haxes = plot(d(Di).outp.Zscale,d(Di).outp.bunching.mean_S,'LineWidth',2,'color','r');
ylabel('bunching\newline|<exp(i \theta)>|');
xlabel('z [m]');
xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
    
fig2=figure(2);
set(fig2,'name',['Photons ',nm_p{Di}],'numbertitle','off');

h1=subplot(3,1,1);
semilogy(d(Di).outp.Zscale,d(Di).outp.power.max_S,'LineWidth',1,'linestyle','--','color','k');
hold all
semilogy(d(Di).outp.Zscale,d(Di).outp.power.mean_S,'LineWidth',2,'color','r');
hold off
ylim([min(d(Di).outp.power.mean_S) max(d(Di).outp.power.max_S)]);
set(gca,'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
ylabel('P [W]');
%xlabel('z [m]');
xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
legend('max','mean','location','NorthWest');

h2=subplot(3,1,2);
%semilogy(d(Di).outp.Zscale,d(Di).outp.spectrum.max_S,'LineWidth',2);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,outp.spectrum.max_S,d(Di).outp.Zscale,d(Di).outp.spectrum.rms_lamd,'semilogy','plot');

%plot(d(Di).outp.Zscale,d(Di).outp.spectrum.rms_lamd,'LineWidth',2);
%hold all
ylabel(haxes(1),'P_{max}(\lambda) [a.u.]') % label left y-axis
ylabel(haxes(2),'\sigma_{P(\lambda)} [nm]') % label right y-axis
set(hline(1),'LineWidth',1.5,'color','b','linewidth',2);
set(hline(2),'LineWidth',1.5,'color',[0 0.5 0],'linewidth',2);
%xlabel(haxes(2),'z [m]') % label x-axis
axis tight
set(haxes(1),'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
        %set(haxes(i),'axis','tight');
        axis tight
end


h3=subplot(3,1,3);
plot(d(Di).outp.Zscale,d(Di).outp.r_size.mean_S,'LineWidth',2);
ylabel('\sigma_{rad} [m]');
xlabel('z [m]');
xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
axis tight

%% ################################

Z=0;

for Di=1:DiN
    Zi(Di)=find(d(Di).outp.Zscale>=Z,1,'first');
end
    %Z=0;


%     figure(54);
%     %plot(linspace(max(d(Di).outp.Lamdscale),min(d(Di).outp.Lamdscale),intN*inp.Sn),d(Di).outp.spectrum.max_Z);
%     plot(d(Di).outp.Lamdscale,d(Di).outp.spectrum.max_Z);
%     xlim([1.238 1.24]);
% 
%     figure(55);
%     plot(d(Di).outp.Zscale,d(Di).outp.spectrum.max_S);
    
    if DiN==2
        
    figure(56);  
    plot(d(1).outp.Lamdscale,d(1).outp.spectrum.v(:,Zi(1)),'linewidth',2);%/max(d(1).outp.spectrum.v(:,Zi(1)))
    hold all
    plot(d(2).outp.Lamdscale,d(2).outp.spectrum.v(:,Zi(2)),'linewidth',2,'color','r','linestyle','--');%/max(d(2).outp.spectrum.v(:,Zi(2)))
    hold off
    legend('direct\newlineapproach','phenomenological\newlineapproach','Location','northwest')
    xlim([1.2385 1.2395]);
    %xlim([1.764 1.766]);
    ylabel('P(\lambda) [a.u.]');
    xlabel('\lambda [nm]');
%     figure(57);
%     plot(d(Di).outp.Sscale,d(Di).outp.power.v(:,Zi(Di)));
%     figure(58);
%     plot(d(Di).outp.Sscale,d(Di).outp.phi_mid.v(:,Zi(Di)));

    figure(59);
    semilogy(d(1).outp.Zscale,d(1).outp.power.mean_S,'LineWidth',2);
    hold all
    semilogy(d(2).outp.Zscale,d(2).outp.power.mean_S,'LineWidth',2,'color','r','linestyle','--');
    hold off
%     legend('direct approach','phenomenological approach','Location','northwest')
    legend(nm_p{1},nm_p{2},'Location','northwest')
    ylabel('P [W]');
    xlabel('z [m]');
    xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
    
    end
    
    %plot pulse evolution
    Zgrid_p=repmat(d(Di).outp.Zscale,d(Di).inp.Sn,1);
    Sgrid=repmat(d(Di).outp.Sscale,1,d(Di).inp.Zn);
    Lamdgrid=repmat(d(Di).outp.Lamdscale,1,d(Di).inp.Zn);
    P=double(d(Di).outp.power.v);
    S=double(d(Di).outp.spectrum.v);
    
%     figure(611);
%     h=surf(Zgrid,Sgrid,P,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
    Zgrid_s=repmat(d(Di).outp.Zscale,numel(d(Di).outp.Lamdscale),1);
%     figure(612);
%     h=surf(Zgrid,Lamdgrid,S,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
     Zm=20; Sm=5;
     [Zgrid_p_n, Sgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).inp.Zn/Zm), linspace(min(d(Di).outp.Sscale),max(d(Di).outp.Sscale),d(Di).inp.Sn/Sm));
     P_n=interp2(Zgrid_p,Sgrid,P,Zgrid_p_n,Sgrid_n);
     figure(621);
     h=surf(Zgrid_p_n,Sgrid_n,P_n,log(P_n),'linestyle','none');
     set(get(h,'Parent'),'ZScale','log');
     title('power');
     axis tight
    
     
     Zm=20; Lm=2;
    [Zgrid_s_n, Lamdgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).inp.Zn/Zm), linspace(min(d(Di).outp.Lamdscale),max(d(Di).outp.Lamdscale),numel(d(Di).outp.Lamdscale)/Lm));
    S_n=interp2(Zgrid_s,Lamdgrid,S,Zgrid_s_n,Lamdgrid_n);
    figure(622);
    h=surf(Zgrid_n,Lamdgrid_n,S_n,log(S_n),'linestyle','none');
    set(get(h,'Parent'),'ZScale','log');
    title('spectrum');
    axis tight
    %xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);

% [haxes,hline(1),hline(2)] = plotyy(outp.Zscale,outp.power.mean_S,outp.Zscale,outp.power.mean_S,'semilogy','plot');
% %set(haxes,{'ycolor'},{'b';'r'})
% ylabel(haxes(1),'P [W] (log)');
% ylabel(haxes(2),'P [W] (lin)');
% xlabel(haxes(2),'z [m]');
% set(hline(1),'LineWidth',2);
% set(hline(2),'LineWidth',1,'color',[0 0.5 0]);
% for i=1:2
%     ax(i) = get(hline(i),'Parent');
%     set(ax(i),'XLim',[outp.Zscale(1) outp.Zscale(end)]);
% end
% clear ax


break
%%
h2=subplot(3,1,2);
[haxes,hline(1),hline(2)] = plotyy(outp.Zscale,outp.xrms.mean_S,outp.Zscale,outp.x.mean_S);
hline(3)=line(outp.Zscale,outp.yrms.mean_S, 'Parent', haxes(1), 'linestyle','--','color','b','linewidth',2);
hline(4)=line(outp.Zscale,outp.y.mean_S, 'Parent', haxes(2), 'linestyle','--','color',[0 .5 0],'linewidth',2);
ylabel(haxes(1),'\sigma [m]');
ylabel(haxes(2),' -x --y [m]');
xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',1.5,'color','b','linewidth',2);
set(hline(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',2);
for i=1:2
    ax(i) = get(hline(i),'Parent');
    set(ax(i),'XLim',[outp.Zscale(1) outp.Zscale(end)]);
end
clear ax

h3=subplot(3,1,3);
[haxes,hline(1),hline(2)] = plotyy(outp.Zscale,outp.energy.mean_S,outp.Zscale,outp.e_spread.mean_S);
ylabel(haxes(1),'energy (\Delta\gamma)');
ylabel(haxes(2),'\sigma_\gamma');
xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2);
for i=1:2
    ax(i) = get(hline(i),'Parent');
    set(ax(i),'XLim',[outp.Zscale(1) outp.Zscale(end)]);
end
clear ax

break

subfigure(1,3,1);
hold all
plot(outp.Zscale,outp.aw,linewidth,'2');
plot(outp.Zscale,outp.qfld,linewidth,'2');
hold off