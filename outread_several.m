% reads .out file
fclose all;
clear all
%N=10;'
n1=1;
n2=100;
%n2=5;
n=n1:n2;
%     n=[1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
%  n=[6 7 8 9];
s1=1.8e-6;
s2=2.3e-6; % XFEL

% s1=10e-6;
% s2=15e-6;% LCLS

dflread=0;



%for i=n1:n2
for i=n
      %nm_p{i}=['C:\-D-\Work\SASE3_chicane\run6\U1.',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\SASE3_chicane\run6\U2.',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\530_tdp\530_u1_tdp_',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u2_28_',num2str(i),'.out'];
%       nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_fc_nn_ps_',num2str(i),'_sc10.out'];
%       nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_old\930_sl_u2_',num2str(i),'_sc10.out'];
    
    nm_p{i}=['C:\-D-\Work\SASE3_chicane\630_250\statrun1\U1.',num2str(i),'.out'];
    nm_p{i}=['C:\-D-\Work\SASE3_chicane\630_250\statrun1\U2.',num2str(i),'.out'];


    %nm_p{1}='C:\-D-\Work\SASE3_chicane\630_250\statrun\U2.1.out';
%      nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_',num2str(i),'.out'];
%        nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u2_',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u1_',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated\930_u1_',num2str(i),'.out'];
%         nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_mm\930_sl_u2_',num2str(i),'_sc.out'];
%         nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_old\930_sl_u2_',num2str(i),'_sc10.out'];
%       nm_p{i}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\700_u2_',num2str(i),'_n.out'];
%      %  nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_4\U2_t.',num2str(i),'.out'];
%         nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_3\U1.',num2str(i),'.out'];
%         nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_4\U2_short.',num2str(i),'.out'];
% %       nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_3\U2_t_noise.',num2str(i),'.out'];
% %            nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_4\U2_t.',num2str(i),'.out'];
% %         nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_3\U1_sase.',num2str(i),'.out'];
%        % nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_3_tap\U1.',num2str(i),'.out'];
%  nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_3\U1_sase.',num2str(i),'.out'];
% % nm_p{i}=['C:\-D-\Work\SASE3_SXRSS\tdp_5\U1.',num2str(i),'.out'];
end

%nm_p{1}=['C:\-D-\Work\SASE3_chicane\run5\U2.1.out'];
%nm_p{1}=['C:\-D-\Work\LCLS\tmp\930_u2_tdp_40.out'];

   %break    
%  DiN=size(nm_p,2);
% DiN=1;    

% for Di=n1:n2 %Data index
index=0;
for Di=n
    index=index+1;
   %try load([nm_p{Di},'.mat']); 
   %catch error
%     try 
disp([' ']);
disp(['-file ',num2str(index),' of ',num2str(numel(n))]);
d(Di)=outread(nm_p{Di},1,0,3); 
%     %    d(Di).outp.increment.v_last=d(Di).outp.increment.v(:,end); d(Di).outp.increment.v(:,1:end-1)=[];
%     %    d(Di).outp.p_mid.v_last=d(Di).outp.p_mid.v(:,end); d(Di).outp.p_mid.v(:,1:end-1)=[];
%     %    d(Di).outp.phi_mid.v_last=d(Di).outp.phi_mid.v(:,end); d(Di).outp.phi_mid.v(:,1:end-1)=[];
%     %    d(Di).outp.r_size.v_last=d(Di).outp.r_size.v(:,end); d(Di).outp.r_size.v(:,1:end-1)=[];
%        catch  error
%        end
%        save([nm_p{Di},'.mat'],'d');

end
%%

color1_shot=[0.8 0.8 0.8];
color1_select=[0.5 0.5 0.5]; %grey
color1_select=color1_shot;
color1_mean=[0.2 0.2 0.2];

% color1_shot=[0.6 0.8 0.6];
% color1_select=[0.6 0.8 0.6]; %green
% color1_mean=[0. 0.3 0.];
%
% color1_shot=[0.8 0.6 0.6];
% color1_select=[0.8 0.6 0.6]; %red
% color1_mean=[0.3 0. 0.];

% color1_select=[1 1 1];
% color1_shot=[1 1 1];

%  color1_mean=[0.1 0.7 0.1]; %green
%  color1_mean=[0.7 0.1 0.1]; %red

dohold=0;
Z=50;
for i=[880 887 886 888 889 8900 890 891 892 8875 887 8890 8910 8861 8881 9001 9002 9003]
    if dohold
        figure(i);
        hold on
    else
        figure(i);
        hold off
    end
end

for Di=n
    
    Zi(Di)=find(d(Di).outp.Zscale<=Z,1,'last');
%     if isempty(Zi(Di))
%         Zi(Di)=d(Di).outp.Zn;
%     end

Ex(Di,:)=sqrt(d(Di).outp.p_mid.v(:,Zi(Di))).*exp(1i*d(Di).outp.phi_mid.v(:,Zi(Di)));
% addpath('c:\-D-\Work\LCLS\New_matlab\tftb-0.2\mfiles');

% [W,~,Freqscale] = tfrwv(Ex,d(Di).outp.Zscale*1e-6/3e8);
% break
% figure(6996)
% imagesc(d(Di).outp.Zscale,Freqscale,W);

% W(Di,:,:)=mywigner(Ex(Di,:)');
% disp(Di)


    Pse(Di-n1+1,:)=d(Di).outp.power.v(:,Zi(Di)); %Power over S at Z    
    Ez(Di-n1+1,:)=d(Di).outp.power.E;
    Sscale(Di-n1+1,:)=d(Di).outp.Sscale;
    Zscale(Di-n1+1,:)=d(Di).outp.Zscale;
    Se_all(Di-n1+1,:,:)=d(Di).outp.spectrum_mid.v; %%%%%% remove
    Se(Di-n1+1,:)=Se_all(Di-n1+1,:,Zi(Di));
    %Se(Di-n1+1,:)=d(Di).outp.spectrum_mid.v(:,Zi(Di));
    Lscale(Di-n1+1,:)=d(Di).outp.Lamdscale;
    
    
    try
    if d(Di).parm.dflpresent
    Se_dfl(Di-n1+1,:)=d(Di).dfl.spectrum; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% remove scaling!!!!!
    Pse_dfl(Di-n1+1,:)=d(Di).dfl.power;
    %Se(Di-n1+1,:)=d(Di).outp.spectrum_mid.v(:,Zi(Di));
    Lscale_dfl(Di-n1+1,:)=d(Di).dfl.Lamdscale;
    Sscale_dfl=d(Di).dfl.Sscale;
    Xt(:,:,Di-n1+1)=d(Di).dfl.X;
    Xfz(:,:,Di-n1+1)=d(Di).dfl.Xfz;
    end
    catch error
    end;

s2i=find(d(Di).outp.Sscale<=s2,1,'last');
s1i=find(d(Di).outp.Sscale>=s1,1,'first');

E_beam_energy_mean(Di-n1+1,:)=mean(d(Di).outp.energy.v(s1i:s2i,:)+d(Di).inp.gamma0,1);
E_beam_energy_spread(Di-n1+1,:)=mean(d(Di).outp.e_spread.v(s1i:s2i,:),1);

Rad_size_transv(Di-n1+1,:)=d(Di).outp.r_size.mean_S_norm.*2*1e6;
Rad_size_long(Di-n1+1,:)=d(Di).outp.power.std.*2*1e6;

% %Pz_s(Di-n1+1,:)=mean(d(Di).outp.power.v(s1i:s2i,:),1);
%Pz_s(Di,:)=mean(d(Di).outp.power.v(s1i:s2i,:),1); %!!!!!!!!!! calculates within a window!!!!!!!
 Pz_s(Di,:)=max(d(Di).outp.power.v,[],1); %check if is relevant
Pl(Di,:)=d(Di).outp.spectrum_mid.max_S;

    if exist([nm_p{Di},'.dfl'],'file') && dflread

        disp('-importing *.dfl files');
        disp(['-file ',num2str(Di),' of ',num2str(numel(n))]);
        disp(' ');

        Nn=floor(d(Di).inp.nslice+d(Di).inp.zstop/d(Di).inp.xlamd/d(Di).inp.zsep);
        [XX,N]=fieldimport_all([nm_p{Di},'.dfl'],d(Di).inp.ncar,1);
        %XX=prop_TF(XX,d(Di).inp.dgrid*2,d(Di).inp.xlamds,0);
        dflSscale=linspace(0,d(Di).inp.xlamds*d(Di).inp.zsep*N,N);

        P0=sum(sum(sum(abs(XX).^2)))/N*abs(dflSscale(end)-dflSscale(1))/3e8;

    %    fieldplot3d(110,XX,d(Di).inp.leng,1,dflSscale,[nm_p{Di},'.dfl'],1);
        [X2,dflLamdscale]=dfl_time2freq(XX,dflSscale,d(Di).inp.xlamds);
        %Se_dfl(Di-n1+1,:)=reshape(mean(mean(abs(X2).^2,1),2),1,[]);
        Se_dfl(Di-n1+1,:)=reshape(sum(sum(abs(X2).^2,1),2),1,[]);
        Pse_dfl(Di-n1+1,:)=reshape(sum(sum(abs(XX).^2,1),2),1,[]);
        Lscale_dfl=dflLamdscale;
        Sscale_dfl=dflSscale;
        clear X2 XX
    end
    
dfl_leng(Di-n1+1)=d(Di).inp.leng;
xlamds(Di-n1+1)=d(Di).inp.xlamds;
end

    Zi(Zi==0)=[];
    Pse(sum(Pse,2)==0,:)=[];
    Pz_s(sum(Pz_s,2)==0,:)=[];
    Pl(sum(Pl,2)==0,:)=[];
    Ez(sum(Ez,2)==0,:)=[];
    Se(sum(Se,2)==0,:)=[];
    Lscale(sum(Lscale,2)==0,:)=[];
    Sscale(sum(Sscale,2)==0,:)=[];
    Zscale(sum(Zscale,2)==0,:)=[];
    Se_all(sum(sum(Se_all,3),2)==0,:,:)=[];

    Se_sum_all=sum(permute(Se_all,[3 2 1]),3);

    Lscale=Lscale(1,:);
    Sscale=Sscale(1,:);
    Zscale=Zscale(1,:);
    
%     if mean(dfl_leng(n(1):n(end)))-dfl_leng(n(1))>2e-5
%         error('transverse mesh sie is dirrefent');
%     else
%         dfl_leng=(dfl_leng(n(1)));
%     end
    
    if mean(xlamds)-xlamds(1)>1e-5
        error('transverse mesh sie is dirrefent');
    else
        xlamds=xlamds(1);
    end
    
    if exist([nm_p{Di},'.dfl'],'file') && dflread
        Xf=sqrt(sum(abs(Xfz).^2,3)).*exp(1i.*angle(mean(Xfz,3)));
        X=sqrt(sum(abs(Xt).^2,3)).*exp(1i.*angle(mean(Xt,3)));
    end
%  %%
if numel(Se)~=0 && 1
    figure(887)
    plot(Lscale.*1e9,Se,'linewidth',1,'color',color1_shot);
    hold on
    plot(Lscale.*1e9,Se(1,:),'linewidth',0.5,'color',color1_select);
    plot(Lscale.*1e9,mean(Se,1),'linewidth',2,'color',color1_mean);
    hold off
    axis tight
    ylabel('P(\lambda) [arb.units]');
    xlabel('\lambda [nm]');

    xlamds./findFWHM(Lscale,mean(Se,1));
    
    L0=1239.8/Lscale(mean(Se,1)==max(mean(Se,1)))*1e-9;
    %L0=699;
    %L0=500;
    figure(886)
    %plot(1239.8./(Lscale.*1e9)-L0,Se,'linewidth',1,'color',color1_shot); %uncomment!!
    hold on
    %plot(1239.8./(Lscale.*1e9)-L0,Se(1,:),'linewidth',0.5,'color',color1_select);
    plot(1239.8./(Lscale.*1e9)-L0,mean(Se,1),'linewidth',2,'color',color1_mean);
    hold off
    axis tight

    try
        spec_resol=L0/findFWHM(1239.8./(Lscale.*1e9)-L0,mean(Se,1));
            findFWHM(1239.8./Lscale_dfl(1,:)*1e-9,mean(Se_dfl,1));
    catch
        spec_resol=0;
    end
%findFWHM(1239.8./Lscale_dfl(1,:),mean(Se_dfl,1));

    
else
    try
        close 886 
    end
    try
        close 887
    end
end



ylabel('P(\lambda) [arb. units]');
xlabel('Photon energy [eV]');
xlim([-2 2]);




if exist('Lscale_dfl','var')
    fig=figure(8861);
    set(fig,'name','Lscale_dfl','numbertitle','on');
    L0_dfl=1239.8/Lscale_dfl(1,mean(Se_dfl,1)==max(mean(Se_dfl,1)))*1e-9;
    plot(1239.8./(Lscale_dfl(1,:).*1e9)-L0_dfl,Se_dfl,'linewidth',1,'color',color1_shot);
    hold on
    plot(1239.8./(Lscale_dfl(1,:).*1e9)-L0_dfl,Se_dfl(1,:),'linewidth',0.5,'color',color1_select);
    plot(1239.8./(Lscale_dfl(1,:).*1e9)-L0_dfl,mean(Se_dfl,1),'linewidth',2,'color',color1_mean);
%     plot(1239.8./(Lscale_dfl(1,:).*1e9)-L0_dfl,mean(Se_dfl,1),'linewidth',2,'color','b','linestyle','--');
    hold off
    axis tight
    ylabel('P(\lambda) [arb. units]');
    xlabel('Photon energy [eV]');
    %xlim([2.3 2.45]); %500eV
    %xlim([0.81 0.84]); %1500eV
    %xlim([0.76 0.79]); %2000eV
    %xlim([1.85 2.05]); %500eV
    %xlim([-1 1]);
    try
        spec_resol_dfl=L0_dfl/findFWHM(1239.8./(Lscale_dfl.*1e9)-L0_dfl,mean(Se_dfl,1));
    catch
        spec_resol_dfl=0;
    end
    
    fig=figure(8881);
    set(fig,'name','Sscale_dfl','numbertitle','on');
    plot(Sscale_dfl.*1e6,Pse_dfl,'linewidth',1,'color',color1_shot);
    hold on
    plot(Sscale_dfl.*1e6,Pse_dfl(1,:),'linewidth',0.5,'color',color1_select);
    plot(Sscale_dfl.*1e6,mean(Pse_dfl,1),'linewidth',2,'color',color1_mean);
 %   plot(Sscale_dfl.*1e6,mean(Pse_dfl,1),'linewidth',2,'color','b','linestyle','--');
    hold off
    axis tight
    ylabel('P [W]');
    xlabel('s [\mum]');
    %xlim([2.3 2.45]); %500eV
    %xlim([0.81 0.84]); %1500eV
    %xlim([0.76 0.79]); %2000eV
    %xlim([1.85 2.05]); %500eV

    
%     fieldplot(66,X,dfl_leng,'x_dist',1);
%     fieldplot_ff(67,Xf,dfl_leng,xlamds,'phi_dist',1);
    
    if exist([nm_p{Di},'.dfl'],'file') && dflread
        Ix=sum(sum(abs(X).^2,3),1)'./sum(sum(sum(abs(X).^2)));
        Iy=sum(sum(abs(X).^2,3),2)./sum(sum(sum(abs(X).^2)));
        Ix=Ix./max(Ix);
        Iy=Iy./max(Iy);
    
        Ifx=sum(sum(abs(Xf).^2,3),1)'./sum(sum(sum(abs(Xf).^2)));
        Ify=sum(sum(abs(Xf).^2,3),2)./sum(sum(sum(abs(Xf).^2)));
        Ifx=Ifx./max(Ifx);
        Ify=Ify./max(Ify);

    
    [Mx,~,~]=size(X);
    dx=dfl_leng/Mx;
    xscale=((Mx-1)/2+1-(1:Mx))*dx*1000;
    xfscale=linspace(-1/(2*dx)+1/dfl_leng/2,1/(2*dx)-1/dfl_leng/2,Mx).*xlamds*1000;
    
    fig=figure(1111);
    set(fig,'name','trans_dist_xy','numbertitle','on');
    plot(xscale,Ix,'linewidth',2,'color','b');
    hold on
    plot(xscale,Iy,'linewidth',2,'color','r','linestyle','--');
    hold off
    axis tight
    ylabel('P [arb.units]');
    xlabel('x,y [mm]');
    h_leg=legend('X profile','Y profile');
    set(h_leg,'FontSize',8);
    xlim([-0.6 0.6]);
    
    fig=figure(1112);
    set(fig,'name','trans_dist_fxfy','numbertitle','on');
    plot(xfscale,Ifx,'linewidth',2,'color','b');
    hold on
    plot(xfscale,Ify,'linewidth',2,'color','r','linestyle','--');
    hold off
    axis tight
    ylabel('P [arb.units]');
    xlabel('\theta [mrad]');
    h_leg=legend('X divergence','Y divergence');
    set(h_leg,'FontSize',8);
    xlim([-0.03 0.03]);
    
    end
    
    %Xf=fftshift(fft2(ifftshift(X)));
    %fieldplot_ff(66,Xf,dfl_leng,xlamds,'far field',1)
    %fieldplot_ff(67,Xfz,dfl_leng,xlamds,'true far field',1)
end

fig=figure(888);
set(fig,'name','axis Power along beam (Pse)','numbertitle','on');
plot(Sscale(1,:).*1e6,Pse,'linewidth',1,'color',color1_shot);   %uncomment!!
hold on
plot(Sscale(1,:).*1e6,Pse(1,:),'linewidth',0.5,'color',color1_select);
plot(Sscale(1,:).*1e6,mean(Pse,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('P [W]');
xlabel('s [\mum]');
%xlim([1 5]);


fig=figure(889);
set(fig,'name','Ez_log','numbertitle','on');
semilogy(Zscale(1,:),Ez,'linewidth',1,'color',color1_shot);
hold on
semilogy(Zscale(1,:),Ez(1,:),'linewidth',2,'color',color1_select);
semilogy(Zscale(1,:),mean(Ez,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('E [J]');
xlabel('z [m]');

fig=figure(880);
set(fig,'name','Ez_var','numbertitle','on');
Ez1=Ez./(ones(size(Ez))*diag(mean(Ez,1)));
Ez_var=mean(Ez1-(ones(size(Ez1))*diag(mean(Ez1,1))).^2,1);
Ez_var=var(Ez1,1);
plot(Zscale(1,:),Ez_var,'linewidth',2,'color',color1_mean);
axis tight
ylabel('Variance of energy fluctuations');
xlabel('z [m]');

fig=figure(8890);
set(fig,'name','Ez','numbertitle','on');
plot(Zscale(1,:),Ez,'linewidth',1,'color',color1_shot);
 hold on
% plot(SScale(1,:),Pse(1,:),'linewidth',0.5,'color',color1_select);
plot(Zscale(1,:),mean(Ez,1),'linewidth',2,'color',color1_mean);
 hold off
 axis tight
 ylabel('E [J]');
 xlabel('z [m]');

 if max(Ez(:,1))~=0
     fig=figure(8900);
     set(fig,'name','Ez_0_run','numbertitle','on');
     plot(Ez(:,1));
     ylabel('E [J]');
     xlabel('Run #');
     ylim([0 max(Ez(:,1))]);
 end
 
 fig=figure(890);
 set(fig,'name','Ez_Z_run','numbertitle','on');
 plot(Ez(:,Zi(1)));
 ylabel('E [J]');
 xlabel('Run #');
 ylim([0 max(Ez(:,Zi(1)))]);
 
 fig=figure(8901);
 set(fig,'name','Pz_s_Z_run','numbertitle','on');
 plot(Pz_s(:,Zi(1)));
 ylabel('P [W]');
 xlabel('Run #');
 ylim([0 max(Pz_s(:,Zi(1)))]);
 
 fig=figure(891);
 set(fig,'name','Pz_s_log','numbertitle','on');
 semilogy(Zscale(1,:),Pz_s,'linewidth',1,'color',color1_shot); %uncomment!!
 hold on
 semilogy(Zscale(1,:),Pz_s(1,:),'linewidth',1,'color',color1_select);
 semilogy(Zscale(1,:),mean(Pz_s,1),'linewidth',2,'color',color1_mean);
 hold off
 axis tight
 ylabel('P [W]');
 xlabel('z [m]');
 
  fig=figure(8910);
  set(fig,'name','Pz_s','numbertitle','on');
 plot(Zscale(1,:),Pz_s,'linewidth',1,'color',color1_shot); %uncomment!!
 hold on
 plot(Zscale(1,:),Pz_s(1,:),'linewidth',1,'color',color1_select);
 plot(Zscale(1,:),mean(Pz_s,1),'linewidth',2,'color',color1_mean);
 hold off
 axis tight
 ylabel('P [W]');
 xlabel('z [m]');
 
  fig=figure(892);
  set(fig,'name','Pl','numbertitle','on');
semilogy(Zscale(1,:),Pl,'linewidth',1,'color',color1_shot);
 hold on
 plot(Zscale(1,:),Pl(1,:),'linewidth',1,'color',color1_select);
 semilogy(Zscale(1,:),mean(Pl,1),'linewidth',2,'color',color1_mean);
 hold off
 axis tight
 ylabel('P_{max}(\lambda) [arb. units]');
 xlabel('z [m]');

%  Li=find(d(Di).outp.Lamdscale<=1.3341e-9,1,'first');
%  Se_variance=reshape(std(Se_all(:,Li,:),1),1,[])./reshape(mean(Se_all(:,Li,:),1),1,[]);
%  E_variance=std(Pz_s,1)./mean(Pz_s,1);
%  figure(8875)
%  plot(Zscale,Se_variance);
%  ylim([0,1]);
 figure(4521)
 hist(Ez(:,Zi(1)),10);
 xlabel('shot energy [J]');
 ylabel('number of shots');
 %
 %figure(886)

 fig=figure(9001);
set(fig,'name','average electron beam energy within a window','numbertitle','on');
plot(Zscale(1,:).*1e6,E_beam_energy_mean,'linewidth',1,'color',color1_shot);   %uncomment!!
hold on
plot(Zscale(1,:).*1e6,E_beam_energy_mean(1,:),'linewidth',0.5,'color',color1_select);
plot(Zscale(1,:).*1e6,mean(E_beam_energy_mean,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('<\gamma>');
xlabel('z [m]');
%xlim([1 5]);

 fig=figure(9002);
set(fig,'name','average electron beam energy spread within a window','numbertitle','on');
plot(Zscale(1,:).*1e6,E_beam_energy_spread,'linewidth',1,'color',color1_shot);   %uncomment!!
hold on
plot(Zscale(1,:).*1e6,E_beam_energy_spread(1,:),'linewidth',0.5,'color',color1_select);
plot(Zscale(1,:).*1e6,mean(E_beam_energy_spread,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('<\delta\gamma>');
xlabel('z [m]');
%xlim([1 5]);
 
fig=figure(9003);
set(fig,'name','correlation','numbertitle','on');
scatter(E_beam_energy_spread(:,1),Pz_s(:,end));   %uncomment!!
axis tight
xlabel('<\delta\gamma_{init}>');
ylabel('<Power_{final}>');
%xlim([1 5]);

fig=figure(9101);
set(fig,'name','Radiaton_size_long','numbertitle','on');
plot(Zscale(1,:).*1e6,Rad_size_long,'linewidth',1,'color',color1_shot);   %uncomment!!
hold on
plot(Zscale(1,:).*1e6,Rad_size_long(1,:),'linewidth',0.5,'color',color1_select);
plot(Zscale(1,:).*1e6,mean(Rad_size_long,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('<\delta s_l [\mum]>');
xlabel('z [m]');

fig=figure(9102);
set(fig,'name','Radiaton_size_transv','numbertitle','on');
plot(Zscale(1,:).*1e6,Rad_size_transv,'linewidth',1,'color',color1_shot);   %uncomment!!
hold on
plot(Zscale(1,:).*1e6,Rad_size_transv(1,:),'linewidth',0.5,'color',color1_select);
plot(Zscale(1,:).*1e6,mean(Rad_size_transv,1),'linewidth',2,'color',color1_mean);
hold off
axis tight
ylabel('<\delta s_t [\mum]>');
xlabel('z [m]');


findFWHM(3e8./Lscale,Se_sum_all(end,:))*findFWHM(Sscale(1,:).*3.*1e-15*1e6,mean(Pse,1)); %coefficient of transform-limitation ?

set(figure(891), 'Position', [300, 300, 450, 350]);  
set(figure(889), 'Position', [300, 300, 450, 350]);
set(figure(892), 'Position', [300, 300, 450, 350]);  

set(figure(880), 'Position', [300, 200, 450, 350]);  
set(figure(887), 'Position', [300, 200, 450, 350]);  

set(figure(886), 'Position', [100, 100, 450, 350]);
set(figure(888), 'Position', [100, 100, 450, 350]);

set(figure(8890), 'Position', [400, 100, 450, 350]);
set(figure(8910), 'Position', [400, 100, 450, 350]);

set(figure(4521), 'Position', [450, 100, 500, 350]);

if exist('Lscale_dfl','var')
    set(figure(8861), 'Position', [200, 200, 450, 350]);
    set(figure(8881), 'Position', [200, 200, 450, 350]);
    
    set(figure(1111), 'Position', [200, 300, 450, 350]);
    set(figure(1112), 'Position', [200, 300, 450, 350]);
end
%%% Wigner

% Zi=find(d(Di).outp.Zscale<=Z,1,'last');
% Ex=sqrt(d(Di).outp.p_mid.v(:,Zi)).*exp(1i*d(Di).outp.phi_mid.v(:,Zi));
% figure(6995)
% plot(d(Di).outp.Sscale,abs(Ex).^2);
% % addpath('c:\-D-\Work\LCLS\New_matlab\tftb-0.2\mfiles');
% 
% % [W,~,Freqscale] = tfrwv(Ex,d(Di).outp.Zscale*1e-6/3e8);
% % break
% % figure(6996)
% % imagesc(d(Di).outp.Zscale,Freqscale,W);
% 
% W=mywigner(Ex);

% % % % % % % Sn=d(Di).outp.Sn;
% % % % % % % sc=linspace(-Sn/2,Sn/2,Sn);%original
% % % % % % % %sc=-(outp.Sn-1)/2:1:(outp.Sn-1)/2;
% % % % % % % 
% % % % % % % k0=2*pi/d(Di).inp.xlamds;
% % % % % % % dk=2*pi/(d(Di).outp.Sn*d(Di).inp.xlamds*d(Di).inp.zsep);
% % % % % % % Kscale=k0+dk*sc';
% % % % % % % Lamdscale=flipud(2*pi./Kscale);
% % % % % % % clear Kscale
% % % % % % % mW=mean(W,1);
% % % % % % % mW=permute(mW,[3,2,1]);
% % % % % % % figure(6996)
% % % % % % % imagesc(d(Di).outp.Sscale,Lamdscale,mW);
% % % % % % % ylim([1.775e-9 1.782e-9]);
% % % % % % % 
% % % % % % % figure(6997)
% % % % % % % plot(Lamdscale,mean(mW,2)/max(sum(mW,2)))
% % % % % % % axis tight
% % % % % % % 
% % % % % % % figure(6998)
% % % % % % % plot(Lscale.*1e9,mean(Se,1)/max(mean(Se,1)),'color','r');

break
%%

figure(5896)
xxx=1239.8./(Lscale_dfl(1,:).*1e9)-L0_dfl;
yyy=mean(Se_dfl,1);
Ph_energy=L0_dfl/6.241e18;
% c=3e8;
% omega=2*pi*c/xlamds;
bw=0.0002;
N_photons_avrg=mean(Ez(:,end))/Ph_energy;
n1=find(xxx<-700*bw/2,1,'last');
n2=find(xxx>700*bw/2,1,'first');
Sp_01=sum(yyy(n1:n2))/sum(yyy);
N_photons_01=N_photons_avrg*Sp_01
plot(xxx,yyy);
hold all
plot([xxx(n1) xxx(n1:n2) xxx(n2)],[0 yyy(n1:n2) 0]);
hold off
break
 %%
 cd('c:\Users\Svitozar\GoogleDrive\Documents\Dissertation\Figures\XFEL\S2E\');
%% General amplification process plots
Di=3;
outpot_e(1,d(Di))
outpot_ph(2,d(Di))
%% radiation longitudinal size and position

% figure
% plot(d(1).outp.Zscale,d(1).outp.power.std.*2)
% hold all
% plot(d(1).outp.Zscale,d(1).outp.power.peakpos)
% hold off

%% Bunch profile parameters at given position

Di=2;
Z=36; %[m]
outplot_z(3,d(Di),Z);

%% dfl file
Di=1;

[XX,N]=fieldimport_all([nm_p{Di},'.dfl'],d(Di).inp.ncar,1);
%XX=prop_TF(XX,d(Di).inp.dgrid*2,d(Di).inp.xlamds,0);
zscale=linspace(0,d(Di).inp.xlamds*d(Di).inp.zsep*N,N);
fieldplot3d(110,XX,d(Di).inp.dgrid*2,1,zscale,[nm_p{Di},'.dfl'],1);
break
%%
figure(67)
plot(reshape(sum(sum(abs(XX).^2,1),2),1,[]));
%%
waistscan_1(XX(:,:,200:700),d(Di).inp.dgrid*2,d(Di).inp.xlamds,[-8:0.5:3]);
%[-8:0.5:3]
%%
sigma=2e-5;
z=-5.5; %650
% sigma=1.4e-5;
% z=-3.5;
X=fieldgaussian(d(Di).inp.ncar,d(Di).inp.dgrid*2,sigma,sigma,z,z,d(Di).inp.xlamds,1);
waistscan_g(X,d(Di).inp.dgrid*2,d(Di).inp.xlamds,[-8:0.5:3]);
%% Comparison plots

% Z1=24;
% Z2=24;
% a1=1;
% a2=2;
% 
% figure(1432)
% plot(d(a1).outp.Sscale,d(a1).outp.e_spread.v(:,find(d(a1).outp.Zscale>=Z1,1,'first')),'linewidth',0.5,'linestyle','--','color','b');
% hold on
% plot(d(a2).outp.Sscale,d(a2).outp.e_spread.v(:,find(d(a2).outp.Zscale>=Z2,1,'first')),'linewidth',0.5,'linestyle','--','color','r');
% hold off
%%
if DiN==2
    
    Z=0;
    for Di=1:DiN
        Zi(Di)=find(d(Di).outp.Zscale>=Z,1,'first');
    end

    figure(56);  
    plot(d(1).outp.Lamdscale,d(1).outp.spectrum_mid.v(:,Zi(1)),'linewidth',2);%/max(d(1).outp.spectrum_mid.v(:,Zi(1)))
    hold all
    plot(d(2).outp.Lamdscale,d(2).outp.spectrum_mid.v(:,Zi(2)),'linewidth',2,'color','r','linestyle','--');%/max(d(2).outp.spectrum_mid.v(:,Zi(2)))
    hold off
    %legend('direct\newlineapproach','phenomenological\newlineapproach','Location','northwest')
    %xlim([1.2385 1.2395]);
    %xlim([1.764 1.766]);
    ylabel('P(\lambda) [arb. units]');
    xlabel('\lambda [nm]');
    
    

    figure(59);
    semilogy(d(1).outp.Zscale,d(1).outp.power.max_S,'LineWidth',2);
    hold all
    semilogy(d(2).outp.Zscale,d(2).outp.power.max_S,'LineWidth',2,'color','r','linestyle','--');
    hold off
    %     legend('direct approach','phenomenological approach','Location','northwest')
    l=legend(nm_p{1},nm_p{2},'Location','northwest');
    set(l, 'Interpreter', 'none');
    ylabel('P [W]');
    xlabel('z [m]');
    %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
    axis tight

    figure(591);
    semilogy(d(1).outp.Zscale,d(1).outp.power.E,'LineWidth',2);
    hold all
    semilogy(d(2).outp.Zscale,d(2).outp.power.E,'LineWidth',2,'color','r','linestyle','--');
    hold off
    %     legend('direct approach','phenomenological approach','Location','northwest')
    l=legend(nm_p{1},nm_p{2},'Location','northwest');
    set(l, 'Interpreter', 'none')
    ylabel('E [J]');
    xlabel('z [m]');
    %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
    axis tight
    clear l
end




%% spectrum and power evolution plots
Di=10;

if d(Di).inp.itdp~=1
    error('steady state: cannot plot spectrum & power')
end
    
    %plot pulse evolution
    Zgrid_p=repmat(d(Di).outp.Zscale,d(Di).outp.Sn,1);
    Sgrid=repmat(d(Di).outp.Sscale,1,d(Di).outp.Zn);
    Lamdgrid=repmat(d(Di).outp.Lamdscale,1,d(Di).outp.Zn);
    P=double(d(Di).outp.power.v);
    S=double(d(Di).outp.spectrum_mid.v);
%     figure(611);
%     h=surf(Zgrid,Sgrid,P,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
    Zgrid_s=repmat(d(Di).outp.Zscale,numel(d(Di).outp.Lamdscale),1);
     figure(612);
%     h=surf(Zgrid,Lamdgrid,S,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
     Zm=10; Sm=5;
     [Zgrid_p_n, Sgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).outp.Zn/Zm), linspace(min(d(Di).outp.Sscale),max(d(Di).outp.Sscale),d(Di).outp.Sn/Sm));
     P_n=interp2(Zgrid_p,Sgrid,P,Zgrid_p_n,Sgrid_n);
     %P_n(P_n<1e5)=NaN;
     H.h5=surf(Zgrid_p_n,Sgrid_n,P_n,log(P_n),'linestyle','none');
     %set(get(H.h3,'Parent'),'ZScale','lin');
     %set(get(H.h3,'Parent'),'ZTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
     title('power');
     axis tight
     xlabel('z [m]');
     ylabel('s [m]');
     zlabel('P [W]');
     
     Zm=10; Lm=2;
    [Zgrid_s_n, Lamdgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).outp.Zn/Zm), linspace(min(d(Di).outp.Lamdscale),max(d(Di).outp.Lamdscale),numel(d(Di).outp.Lamdscale)/Lm));
    S_n=interp2(Zgrid_s,Lamdgrid,S,Zgrid_s_n,Lamdgrid_n);
    S_n(:,1)=NaN;
    figure(622);
    H.h6=surf(Zgrid_s_n,Lamdgrid_n.*1e9,S_n,log(S_n),'linestyle','none');
    %set(get(H.h4,'Parent'),'ZScale','lin');
    %set(get(H.h4,'Parent'),'ZTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
    title('spectrum');
    axis tight
     xlabel('z [m]');
     ylabel('\lambda [nm]');
     zlabel('P(\lambda) [arb. units]');
    %xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
    
    clear Zm Lm Sm Zgrid_p Zgrid_s Sgrid Lamdgrid Zgrid_p_n Zgrid_s_n Sgrid_n Lamdgrid_n P P_n S S_n

    
    
    
    
    

%% Ebeam focusing plot
break
%     Xmin=min(d(Di).outp.x.mean_S-2*d(Di).outp.xrms.mean_S);
%     Xmax=max(d(Di).outp.x.mean_S+2*d(Di).outp.xrms.mean_S);
%     Ymin=min(d(Di).outp.y.mean_S-2*d(Di).outp.yrms.mean_S);
%     Ymax=max(d(Di).outp.y.mean_S+2*d(Di).outp.yrms.mean_S);
%    [X,Y]=meshgrid(linspace(Xmin,Xmax,50),linspace(Ymin,Ymax,50));
Zn=linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),30);

x_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.x.mean_S,Zn);
y_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.y.mean_S,Zn);
xrms_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.xrms.mean_S,Zn);
yrms_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S,Zn);

    figure(710);
    for i=1:100
        t=2*pi/100*i;
    x=x_mean_S+1*xrms_mean_S*cos(t);
    y=y_mean_S+1*yrms_mean_S*sin(t);
    plot3(x,y,Zn,'linewidth',1.5);
    hold all
    end
    hold off
    grid on
    axis tight
    
    clear x y x_mean_S y_mean_S xrms_mean_S yrms_mean_S Zn i t
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

%% tmp

X=d(Di).outp.Sscale*1e6;
E=d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0;
I=d(Di).outp.current;
dE=d(Di).outp.e_spread.v(:,Zi);
dE(E<9240)=NaN;
E(E<9240)=NaN;

figure(3456)
plot(X,E,'linewidth',1.5);
xlabel('s [um]');
ylabel('\gamma');
axis tight
xlim([min(X) max(X)]);
set(figure(3456), 'Position', [100, 100, 550, 450]);

figure(3457)
plot(X,dE,'linewidth',1.5);
xlabel('s [um]');
ylabel('\sigma_\gamma');
axis tight
xlim([min(X) max(X)]);
set(figure(3457), 'Position', [100, 100, 550, 450]);

figure(3458)
plot(X,I,'linewidth',1.5);
xlabel('s [um]');
ylabel('I [A]');
axis tight
xlim([min(X) max(X)]);
set(figure(3458), 'Position', [100, 100, 550, 450]);

%% tdp spectrum comparison
Z=60;
Zi=find(Zscale(1,:)<=Z,1,'last');

fig=figure(5257);

clf;

L0=929.3-0.08;
%L0=500;
%plot(1239.8./(Lscale.*1e9)-L0,Sef_all(:,:,Zi),'linewidth',1,'color',[1 0.8 0.8]);
%hold on
%plot(1239.8./(Lscale.*1e9)-L0,Sef_all(1,:,Zi),'linewidth',0.5,'color',[0.7 0.5 0.5]);
hl1=line(1239.8./(Lscale.*1e9)-L0,mean(Sefc_all(:,:,Zi),1),'linewidth',2,'color',[0.7 0.2 0.2]);
hold off
axis tight
ylabel('P(\lambda) [arb. units]');
xlabel('Photon energy [eV]');
%xlim([-3 3]);

hold on

%plot(1239.8./(Lscale.*1e9)-L0,Sefc_all(:,:,Zi),'linewidth',1,'color',[0.8 1 0.8]);
%hold on
%plot(1239.8./(Lscale.*1e9)-L0,Sefc_all(1,:,Zi),'linewidth',0.5,'color',[0.5 0.7 0.5]);
hl2=line(1239.8./(Lscale.*1e9)-L0,mean(Sefcnn_all(:,:,Zi),1),'linewidth',2,'color',[0.2 0.7 0.2]);
%hold off
axis tight
ylabel('P(\lambda) [arb. units]');
xlabel('Photon energy [eV]');
xlim([-3 3]);


ax1=gca;
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');

%figure(5258)
hl3=line(1239.8./(Lscale.*1e9)-L0,cumsum(mean(Sefc_all(:,:,Zi),1))./sum(mean(Sefc_all(:,:,Zi),1)),'linewidth',1,'color',[0.7 0.2 0.2],'linestyle','--','Parent',ax2);
hold all
hl4=line(1239.8./(Lscale.*1e9)-L0,cumsum(mean(Sefcnn_all(:,:,Zi),1))./sum(mean(Sefcnn_all(:,:,Zi),1)),'linewidth',1,'color',[0.2 0.7 0.2],'linestyle','--','Parent',ax2);
hold off
xlim([-3 3]);
