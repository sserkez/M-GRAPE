% reads .out file
fclose all;
clear all
%DiN=2;
%  nm_p{1}='tdp\500_u1_tdp.out'; %Data
%   nm_p{1}='tdp\1000_u1_tdp.out'; %Data
%  nm_p{1}='tdp\1000_u2_tdp_u8_DA.out'; %Data
%  nm_p{2}='tdp\1000_u2_tdp_u8_KK.out';
%nm_p{1}='C:\-D-\Work\SASE3_chicane\whole\sase3.psase.3.out';
% nm_p{1}='C:\-D-\Work\SASE3_chicane\run3\U1.1.out';
% nm_p{2}='C:\-D-\Work\SASE3_chicane\run3\U2.1.out';
%nm_p{1}='C:\-D-\Work\SASE3_chicane\U2.1_2.out';
%nm_p{2}='C:\-D-\Work\SASE3_chicane\U2.1_3.out';
%   nm_p{1}='C:\-D-\Work\SASE3_chicane\run2\U2.1_4u_3um_5u.out';
%   nm_p{2}='C:\-D-\Work\SASE3_chicane\run2\U3.1_4u_3um_5u.out';
  
%   nm_p{1}='C:\-D-\Work\SASE3_chicane\run2\U1.1_5u.out';
%    nm_p{1}='C:\-D-\Work\SASE3_chicane\run2\U2.1_5u_3um_5u.out';
%    nm_p{2}='C:\-D-\Work\SASE3_chicane\run2\U3.1_5u_3um_5u.out';
  
%        nm_p{1}='C:\-D-\Work\SASE3_chicane\run4\U1.1.out';
%        nm_p{2}='C:\-D-\Work\SASE3_chicane\run4\U2.1.out';
   
%      nm_p{1}='C:\-D-\Work\LCLS\New_matlab\tdp\1000_u1_tdp.out';
      nm_p{1}='C:\-D-\Work\SASE3_chicane\run5\U1.2.out';
      nm_p{2}='C:\-D-\Work\SASE3_chicane\run5\U2.2.out';
       
 DiN=size(nm_p,2);
% DiN=1;
for Di=1:DiN %Data index
   d(Di)=outread(nm_p{Di},1); 
end

%% General amplification process plots
Di=1;

fig1=figure(1);
clf
set(fig1,'name',['Electrons ',nm_p{Di}],'numbertitle','off');
H.h1(Di).h1=subplot(4,1,1);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.aw,d(Di).outp.Zscale,d(Di).outp.qfld);
set(haxes,{'ycolor'},{'b';'r'})
ylabel(haxes(1),'Undulator\newlinea_W');
ylabel(haxes(2),'Quadrupole\newlinedB/dx [T/m]');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2,'color','r');
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
set(haxes(2),'YLim',[min(d(Di).outp.qfld)*1.1 max(d(Di).outp.qfld)*1.1]);
set(haxes(1),'box','off');
clear ax i

%---------
try
H.h1(Di).h2=subplot(4,1,2);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.xrms.mean_S,d(Di).outp.Zscale,d(Di).outp.x.mean_S);
hline(3)=line(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S, 'Parent', haxes(1), 'linestyle','--','color','b','linewidth',2.5);
hline(4)=line(d(Di).outp.Zscale,d(Di).outp.y.mean_S, 'Parent', haxes(2), 'linestyle','--','color',[0 .5 0],'linewidth',1.5);
ylabel(haxes(1),'Beam size\newline\sigma [m]');
ylabel(haxes(2),'Beam centroid\newline-x - -y [m]');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',1.5,'color','b','linewidth',2.5);
set(hline(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',1.5);
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
set(haxes(1),'box','off');
max_xy_rms=max([d(Di).outp.yrms.mean_S d(Di).outp.xrms.mean_S]);
maxabs_xy=max([abs(d(Di).outp.y.mean_S) abs(d(Di).outp.x.mean_S)]);
 set(haxes(1),'YLim',[0 max_xy_rms*1.1],'YTickMode','auto');
 set(haxes(2),'YLim',[-maxabs_xy*1.1 maxabs_xy*1.1],'YTickMode','auto');

clear ax i max_xy_rms maxabs_xy
catch
    H.h1(Di).h2=subplot(4,1,2);
    plot(d(Di).outp.Zscale,d(Di).outp.xrms.mean_S,'linewidth',1.5)
    hold all
    plot(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S,'linewidth',1.5);
    %hline(3)=line(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S, 'Parent', haxes(1), 'linestyle','--','color','b','linewidth',2.5);
    %hline(4)=line(d(Di).outp.Zscale,d(Di).outp.y.mean_S, 'Parent', haxes(2), 'linestyle','--','color',[0 .5 0],'linewidth',1.5);
    ylabel('Beam size\newline\sigma [m]');
    %ylabel(haxes(2),'Beam centroid\newline-x - -y [m]');
    %xlabel(haxes(2),'z [m]');
    %set(hline(1),'LineWidth',1.5,'color','b','linewidth',2.5);
    %set(hline(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',1.5);
    xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
%    
%     set(haxes(1),'box','off');
    max_xy_rms=max([d(Di).outp.yrms.mean_S d(Di).outp.xrms.mean_S]);
%     maxabs_xy=max([abs(d(Di).outp.y.mean_S) abs(d(Di).outp.x.mean_S)]);
    ylim([0 max_xy_rms*1.1]);
%      set(haxes(2),'YLim',[-maxabs_xy*1.1 maxabs_xy*1.1],'YTickMode','auto');

    clear ax i max_xy_rms maxabs_xy
end

%---------
H.h1(Di).h3=subplot(4,1,3);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.energy.mean_S+d(Di).inp.gamma0,d(Di).outp.Zscale,d(Di).outp.e_spread.mean_S);
ylabel(haxes(1),'energy (\Delta\gamma)');
ylabel(haxes(2),'spread (\sigma_\gamma)');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2);
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
set(haxes(1),'box','off');
% text(0,1,sprintf(' gamma_0 = %.0f', d(Di).inp.gamma0),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized');
clear ax i
%---------
H.h1(Di).h4=subplot(4,1,4);
plot(d(Di).outp.Zscale,d(Di).outp.bunching.max_S,'LineWidth',2,'color','r');
hold on
plot(d(Di).outp.Zscale,d(Di).outp.bunching.mean_S,'LineWidth',1,'color','k','linestyle','--');
ylabel('bunching\newline|<exp(i \theta)>|');
xlabel('z [m]');
xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
    legend('max','mean','location','NorthWest');





fig2=figure(2);
clf;
set(fig2,'name',['Photons ',nm_p{Di}],'numbertitle','off');

H.h2(Di).h1=subplot(3,1,1);
semilogy(d(Di).outp.Zscale,d(Di).outp.power.max_S,'LineWidth',2,'color','r');
hold all
semilogy(d(Di).outp.Zscale,d(Di).outp.power.mean_S,'LineWidth',1,'linestyle','--','color','k');
hold off
ylim([min(d(Di).outp.power.mean_S) max(d(Di).outp.power.max_S)]);
grid on
grid minor
set(gca,'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15],'YMinorTick','on','XMinorGrid','off');
ylabel('P [W]');
%xlabel('z [m]');
xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
legend('max','mean','location','NorthWest');

%YMinorTick on

H.h2(Di).h2=subplot(3,1,2);
if d(Di).inp.itdp==1
   [haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.spectrum_mid.max_S,d(Di).outp.Zscale,2.*d(Di).outp.spectrum_mid.std_lamd./d(Di).inp.xlamds,'semilogy','plot');

    %semilogy(d(Di).outp.Zscale,d(Di).outp.spectrum_mid.max_S,'LineWidth',2);


    %plot(d(Di).outp.Zscale,d(Di).outp.spectrum_mid.std_lamd,'LineWidth',2);
    %hold all
    ylabel(haxes(1),'P_{max}(\lambda) [a.u.]') % label left y-axis
    ylabel(haxes(2),'2\sigma_{P(\lambda)}/\lambda') % label right y-axis
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
    set(haxes(1),'box','off');
    set(haxes(2),'YTickMode','auto');
    text(0,1,sprintf(' xlamds = %.3e', d(Di).inp.xlamds),...
                        'HorizontalAlignment','left','VerticalAlignment',...
                        'top','FontSize',10,'units','normalized');
    text(0.5,1,sprintf('(on axis)'),...
                        'HorizontalAlignment','center','VerticalAlignment',...
                        'top','FontSize',8,'units','normalized');
else
    
    text(0.5,0.5,sprintf('steadystate (no data)'),...
                        'HorizontalAlignment','center','VerticalAlignment',...
                        'middle','FontSize',8,'units','normalized');
end
                
                

H.h2(Di).h3=subplot(3,1,3);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Zscale,d(Di).outp.r_size.mean_S.*2*1e6,d(Di).outp.Zscale,d(Di).outp.power.std.*2*1e6);
ylabel(haxes(1),'2\sigma_{rad}^{transverse} [\mum]');
ylabel(haxes(2),'2\sigma_{rad}^{longitudinal} [\mum]');
xlabel(haxes(2),'z [m]');
set(haxes,{'ycolor'},{'b';'r'})
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2,'color','r');
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
end
set(haxes(1),'box','off');

clear haxes hline i fig1 fig2
%%

figure
plot(d(1).outp.Zscale,d(1).outp.power.std.*2)
hold all
plot(d(1).outp.Zscale,d(1).outp.power.peakpos)
hold off

%% Bunch profile parameters at given position

Di=1;
Z=0; %[m]

%  Di=1;
%  Z=0; %[m]

if Z>max(d(Di).outp.Zscale)
    disp('Z parameter exceeds data limits');
    Z=max(d(Di).outp.Zscale);
elseif Z<min(d(Di).outp.Zscale)
    Z=min(d(Di).outp.Zscale);
    disp('Z parameter exceeds data limits');
end
Zi=find(d(Di).outp.Zscale>=Z,1,'first');
Z=d(Di).outp.Zscale(Zi);
minsc=min(d(Di).outp.Sscale);
maxsc=max(d(Di).outp.Sscale);

fig3=figure(3);
set(fig3,'name',['Bunch ',nm_p{Di}],'numbertitle','off');
H.h3(Di).h1=subplot(4,1,1);
[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Sscale,d(Di).outp.current,d(Di).outp.Sscale,d(Di).outp.power.v(:,Zi));
ylabel(haxes(1),'I [A]');
ylabel(haxes(2),'P [W]');
xlabel(haxes(2),'s [m]');
set(haxes,{'ycolor'},{'k';[0 0.5 0]})
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',0.5,'color','k','linestyle','--');
set(hline(2),'LineWidth',1.5,'color',[0 0.5 0]);
for i=1:2
    set(haxes(i),'XLim',[d(Di).outp.Sscale(1) d(Di).outp.Sscale(end)]);
end
set(haxes(1),'box','off');
title(['Z=',num2str(Z),'m']);
text(0,1,sprintf(' q ~ %.3f nC', d(Di).inp.charge*1e9),...
                        'HorizontalAlignment','left','VerticalAlignment',...
                        'top','FontSize',10,'units','normalized');
text(1,1,sprintf('E = %.3e J ', d(Di).outp.power.E(Zi)),...
                        'HorizontalAlignment','right','VerticalAlignment',...
                        'top','FontSize',10,'units','normalized');




% plot(d(Di).outp.Sscale,d(Di).outp.current,'linewidth',2,'color','k');
% %xlabel('s [m]');
% ylabel('I [A]');
% axis tight
% %xlim([minsc maxsc]);
% title(['Z=',num2str(Z),'m']);
% text(0,1,sprintf(' q ~ %.3f nC', d(Di).inp.charge*1e9),...
%                         'HorizontalAlignment','left','VerticalAlignment',...
%                         'top','FontSize',10,'units','normalized');


                    
% H.h3(Di).h2=subplot(4,1,2);
% plot(d(Di).outp.Sscale,d(Di).outp.power.v(:,Zi),'linewidth',1.5,'color',[0 0.5 0]);
% %xlabel('s [m]');
% ylabel('P [W]');
% axis tight
% 
% text(0,1,sprintf(' E = %.3e J', d(Di).outp.power.E(Zi)),...
%                         'HorizontalAlignment','left','VerticalAlignment',...
%                         'top','FontSize',10,'units','normalized');

H.h3(Di).h3=subplot(4,1,2);
%[haxes,hline(1),hline(2)] = plotyy(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0,d(Di).outp.Sscale,d(Di).outp.bunching.v(:,Zi));
plot(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0,'linewidth',1.5);
hold on
% hline(3)=line(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0-d(Di).outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r','parent',haxes(1));
% hline(4)=line(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0+d(Di).outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r','parent',haxes(1));
plot(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0-d(Di).outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r');
plot(d(Di).outp.Sscale,d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0+d(Di).outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r');  hold off
%xlabel('s [m]');
ylabel('\gamma');
%set(haxes(1),'box','off');
axis tight
legend('<\gamma>','+/- \sigma_\gamma')
%set(gca,'XTickLabel',[]);

H.h3(Di).h4=subplot(4,1,3);
plot(d(Di).outp.Sscale,d(Di).outp.bunching.v(:,Zi));
axis tight
ylabel('bunching\newline|<exp(i \theta)>|');
xlabel('s [m]');


H.h3(Di).h4=subplot(4,1,4);
plot(d(Di).outp.Lamdscale*1e9,d(Di).outp.spectrum_mid.v(:,Zi),'linewidth',1.5,'color','r');
axis tight
ylabel('P(\lambda) [a.u.]');
xlabel('\lambda [nm]');
% figure(55555);
% subplot(2,1,1)
% plot(d(Di).outp.Sscale,d(Di).outp.power.v(:,Zi),'linewidth',1.5)
% subplot(2,1,2)
% [haxes,hline(1),hline(2)]=plotyy(d(Di).outp.Sscale,unwrap(d(Di).outp.phi_mid.v(:,Zi)),d(Di).outp.Sscale,[0; diff(unwrap(d(Di).outp.phi_mid.v(:,Zi)))]);
% hline(3)=line(d(Di).outp.Sscale,linspace(0,0,d(Di).outp.Sn),'parent',haxes(2),'color','k','linestyle','--');

clear fig3 haxes hline
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
    ylabel('P(\lambda) [a.u.]');
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
Di=1;

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
     zlabel('P(\lambda) [a.u.]');
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