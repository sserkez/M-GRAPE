function Handle=outpot_e(figN,d)
%figN=550;
fig1=figure(figN);
clf
set(fig1,'name',['Electrons ',d.nm_p],'numbertitle','off');
H.h1.h1=subplot(4,1,1);
[haxes1,hline1(1),hline1(2)] = plotyy(d.outp.Zscale,d.outp.aw,d.outp.Zscale,d.outp.qfld);
set(haxes1,{'ycolor'},{'b';'r'})
ylabel(haxes1(1),'Undulator\newlinea_W');
ylabel(haxes1(2),'Quadrupole\newlinedB/dx [T/m]');
%xlabel(haxes(2),'z [m]');
set(hline1(1),'LineWidth',2);
set(hline1(2),'LineWidth',2,'color','r');
for i=1:2
    set(haxes1(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
end
if min(d.outp.qfld)~=max(d.outp.qfld)
    set(haxes1(2),'YLim',[min(d.outp.qfld)*1.1 max(d.outp.qfld)*1.1]);
end

Ku=d.outp.aw(d.outp.aw~=0);
Ku=Ku-max(Ku);
if sum(Ku)==0
    set(haxes1(1),'box','off');
else
    set(haxes1(1),'box','off','Ylim',[min(d.outp.aw(d.outp.aw~=0)) max(d.outp.aw)]);
    %set(haxes1(1), 'XTickMode', 'auto', 'XTickLabelMode', 'auto');
    set(haxes1(1),'YTick',0:10^floor(log10(max(d.outp.aw)-min(d.outp.aw(d.outp.aw~=0)))):max(d.outp.aw),'YMinorTick','on','XMinorGrid','off');
end

clear ax i
%linkprop([hline1(1) hline1(2)],'Xlim');

%---------
try
H.h1.h2=subplot(4,1,2);
[haxes2,hline2(1),hline2(2)] = plotyy(d.outp.Zscale,d.outp.xrms.mean_S,d.outp.Zscale,d.outp.x.mean_S);
hline2(3)=line(d.outp.Zscale,d.outp.yrms.mean_S, 'Parent', haxes2(1), 'linestyle','--','color','b','linewidth',2.5);
hline2(4)=line(d.outp.Zscale,d.outp.y.mean_S, 'Parent', haxes2(2), 'linestyle','--','color',[0 .5 0],'linewidth',1.5);
ylabel(haxes2(1),'Beam size\newline\sigma [m]');
ylabel(haxes2(2),'Beam centroid\newline-x - -y [m]');
%xlabel(haxes(2),'z [m]');
set(hline2(1),'LineWidth',1.5,'color','b','linewidth',2.5);
set(hline2(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',1.5);
for i=1:2
    set(haxes2(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
end
set(haxes2(1),'box','off');
max_xy_rms=max([d.outp.yrms.mean_S d.outp.xrms.mean_S]);
maxabs_xy=max([abs(d.outp.y.mean_S) abs(d.outp.x.mean_S)]);
 set(haxes2(1),'YLim',[0 max_xy_rms*1.1],'YTickMode','auto');
 set(haxes2(2),'YLim',[-maxabs_xy*1.1 maxabs_xy*1.1],'YTickMode','auto');

clear ax i max_xy_rms maxabs_xy
%hlinkx=linkprop([haxes1(1) haxes1(2) haxes2(1) haxes2(2)],'XLim');
catch
    H.h1.h2=subplot(4,1,2);
    plot(d.outp.Zscale,d.outp.xrms.mean_S,'linewidth',1.5)
    hold all
    haxes2=plot(d.outp.Zscale,d.outp.yrms.mean_S,'linewidth',1.5);
    %hline(3)=line(d.outp.Zscale,d.outp.yrms.mean_S, 'Parent', haxes(1), 'linestyle','--','color','b','linewidth',2.5);
    %hline(4)=line(d.outp.Zscale,d.outp.y.mean_S, 'Parent', haxes(2), 'linestyle','--','color',[0 .5 0],'linewidth',1.5);
    ylabel('Beam size\newline\sigma [m]');
    %ylabel(haxes(2),'Beam centroid\newline-x - -y [m]');
    %xlabel(haxes(2),'z [m]');
    %set(hline(1),'LineWidth',1.5,'color','b','linewidth',2.5);
    %set(hline(2),'LineWidth',1.5,'color',[0 .5 0],'linewidth',1.5);
    xlim([d.outp.Zscale(1) d.outp.Zscale(end)]);
%    
%     set(haxes(1),'box','off');
    max_xy_rms=max([d.outp.yrms.mean_S d.outp.xrms.mean_S]);
%     maxabs_xy=max([abs(d.outp.y.mean_S) abs(d.outp.x.mean_S)]);
    ylim([0 max_xy_rms*1.1]);
%      set(haxes(2),'YLim',[-maxabs_xy*1.1 maxabs_xy*1.1],'YTickMode','auto');
%hlinkx=linkprop([haxes1(1) haxes1(2) haxes2],'XLim');
    clear ax i max_xy_rms maxabs_xy
end

%---------
H.h1.h3=subplot(4,1,3);
[haxes3,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.energy.mean_S+d.inp.gamma0,d.outp.Zscale,d.outp.e_spread.mean_S);
ylabel(haxes3(1),'energy (\Delta\gamma)');
ylabel(haxes3(2),'spread (\sigma_\gamma)');
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2);
for i=1:2
    set(haxes3(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
end
set(haxes3(1),'box','off');
% text(0,1,sprintf(' gamma_0 = %.0f', d.inp.gamma0),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized');
clear ax i
%---------
H.h1.h4=subplot(4,1,4);
plot(d.outp.Zscale,d.outp.bunching.max_S,'LineWidth',2,'color','r');
hold on
plot(d.outp.Zscale,d.outp.bunching.mean_S,'LineWidth',1,'color','k','linestyle','--');
ylabel('bunching\newline|<exp(i \theta)>|');
xlabel('z [m]');
xlim([d.outp.Zscale(1) d.outp.Zscale(end)]);
legend('max','mean','location','NorthWest');

if isfield(d.outp,'x') && isfield(d.outp,'y')
%Handle=linkprop([H.h1.h1 H.h1.h2 H.h1.h3 H.h1.h4], 'XLim');
Handle=linkprop([haxes1(1) haxes1(2) haxes2(1) haxes2(2) haxes3(1) haxes3(2) H.h1.h4], 'XLim');
else  
Handle=linkprop([haxes1(1) haxes1(2) haxes2(1) haxes3(1) haxes3(2) H.h1.h4], 'XLim');
end
return