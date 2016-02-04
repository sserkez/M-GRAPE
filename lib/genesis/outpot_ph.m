function Handle=outpot_ph(figN,d)
%figN=550;
fig2=figure(figN);
clf
set(fig2,'name',['Photons ',d.nm_p],'numbertitle','off');

 H.h2.h1=subplot(3,1,1);
 [haxes1,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.power.max_S,d.outp.Zscale,d.outp.power.E,'semilogy','semilogy');

axis tight
set(haxes1(1),'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15],'YMinorTick','on','XMinorGrid','off');
set(haxes1(2),'YTick',[1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1],'YMinorTick','on','XMinorGrid','off');
set(haxes1(1), 'YMinorTick','on');
set(haxes1(2), 'YMinorTick','on');
for i=1:2
    set(haxes1(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
        %set(haxes(i),'axis','tight');
    axis tight
end

if min(d.outp.power.E)~=max(d.outp.power.E)
    set(haxes1(2),'YLim',[min(d.outp.power.E) max(d.outp.power.E)]);
end
    set(haxes1(1),'box','off');
    set(hline(1),'LineWidth',2,'color','r');
    set(hline(2),'LineWidth',2,'linestyle','--','color','b');
    set(haxes1,{'ycolor'},{'r';'b'})
    %set(haxes(2),'YTickMode','auto');
    text(1,0,sprintf(' P_{max}^{end} = %.2e W \n E_{pulse}^{end} = %.2e J ', d.outp.power.max_S(end), d.outp.power.E(end)),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'bottom','FontSize',10,'units','normalized');
    ylabel(haxes1(1),'P_{max} [W]') % label left y-axis
    ylabel(haxes1(2),'E_{pulse} [J]') % label right y-axis
% semilogy(d.outp.Zscale,d.outp.power.max_S,'LineWidth',2,'color','r');
% hold all
% semilogy(d.outp.Zscale,d.outp.power.mean_S,'LineWidth',1,'linestyle','--','color','k');
% hold off
% ylim([min(d.outp.power.mean_S) max(d.outp.power.max_S)]);
% grid on
% grid minor
% set(gsa,'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15],'YMinorTick','on','XMinorGrid','off');
% ylabel('P [W]');
% %xlabel('z [m]');
% xlim([d.outp.Zscale(1) d.outp.Zscale(end)]);
% legend('max','mean','location','NorthWest');

%YMinorTick on

H.h2.h2=subplot(3,1,2);
if d.inp.itdp==1
%    [haxes,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.spectrum_mid.max_S,d.outp.Zscale,2.*d.outp.spectrum_mid.std_lamd./d.inp.xlamds,'semilogy','plot');
% 
%     %semilogy(d.outp.Zscale,d.outp.spectrum_mid.max_S,'LineWidth',2);
% 
% 
%     %plot(d.outp.Zscale,d.outp.spectrum_mid.std_lamd,'LineWidth',2);
%     %hold all
%     ylabel(haxes(1),'P_{max}(\lambda) [a.u.]') % label left y-axis
%     ylabel(haxes(2),'2\sigma_{P(\lambda)}/\lambda') % label right y-axis
%     set(hline(1),'LineWidth',1.5,'color','b','linewidth',2);
%     set(hline(2),'LineWidth',1.5,'color',[0 0.5 0],'linewidth',2);
%     %xlabel(haxes(2),'z [m]') % label x-axis
%     axis tight
%     set(haxes(1),'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
%     for i=1:2
%         set(haxes(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
%             %set(haxes(i),'axis','tight');
%             axis tight
%     end

% [haxes,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.spectrum_mid.max_S,d.outp.Zscale,2.*d.outp.spectrum_mid.std_lamd2./d.inp.xlamds,'semilogy','plot');

[haxes2,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.spectrum_mid.max_S,d.outp.Zscale,2.*d.outp.spectrum_mid.std_lamd2./d.inp.xlamds,'semilogy','plot');
if isfield(d.outp.spectrum_mid,'s_gfit')
    hline(3)=line(d.outp.Zscale,2.*d.outp.spectrum_mid.s_gfit./d.inp.xlamds, 'Parent', haxes2(2), 'linestyle','-','color',[0 .5 0],'linewidth',1.5);
'should_plot'
end

    ylabel(haxes2(1),'P_{max}(\lambda) [a.u.]') % label left y-axis
    ylabel(haxes2(2),'2\sigma_{P(\lambda)}/\lambda [nm]\newline-fit - -calc') % label right y-axis
    set(hline(1),'LineWidth',1.5,'color','b','linewidth',2);
    set(hline(2),'LineWidth',1.5,'color',[0 0.5 0],'linewidth',2,'linestyle','--');
    %xlabel(haxes(2),'z [m]') % label x-axis
    axis tight
    set(haxes2(1),'YTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
    for i=1:2
        set(haxes2(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
            %set(haxes(i),'axis','tight');
            axis tight
    end
    
    set(haxes2(1),'box','off');
    set(haxes2(2),'YTickMode','auto');
    text(0,1,sprintf(' xlamds = %.3f nm ', d.inp.xlamds*1e9),...
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
                
                

H.h2.h3=subplot(3,1,3);
[haxes3,hline(1),hline(2)] = plotyy(d.outp.Zscale,d.outp.r_size.mean_S_norm.*2*1e6,d.outp.Zscale,d.outp.power.std.*2*1e6);
ylabel(haxes3(1),'2\sigma_{rad}^{transverse} [\mum]');
ylabel(haxes3(2),'2\sigma_{rad}^{longitudinal} [\mum]');
xlabel(haxes3(2),'z [m]');
set(haxes3,{'ycolor'},{'b';'r'})
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',2);
set(hline(2),'LineWidth',2,'color','r');
for i=1:2
    axis tight
    set(haxes3(i),'XLim',[d.outp.Zscale(1) d.outp.Zscale(end)]);
end
set(haxes3(1),'box','off');

if d.inp.itdp
    Handle=linkprop([haxes1(1) haxes1(2) haxes2(1) haxes2(2) haxes3(1) haxes3(2)], 'XLim');
else
    Handle=linkprop([haxes1(1) haxes1(2) haxes3(1) haxes3(2)], 'XLim');
end
clear hline i fig1 fig2
return