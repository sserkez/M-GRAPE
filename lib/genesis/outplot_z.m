function H=outplot_z(figN,d,Z)

Nplots=4;

if Z>max(d.outp.Zscale)
    disp('Z parameter exceeds data limits');
    Z=max(d.outp.Zscale);
elseif Z<min(d.outp.Zscale)
    Z=min(d.outp.Zscale);
    disp('Z parameter exceeds data limits');
end
Zi=find(d.outp.Zscale>=Z,1,'first');
Z=d.outp.Zscale(Zi);
minsc=min(d.outp.Sscale);
maxsc=max(d.outp.Sscale);

fig3=figure(figN);
set(fig3,'name',['Bunch ',d.nm_p],'numbertitle','off');



H.h3.h1=subplot(Nplots,1,1);
[haxes1,hline(1),hline(2)] = plotyy(d.outp.Sscale,d.outp.current,d.outp.Sscale,d.outp.power.v(:,Zi));
ylabel(haxes1(1),'I [A]');
ylabel(haxes1(2),'P [W]');
xlabel(haxes1(2),'s [m]');
set(haxes1,{'ycolor'},{'k';[0 0.5 0]})
%xlabel(haxes(2),'z [m]');
set(hline(1),'LineWidth',0.5,'color','k','linestyle','--');
set(hline(2),'LineWidth',1.5,'color',[0 0.5 0]);
for i=1:2
    set(haxes1(i),'XLim',[minsc maxsc]);
end
set(haxes1(1),'box','off','Ylim',[0 max(d.outp.current)]);
set(haxes1(2),'Ylim',[0 max(d.outp.power.v(:,Zi))]);
title(['Z=',num2str(Z),'m']);
text(0,1,sprintf(' q ~ %.3f nC', d.inp.charge*1e9),...
                        'HorizontalAlignment','left','VerticalAlignment',...
                        'top','FontSize',10,'units','normalized');
text(1,1,sprintf('E = %.3e J ', d.outp.power.E(Zi)),...
                        'HorizontalAlignment','right','VerticalAlignment',...
                        'top','FontSize',10,'units','normalized');




% plot(d.outp.Sscale,d.outp.current,'linewidth',2,'color','k');
% %xlabel('s [m]');
% ylabel('I [A]');
% axis tight
% %xlim([minsc maxsc]);
% title(['Z=',num2str(Z),'m']);
% text(0,1,sprintf(' q ~ %.3f nC', d.inp.charge*1e9),...
%                         'HorizontalAlignment','left','VerticalAlignment',...
%                         'top','FontSize',10,'units','normalized');


                    
% H.h3.h2=subplot(4,1,2);
% plot(d.outp.Sscale,d.outp.power.v(:,Zi),'linewidth',1.5,'color',[0 0.5 0]);
% %xlabel('s [m]');
% ylabel('P [W]');
% axis tight
% 
% text(0,1,sprintf(' E = %.3e J', d.outp.power.E(Zi)),...
%                         'HorizontalAlignment','left','VerticalAlignment',...
%                         'top','FontSize',10,'units','normalized');

H.h3.h2=subplot(Nplots,1,2);
%[haxes,hline(1),hline(2)] = plotyy(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0,d.outp.Sscale,d.outp.bunching.v(:,Zi));
haxes2=plot(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0,'linewidth',1.5);
hold on
% hline(3)=line(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0-d.outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r','parent',haxes(1));
% hline(4)=line(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0+d.outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r','parent',haxes(1));
plot(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0-d.outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r');
plot(d.outp.Sscale,d.outp.energy.v(:,Zi)+d.inp.gamma0+d.outp.e_spread.v(:,Zi),'linewidth',0.5,'linestyle','--','color','r');  hold off
%xlabel('s [m]');
ylabel('\gamma');
%set(haxes(1),'box','off');
axis tight
legend('<\gamma>','+/- \sigma_\gamma')
%set(gca,'XTickLabel',[]);

% H.h3.h3=subplot(Nplots,1,3);
% plot(d.outp.Sscale,d.outp.bunching.v(:,Zi));
% axis tight
% ylabel('bunching\newline|<exp(i \theta)>|');
% %xlabel('s [m]');


H.h3.h5=subplot(Nplots,1,4);
haxes3=plot(d.outp.Lamdscale*1e9,d.outp.spectrum_mid.v(:,Zi),'linewidth',1.5,'color','r');
axis tight
ylabel('P(\lambda) [a.u.]');
xlabel('\lambda [nm]');

text(1,1,sprintf('(on axis)  '), ...
                        'HorizontalAlignment','right','VerticalAlignment',...
                        'top','FontSize',8,'units','normalized');

%                     
% H.h3.h4=subplot(5,1,4);
% [haxes,hline(1),hline(2)]=plotyy(d.outp.Sscale,unwrap(d.outp.phi_mid.v(:,Zi)),d.outp.Sscale,[0; diff(unwrap(d.outp.phi_mid.v(:,Zi)))]);
% hline(3)=line(d.outp.Sscale,linspace(0,0,d.outp.Sn),'parent',haxes(2),'color','k','linestyle','--');
% for i=1:2
%     set(haxes(i),'XLim',[minsc maxsc]);
% end
% set(haxes(2),'YLim',[-pi/4 pi/4]);
% ylabel(haxes(1),'\Phi [rad]');
% ylabel(haxes(2),'\Delta\Phi [rad]');
% xlabel(haxes(2),'s [m]');

H.h3.h4=subplot(Nplots,1,3);
haxes4=plot(d.outp.Sscale,[0; diff(unwrap(d.outp.phi_mid.v(:,Zi)))],'color','k');
hold on
line(d.outp.Sscale,linspace(0,0,d.outp.Sn),'color','k','linestyle','--');
hold off

xlim([minsc maxsc]);
ylim([-pi/4 pi/4]);
%ylabel(haxes(1),'\Phi [rad]');
ylabel('\Delta\Phi [rad]');
xlabel('s [m]');
    text(0.5,1,sprintf('(on axis)'),...
                        'HorizontalAlignment','center','VerticalAlignment',...
                        'top','FontSize',8,'units','normalized');

%H=linkprop([haxes1 haxes2 haxes3 haxes4], 'XLim');
%H=linkaxes([H.h3.h1 H.h3.h2 H.h3.h5 H.h3.h4], 'x');
clear fig3 haxes hline
return