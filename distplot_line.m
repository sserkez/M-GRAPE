%function distplot(fignum,dist,index_range)
clear all
%dist=distread('C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_mm\slotted_mw.dist');
dist=distread('D:\Work\!PROJECTS\Phase_controlled_harmonics\SASE3_v3\beam_0.1nC_sase1_12kev_fresh1');
%%
%dist='c:\-D-\Work\LCLS\tmp\dist_1000.dat';
fignum=2631;

lengthscale=1; %if '0' - time scale, if '0' - length scale
umscale=1; % if '1' lengths 'x' and 'y' displayed in [um], is '0' - in [m]

if exist('index_range','var')
    
    X=dist.X(index_range);
    Y=dist.Y(index_range);
    PX=dist.PX(index_range);
    PY=dist.PY(index_range);
    T=dist.T(index_range);
    G=dist.G(index_range);
else
    X=dist.X;
    Y=dist.Y;
    PX=dist.PX;
    PY=dist.PY;
    T=dist.T;
    G=dist.G;
end
    size_dist=numel(T);
 
% if T(1)<T(end)
%     timedir=1;
% elseif T(1)>T(end)
%     timedir=-1;
% else
% end

N=100;
NT=300;
%ny=200;
T_0=linspace(min(T),max(T),NT);
if lengthscale==1
    T=T.*3e8;
else
    T=T.*1e15;
end

if umscale==1
    X=X.*1e6;
    Y=Y.*1e6;
end


T_=linspace(min(T),max(T),NT);
X_=linspace(min(X),max(X),N);
PX_=linspace(min(PX),max(PX),N);
Y_=linspace(min(Y),max(Y),N);
PY_=linspace(min(PY),max(PY),N);
G_=linspace(min(G),max(G),N);

It=histc(T,T_);
It=It./sum(It.*(T_0(2)-T_0(1))).*dist.charge;
% %Xt = hist3([T X],[NT N]); 
% Xt = hist3([X T],[N NT]); 
% Yt = hist3([Y T],[N NT]);
% % PXt = hist3([PX T],[N NT]); 
% % PYt = hist3([PY T],[N NT]);
% Gt = hist3([G T],[N NT]);
% Yx = hist3([Y X],[N N]); 
% % PXx = hist3([PX X],[N N]); 
% % PYy = hist3([PY Y],[N N]);
% Xpx = hist3([X PX],[N N]); 
% Ypy = hist3([Y PY],[N N]);
% PYpx = hist3([PY PX],[N N]);


for i=1:NT-1
    index=find(T>=T_(i) & T<T_(i+1));

    meanG(i,:)=mean(G(index));
    stdG(i,:)=std(G(index));
    emitx(i,:)=sqrt(mean(X(index).^2).*mean(PX(index).^2)-mean(X(index).*PX(index)).^2).*mean(G(index));
    emity(i,:)=sqrt(mean(Y(index).^2).*mean(PY(index).^2)-mean(Y(index).*PY(index)).^2).*mean(G(index));
    %stdG(i,:)=sqrt(meanG(i,:).^2-(G(index).^2-meanG(i,:)).^2);
    
end
    emitx=real(emitx);
    emity=real(emity);
    


% meanG=Gt'*G_'./sum(Gt,1)';
% stdG=Gt-repmat(meanG',size(Gt,1),[]);
% stdG=stdG'*G_';
%%
fig=figure(fignum);
clf;
%set(fig, 'Position', [100, 100, 600, 800]);
set(fig, 'Position', [100, 100, 1300, 300]);


subplot(1,3,1)
plot(T_,It,'linewidth',2,'color','k');
grid on;
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
ylabel('I[A]');
xlim([min(T_) max(T_)]);
text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','k');

subplot(1,3,2)
%plot(T_,[(meanG-1).*511000 ./1e9 ; NaN],'linewidth',2,'color','b');
[h,h1,h2]=plotyy(T_,[(meanG-1).*511000 ./1e9 ; NaN],T_,[(stdG-1).*511000 ./1e6 ; NaN]);
grid on;
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end

set(h1,'LineWidth',2,'color','b');
set(h2,'LineWidth',2,'color','m','linestyle','--');

set(h,{'ycolor'},{'b';'m'})

ylabel(h(1),'<E> [GeV]');
ylabel(h(2),'\sigma_E [MeV]');
set(h(1),'XLim',[min(T_) max(T_)]);
set(h(2),'XLim',[min(T_) max(T_)]);
%text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
%                    'HorizontalAlignment','left','VerticalAlignment',...
%                    'top','FontSize',10,'units','normalized','color','k');

% subplot(2,2,4)
% plot(T_,[(stdG-1).*511000 ./1e6 ; NaN],'linewidth',2,'color','m');
% grid on;
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
% ylabel('\sigma_E [MeV]');
% xlim([min(T_) max(T_)]);
% %text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
% %                    'HorizontalAlignment','left','VerticalAlignment',...
% %                    'top','FontSize',10,'units','normalized','color','k');

subplot(1,3,3)
plot(T_,[emitx ; NaN],'linewidth',2,'color',[0 0.5 0]);
hold on
plot(T_,[emity ; NaN],'linewidth',2,'color','r','linestyle','--');
hold off
grid on;
legend('x','y')
if lengthscale
    xlabel('s [m]');
else
    xlabel('t [fs]');
end
ylabel('\epsilon [\mum]');
xlim([min(T_) max(T_)]);
%text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
%                    'HorizontalAlignment','left','VerticalAlignment',...
%                    'top','FontSize',10,'units','normalized','color','k');




return











% subplot(2,4,1)
% %[haxes,hline1,hline1] = plotyy([],[],[],[]);
% %break
% imagesc(T_,G_,Gt);
% 
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
% ylabel('\gamma');
% set(gca,'YDir','normal')
% text(0,1,sprintf('  head'),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');
% text(1,1,sprintf('tail  '),...
%                     'HorizontalAlignment','right','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');
% axis tight
% % hold on
% % axes;
% % plot(T_,It);
% % hold off
% 
% subplot(2,4,4)
% imagesc(X_,Y_,Yx)
% if umscale
%     xlabel('x [um]');
%     ylabel('y [um]');
% else
%     xlabel('x [m]');
%     ylabel('y [m]');
% end
% %axis equal
% set(gca,'YDir','normal')
% 
% subplot(2,4,8)
% imagesc(PX_,PY_,PYpx)
%     xlabel('x''');
%     ylabel('y''');
% %axis equal
% set(gca,'YDir','normal')
%                 
% 
% subplot(2,4,2)
% imagesc(T_,X_,Xt)
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
% if umscale
%     ylabel('x [um]');
% else
%     ylabel('x [m]');
% end
% set(gca,'YDir','normal')
% text(0,1,sprintf('  head'),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');
% text(1,1,sprintf('tail  '),...
%                     'HorizontalAlignment','right','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');
%                 
% subplot(2,4,6)
% imagesc(T_,Y_,Yt)
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
% 
% if umscale
%     ylabel('y [um]');
% else
%     ylabel('y [m]');
% end
% 
% set(gca,'YDir','normal')
% text(0,1,sprintf('  head'),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');
% text(1,1,sprintf('tail  '),...
%                     'HorizontalAlignment','right','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','w');               
%                 
% 
% subplot(2,4,3)
% %imagesc(T_,fliplr(X_'),fliplr(Xt'))
% imagesc(PX_,X_,Xpx)
% if umscale
%     ylabel('x [um]');
% else
%     ylabel('x [m]');
% end
% xlabel('x''');
% set(gca,'YDir','normal')
% 
% subplot(2,4,7)
% %imagesc(T_,fliplr(X_'),fliplr(Xt'))
% imagesc(PY_,Y_,Ypy)
% if umscale
%     ylabel('y [um]');
% else
%     ylabel('y [m]');
% end
% xlabel('y''');
% set(gca,'YDir','normal')
% 
% subplot(2,4,5)
% plot(T_,It,'linewidth',2,'color','b');
% grid on;
% if lengthscale
%     xlabel('s [m]');
% else
%     xlabel('t [fs]');
% end
% ylabel('I[A]');
% xlim([min(T_) max(T_)]);
% text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'top','FontSize',10,'units','normalized','color','k');

return