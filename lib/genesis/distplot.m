function distplot(fignum,dist,index_range)
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
    size=numel(T);
 
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
if lengthscale
    if umscale
        T=T.*3e8*1e6;
    else
        T=T.*3e8;
    end
else
    T=T.*1e15;
end

if umscale
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
%Xt = hist3([T X],[NT N]); 
Xt = hist3([X T],[N NT]); 
Yt = hist3([Y T],[N NT]);
% PXt = hist3([PX T],[N NT]); 
% PYt = hist3([PY T],[N NT]);
Gt = hist3([G T],[N NT]);
Yx = hist3([Y X],[N N]); 
% PXx = hist3([PX X],[N N]); 
% PYy = hist3([PY Y],[N N]);
Xpx = hist3([X PX],[N N]); 
Ypy = hist3([Y PY],[N N]);
PYpx = hist3([PY PX],[N N]);


fig=figure(fignum);
clf;
set(fig, 'Position', [100, 100, 1200, 600]);
subplot(2,4,1)
%[haxes,hline1,hline1] = plotyy([],[],[],[]);
%break
imagesc(T_,G_,Gt);

if lengthscale
    if umscale
        xlabel('s [\mum]');
    else
        xlabel('s [m]');
    end
else
    xlabel('t [fs]');
end
ylabel('\gamma');
set(gca,'YDir','normal')
text(0,1,sprintf('  head'),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');
text(1,1,sprintf('tail  '),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');
axis tight
% hold on
% axes;
% plot(T_,It);
% hold off

subplot(2,4,4)
imagesc(X_,Y_,Yx)
if umscale
    xlabel('x [um]');
    ylabel('y [um]');
else
    xlabel('x [m]');
    ylabel('y [m]');
end
%axis equal
set(gca,'YDir','normal')

subplot(2,4,8)
imagesc(PX_,PY_,PYpx)
    xlabel('x''');
    ylabel('y''');
%axis equal
set(gca,'YDir','normal')
                

subplot(2,4,2)
imagesc(T_,X_,Xt)
if lengthscale
    if umscale
        xlabel('s [\mum]');
    else
        xlabel('s [m]');
    end
else
    xlabel('t [fs]');
end
if umscale
    ylabel('x [um]');
else
    ylabel('x [m]');
end
set(gca,'YDir','normal')
text(0,1,sprintf('  head'),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');
text(1,1,sprintf('tail  '),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');
                
subplot(2,4,6)
imagesc(T_,Y_,Yt)

if lengthscale
    if umscale
        xlabel('s [\mum]');
    else
        xlabel('s [m]');
    end
else
    xlabel('t [fs]');
end

if umscale
    ylabel('y [um]');
else
    ylabel('y [m]');
end

set(gca,'YDir','normal')
text(0,1,sprintf('  head'),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');
text(1,1,sprintf('tail  '),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','w');               
                

subplot(2,4,3)
%imagesc(T_,fliplr(X_'),fliplr(Xt'))
imagesc(PX_,X_,Xpx)
if umscale
    ylabel('x [um]');
else
    ylabel('x [m]');
end
xlabel('x''');
set(gca,'YDir','normal')

subplot(2,4,7)
%imagesc(T_,fliplr(X_'),fliplr(Xt'))
imagesc(PY_,Y_,Ypy)
if umscale
    ylabel('y [um]');
else
    ylabel('y [m]');
end
xlabel('y''');
set(gca,'YDir','normal')

subplot(2,4,5)
plot(T_,It,'linewidth',2,'color','b');
grid on;
if lengthscale
    if umscale
        xlabel('s [\mum]');
    else
        xlabel('s [m]');
    end
else
    xlabel('t [fs]');
end
ylabel('I[A]');
xlim([min(T_) max(T_)]);
text(0,1,sprintf('  Q=%.1f pC', dist.charge*1e12),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','k');

return