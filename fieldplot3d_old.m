%plots amplitude and phase
%H is link object and needs to be stored as a value outside  object so
%plots are rescaled properly

function [H]=fieldplot3d(Nf,X,xlen,tdomain,zscale,name,showpictures)
if showpictures

%         function FWHM = findFWHM(x, fx)
%         %findFWHM: Finds the full width half max (FWHM) of a function
%         %
%         %	FWHM = findFWHM(x, fx);
% 
% 
%         [m, n] = max(fx);		%	Find maximum value and index
%         %FWHM = interp1(fx(end:-1:n), x(end:-1:n), m/2, 'spline') - interp1(fx(1:n), x(1:n), m/2, 'spline');
%         ind = find(fx>m/2);	%	Find indicies where I>=max(I)/2
%         nl = min(ind);			%	Leftmost index
%         nr = max(ind);			%	Rightmost index
% 
%         %	Linear interpolate x positions
%         xl = (x(nl)-x(nl-1))*(m/2-fx(nl-1))/(fx(nl)-fx(nl-1)) + x(nl-1);
%         xr = (x(nr)-x(nr-1))*(m/2-fx(nr-1))/(fx(nr)-fx(nr-1)) + x(nr-1);
% 
%         %	Get FWHM
%         FWHM = abs(xr-xl);
%         end
    ylen=xlen;
    [Mx,My,N]=size(X);
    dx=xlen/Mx;
    dy=ylen/My;
    xscale=((Mx-1)/2+1-(1:Mx))*dx;
    yscale=((My-1)/2+1-(1:My))*dy;

    if tdomain
        P0=sum(sum(sum(abs(X).^2)))/N*abs(zscale(end)-zscale(1))/3e8;
    else
        i1=round(N/2);
        i2=i1+1;
        P0=sum(sum(sum(abs(X).^2)))/N*abs((zscale(i1)*zscale(i2))/(zscale(i1)-zscale(i2)))/3e8;
        P0=0;
    end
    %P0=sum(sum(sum(abs(X).^2)));%remove later
    
    adc=0; %add column
    if N~=1
        Xo=X;
        %size(X)
        %size(repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My))
        %aX=angle(X)-repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My);
        X=sqrt(sum(abs(X).^2,3)).*exp(1i.*angle(mean(X,3)));
        adc=1;
    end
    
        fig=figure(Nf);
        %set(fig, 'Position', [100, 100, 1049, 895]);
                clf;
                set(fig,'name',[name,' (mesh=',num2str(round(Mx)),')'],'numbertitle','off');
                h1=subplot(2,2+adc,1);
                
                imagesc(xscale,yscale,abs(X).^2);
                drawnow
                axis equal tight;
                if P0~=0
                    text(0,0,sprintf(' Power = %.2e J', P0),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'bottom','FontSize',10,'units','normalized','color','w');
                end
                title(h1,'Intensity');
                xlabel('x [m]');
                ylabel('y [m]');
                %colorbar('SouthOutside');
                
                h2=subplot(2,2+adc,2);
                %plot(yscale,sum(abs(X),2)');
                %view([90 -90])
                plot(sum(abs(X).^2,2)',yscale, 'linewidth',2);
                set(h2,'Ydir','reverse');
                axis square
                try
                    FWHM_y=findFWHM(yscale,sum(abs(X).^2,2)');
                catch
                    FWHM_y=[];
                end
                text(1,1,sprintf('FWHM = %.2e ', FWHM_y),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h2,'Y projection');   
                ylabel('[m]');
                xlabel('I [a.u.]');
                
                h3=subplot(2,2+adc,3+adc);
                plot(xscale,sum(abs(X).^2,1)', 'linewidth',2);
                axis square
                try
                    FWHM_x=findFWHM(yscale,sum(abs(X).^2,1)');
                catch
                    FWHM_x=[];
                end
                text(1,1,sprintf('FWHM = %.2e ', FWHM_x),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h3,'X projection');   
                xlabel('[m]');
                ylabel('I [a.u.]');
              
                hlinkx=linkprop([h1 h2], 'YLim');
                hlinky=linkprop([h1 h3], 'XLim');               

                    h4=subplot(2,2+adc,4+adc);
                    imagesc(xscale,yscale,angle(X));
                    axis equal tight;
                    %colormap('hsv(256)');
                    %colorbar('SouthOutside');
                    title(h4,'Phase');                
                    linkaxes([h1,h4], 'xy');
                    %linkaxes([h1,h2], 'x');
                    %get(h1);
                
                if N~=1
                    
                    
                    h5=subplot(2,2+adc,3);
                    imagesc(zscale,xscale,reshape(mean(abs(Xo).^2,2),Mx,N));
                    if tdomain==1
                        xlabel('z length [m]');
                    else
                        xlabel('wavelength [m]');
                    end
                    ylabel('y length [m]');
                    axis square
                    drawnow;
                    %title(['z=',num2str(prop_leng),'m']);

                    h6=subplot(2,2+adc,6);
                    imagesc(zscale,yscale,reshape(mean(abs(Xo).^2,1),My,N)); % <- problem

                    if tdomain==1
                        xlabel('z length [m]');
                    else
                        xlabel('wavelength [m]');
                    end
                    ylabel('x length [m]');
                    axis square
                    drawnow;
                    
                    hlinkz=linkprop([h5 h6], 'XLim');
                    
%                     h4=subplot(2,2+adc,4+adc);
%                     %plot(zscale,reshape(mean(mean(abs(Xo).^2,1),2),1,[]));
%                     plot(zscale,mean(reshape(mean(abs(Xo).^2,1),My,N),1));

                    H=[hlinkx,hlinky,hlinkz];
                else 
                    H=[hlinkx,hlinky];
                end

else
H=[];
end
                
    
