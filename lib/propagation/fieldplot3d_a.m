%plots amplitude and phase
%H is link object and needs to be stored as a value outside  object so
%plots are rescaled properly

function [H]=fieldplot3d_a(Nf,X,Lx,Ly,tdomain,zscale,name,showpictures)
%X=X';
if showpictures
    plotphase=0;
    
    if size(zscale,2)~=size(X,3)
        error('zscale vector length not consistent with radation mesh size');
    end
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
    %Ly=Lx;
    [My,Mx,N]=size(X);
    dx=Lx/Mx;
    dy=Ly/My;
    xscale=((Mx-1)/2+1-(1:Mx))*dx;
    yscale=((My-1)/2+1-(1:My))*dy;

if N==1
    P0=sum(sum(sum(abs(X).^2)));
else
    if tdomain
        i1=round(N/2);
        i2=i1+1;
        zsep=zscale(i2)-zscale(i1);
        tsep=zsep/3e8;
        P0=sum(sum(sum(abs(X).^2)))*tsep;
        %P0=sum(sum(sum(abs(X).^2)));%remove!
    else
        i1=round(N/2);
        i2=i1+1;
        dk=zscale(i1)-zscale(i2);
        xlamds=1.77e-9;
        zsep=2*pi/N/xlamds/dk;
        tsep=zsep/3e8;
        P0=sum(sum(sum(abs(X).^2)))*abs((zscale(i1)*zscale(i2))/(zscale(i1)-zscale(i2)))/3e8;
        P0=sum(sum(sum(abs(X).^2)))*tsep;
        %P0=0;
        %P0=sum(sum(sum(abs(X).^2)));%remove!
    end
    
end
    %P0=sum(sum(sum(abs(X).^2)));%remove later
    
    adc=0; %add column
    if N~=1
        Xo=X;
        Io2=abs(Xo).^2; %3D intensity
        %size(X)
        %size(repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My))
        %aX=angle(X)-repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My);
        X=sqrt(sum(abs(X).^2,3)).*exp(1i.*angle(mean(X,3))); %2D field
        adc=1;
    end
    
        fig=figure(Nf);
        if N~=1
            set(gcf,'PaperPositionMode','auto');
             set(fig, 'Position', [100, 100, 1000, 620]);
        else
            set(gcf,'PaperPositionMode','auto');
             set(fig, 'Position', [100, 100, 620, 620]);
        end
                clf;
                set(fig,'name',[name,' (mesh=',num2str(round(Mx)),'x',num2str(round(My)),'x',num2str(round(N)),')'],'numbertitle','off');
                h1=subplot(2,2+adc,1);
                
                imagesc(xscale,yscale,abs(X).^2);
                drawnow
                axis tight square;
                if P0~=0
                    text(0,0,sprintf(' E = %.2e J', P0),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'bottom','FontSize',10,'units','normalized','color','w');
                end
                title(h1,'Intensity');
                xlabel('x [m]');
                ylabel('y [m]');
                drawnow;
                %colorbar('SouthOutside');
                
                h2=subplot(2,2+adc,2);
                %plot(yscale,sum(abs(X),2)');
                %view([90 -90])
                
                plot(sum(abs(X).^2,2)',yscale, 'linewidth',2);
                set(h2,'Ydir','reverse');
                axis square
                try
                    FWHM_y=findFWHM(yscale,sum(abs(X).^2,2));
                catch
                    FWHM_y=[];
                end
                text(1,1,sprintf('fwhm = %.2em', FWHM_y),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h2,'Y projection');   
                ylabel('y [m]');
                xlabel('I [arb. units]');
                
                h3=subplot(2,2+adc,3+adc);
                plot(xscale,sum(abs(X).^2,1)', 'linewidth',2);
                axis square
                try
                    FWHM_x=findFWHM(xscale,sum(abs(X).^2,1)');
                catch
                    FWHM_x=[];
                end
                text(1,1,sprintf('fwhm = %.2em ', FWHM_x),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h3,'X projection');   
                xlabel('x [m]');
                ylabel('I [arb. units]');
              
                hlinkx=linkprop([h1 h2], 'YLim');
                hlinky=linkprop([h1 h3], 'XLim');               

                if plotphase || N==1
                    h4=subplot(2,2+adc,4+adc);
                    imagesc(xscale,yscale,angle(X));
                    axis tight square;
                    %colormap('hsv(256)');
                    %colorbar('SouthOutside');
                    title(h4,'Phase');                
                    linkaxes([h1,h4], 'xy');
                    %linkaxes([h1,h2], 'x');
                    %get(h1);
                end
                
                if N~=1
                    
                    if ~plotphase
                        h4=subplot(2,2+adc,4+adc);
                        lzprojection=reshape(sum(sum(Io2,1),2),1,[]); %lambda/z projection
                        plot(zscale,lzprojection, 'linewidth',2);
                        axis square
                        xlim([min(zscale) max(zscale)]);
                        if tdomain==1
                            title('Z projection')
                            ylabel('I [arb. units]');
                            xlabel('s [m]');
                        else
                            title('Spectrum')
                            ylabel('I [a.u.]');
                            xlabel('\lambda [m]');
                        end
                    end
%                     size(Io2)
%                     size(sum(Io2,2))
%                     size(permute(sum(Io2,2),[1 3 2]))
                    h5=subplot(2,2+adc,3);
                    %imagesc(zscale,xscale,reshape(sum(Io2,2),Mx,N));
                    imagesc(zscale,yscale,permute(sum(Io2,2),[1 3 2]));

                        
                    
                    if tdomain==1
                        xlabel('s [m]');
                    else
                        xlabel('\lambda [m]');
                    end
                    title('side view');
                    ylabel('y [m]');
                    axis square
                    drawnow;
                    %title(['z=',num2str(prop_leng),'m']);

                    h6=subplot(2,2+adc,6);
                    imagesc(zscale,xscale,permute(sum(Io2,1),[2 3 1])); % <- problem

                    if tdomain==1
                        xlabel('s [m]');
                    else
                        xlabel('\lambda [m]');
                    end
                    title('top view');
                    ylabel('x [m]');
                    axis square
                    drawnow;
                    
                    hlinkz=linkprop([h5 h6], 'XLim');
                    
%                     h4=subplot(2,2+adc,4+adc);
%                     %plot(zscale,reshape(mean(mean(Io2,1),2),1,[]));
%                     plot(zscale,mean(reshape(mean(Io2,1),My,N),1));
                    if ~plotphase 
                        hlinkz=linkprop([h4 h5 h6], 'XLim');
                    end
                    
                    H=[hlinkx,hlinky,hlinkz];
                else 
                    H=[hlinkx,hlinky];
                end

else
H=[];
end
                
    
