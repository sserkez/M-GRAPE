%plots amplitude and phase
%H is link object and needs to be stored as a value outside  object so
%plots are rescaled properly

function [H]=fieldplot_ff(Nf,X,xlen,xlamds,name,showpictures)
if showpictures
%         function FWHM = findFWHM(x, fx)
%         %findFWHM: Finds the full width half max (FWHM) of a function
%         %
%         %	FWHM = findFWHM(x, fx);


        %u1=gpuArray(u1);

        %fx= gpuArray.linspace(-1/(2*dx)+1/L/2,1/(2*dx)-1/L/2,M);
  %      [FX,FY]=meshgrid(fx,fx);
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
    xscale=linspace(-1/(2*dx)+1/xlen/2,1/(2*dx)-1/xlen/2,Mx).*xlamds;
    yscale=linspace(-1/(2*dy)+1/ylen/2,1/(2*dy)-1/ylen/2,My).*xlamds;
%     xscale=((Mx-1)/2+1-(1:Mx))*dx;
%     yscale=((My-1)/2+1-(1:My))*dy;
    P0=sum(sum(abs(X).^2));
    Ix=sum(sum(abs(X).^2,3),1)'./sum(sum(sum(abs(X).^2)));
    Iy=sum(sum(abs(X).^2,3),2)./sum(sum(sum(abs(X).^2)));
    mean_x=sum(xscale.*Ix');
    mean_y=sum(xscale.*Iy');
    std_x=sqrt(sum((xscale-mean_x).^2.*Ix'));
    std_y=sqrt(sum((xscale-mean_y).^2.*Iy'));
    
    adc=0; %add column
    if N~=1
        Xo=X;
        %size(X)
        %size(repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My))
        %aX=angle(X)-repmat(angle(X((Mx-1)/2,(My-1)/2,:)),Mx,My);
        X=sqrt(mean(abs(X).^2,3)).*exp(1i.*angle(mean(X,3)));
        adc=1;
    end
    
        fig=figure(Nf);
                clf;
                set(fig,'name',[name,' (mesh=',num2str(round(Mx)),')'],'numbertitle','off');
                h1=subplot(2,2+adc,1);
                
                imagesc(xscale,yscale,abs(X).^2);
                axis equal tight;
%                 text(0,0,sprintf(' Power = %.2e ', P0),...
%                     'HorizontalAlignment','left','VerticalAlignment',...
%                     'bottom','FontSize',10,'units','normalized','color','w');
                title(h1,'Intensity');
                xlabel('[rad]');
                ylabel('[rad]');
                %colorbar('SouthOutside');
                
                h2=subplot(2,2+adc,2);
                %plot(yscale,sum(abs(X),2)');
                %view([90 -90])
                plot(Iy,yscale, 'linewidth',2);
                set(h2,'Ydir','reverse');
                axis square
                text(1,1,sprintf('FWHM = %.2e \n RMS= %.2e', findFWHM(yscale,Iy),std_y),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h2,'Y projection');   
                ylabel('[rad]');
                xlabel('I [a.u.]');
                
                h3=subplot(2,2+adc,3+adc);
                plot(xscale,Ix, 'linewidth',2);
                axis square
                %text(1,1,sprintf('FWHM = %.2e \n  sigma=%f', findFWHM(xscale,sum(abs(X),1)'), 3),...
                text(1,1,sprintf('FWHM = %.2e \n RMS= %.2e', findFWHM(xscale,Ix),std_x),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized');
                title(h3,'X projection');   
                xlabel('[rad]');
                ylabel('I [a.u.]');
              
                h4=subplot(2,2+adc,4+adc);

                imagesc(xscale,yscale,angle(X));
                axis equal tight;
                xlabel('[rad]');
                ylabel('[rad]');
                %colormap('hsv(256)');
                %colorbar('SouthOutside');
                title(h4,'Phase');                
                linkaxes([h1,h4], 'xy');
                %linkaxes([h1,h2], 'x');
                hlinkx=linkprop([h1 h2], 'YLim');
                hlinky=linkprop([h1 h3], 'XLim');
                %get(h1);
                
                if N~=1
                    h5=subplot(2,2+adc,3);
                    imagesc(1:N,xscale,reshape(mean(abs(Xo).^2,2),Mx,N));
                    xlabel('[a.u.]');
                    ylabel('y length [rad]');
                    axis square
                    drawnow;
                    %title(['z=',num2str(prop_leng),'m']);

                    h6=subplot(2,2+adc,6);
                    imagesc(1:N,yscale,reshape(mean(abs(Xo).^2,1),My,N));
                    xlabel('[a.u.]');
                    ylabel('x length [rad]');
                    axis square
                    drawnow;
                    
                    hlinkz=linkprop([h5 h6], 'XLim');
                    H=[hlinkx,hlinky,hlinkz];
                else 
                    H=[hlinkx,hlinky];
                end

else
H=[];
end
                
    
