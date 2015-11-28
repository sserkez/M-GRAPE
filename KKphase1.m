function [phs,haxes]=KKphase1(k,filt)
% clear all
% close all



%datin = dlmread('profile.dat');
%TAKE AMPLITUDE

% k0=2*pi/xlamds;
% %mult = 5;
% 
% %lam = xlamds + datin(:,1)*xlamds*1.0E-5;
% lam = xlamds + resolshift'*xlamds;
% k   = 2*pi./lam(:);
% dk = k(2)-k(1);
% Ntot = length(k)*mult;
% k = k0 + dk*linspace(-Ntot/2,Ntot/2,Ntot);
% 
% cont1 = zeros(1,Ntot/mult*(mult-1)/2)+I_function(1)';%*1e-10;
% cont2 = zeros(1,Ntot/mult*(mult-1)/2)+I_function(end)';%*1e-10;
% filt = cat(2,cont1,I_function,cont2);


ord=1;
% nn=find(filt<100);
% filt(1:nn(1))=filt(nn(1));
% filt(nn(end):end)=filt(nn(end));
filt0=-log(filt)/2;



%tic
phs=-ord*kkimbook2(k,filt0,0);

%k(1:length(k))=k(length(k):-1:1);
%phs(1:length(k))=phs(length(k):-1:1);
%filt(1:length(k))=filt(length(k):-1:1);

% dataP = [k' phs' ];
% dataM = [k' filt'];
% save('Tpha_in.dat', 'dataP', '-ascii','-double')
% save('Tmod_in.dat', 'dataM', '-ascii','-double')
%toc

figure(1001)
xlamds=4*pi/(k(1)+k(end));
%horzscale=2*pi./k*1e9; %lambda
horzscale=2*pi/xlamds./k-1; %dlpl

fwhm=1/findFWHM(horzscale,filt/max(filt));

[haxes,hline1,hline2] = plotyy(horzscale,filt,horzscale,phs);
ylabel(haxes(1),'Abs[T]') % label left y-axis
ylabel(haxes(2),'Ang[T],[rad]') % label right y-axis
% xlabel(haxes(2),'\lambda,nm') % label x-axis
xlabel(haxes(2),'\Delta\lambda/\lambda') % label x-axis
set(hline1,'LineWidth',2);
set(hline2,'LineStyle','--','LineWidth',2);

text(0.02,0.98,sprintf('FWHM^-^1 = %.0f ', fwhm),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',14,'units','normalized');

ylim(haxes(1), [0 0.06]);
ylim(haxes(2), [-3 3]);
set(haxes(2),'YTick',-3:1:3)
xlim(haxes(1), [min(horzscale) max(horzscale)]);
xlim(haxes(2), [min(horzscale) max(horzscale)]);
% plot(2*pi./k,filt,'-r','linewidth',2);
% hold on
% plot(2*pi./k,phs,'-k','linewidth',2);
% hold off
%title('Transmittance function')
%grid on
drawnow
return


