function [x2,XX,hh] = roughness_1(x1, leng, Theta_i, Theta_d, xlamds, nm)
Theta=Theta_i+Theta_d;

load(nm);
Distance=(Distance-(max(Distance)-min(Distance))/2)*1e-3;
Height=Height*1e-9;

K=2*pi/xlamds;



% if k==0
%     x2=x1;
%     NN0=1000;
%     LMR=1e5;
%     XX=linspace(0,1e-6*LMR,NN0);
%     hh=zeros(1,NN0); 
% else
    
    [M,~,N]=size(x1);
    dx=leng/M;
    [xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);
    xscale=((M-1)/2+1-(1:M))*dx;

%    rng(seed)

% %     LMR=1e5;    %mirror length (mikrons)
% %     NN0=1000;   %number of points in mirror
% %     XX=linspace(0,1e-6*LMR,NN0);    %space along mirror (m)
% %     hh=zeros(1,NN0);    
% %     PP=linspace(1e-5,0.01,NN0);     %PSD argument
% %     Sg0=exp(-1.3*log(PP))*1.4e-8;    %PSD function
% % 
% %  
% %      Sg0=Sg0.*k;
% % 
% %     Sg=2*sqrt(2*Sg0/LMR)*NN0;
% % 
% %     Nc=NN0;
% %     ramp=randn(1,Nc);
% %     rphi=rand(1,Nc);
% % 
% %     hh=imag(ifft(Sg.*ramp.*exp(2*pi*1i*rphi)));
% %     hh=hh/1e6;
% %     XX=XX-max(XX)/2;

    XX=Distance;
    hh=Height;
    %disp(['heigth rms=',num2str(std(hh)*1e9)]);
    % plot(XX,hh,'k')
    %  title(['FFT RMS ' num2str(std(hh)) '   Mean ' num2str(mean(hh))]);
    % grid on
    % drawnow
%    disp(['field density/profile density',num2str(dx/tan(Theta_d)/(1e-6*LMR/NN0))]);

XX1=xscale./tan(Theta_d);

if min(XX1)<min(XX)
    XX=[min(XX1); XX];
    hh=[hh(1); hh];
end
if max(XX1)>max(XX)
    XX=[XX; max(XX1)];
    hh=[hh; hh(end)];
end

    hh0=interp1(XX,hh,xscale./tan(Theta_d),'pchip');

    hh1=hh0([linspace(1,1,numel(hh0))],:);
    
    x2=x1;
    for i=1:N
        x2(:,:,i)=abs(x1(:,:,i)).*exp(1i*(angle(x1(:,:,i))+K*hh1.*sin(Theta)));
    end
    
    I=sum(sum(abs(x1).^2,1),3);
    I=I/max(I);
    
    if N==1 && 1 || 1
        figure
        set(gcf,'name',nm);
        %title(nm);
        plot(xscale,I);
        hold all
        plot(xscale,hh0/max(hh0));
        hold off
    end
%     %./tan(Theta_d)
end