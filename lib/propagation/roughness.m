function [x2,XX,hh] = roughness(x1, leng, teta_i, teta_d, K, seed)
teta=teta_i+teta_d;

%k=0;    % 0 nm
k=0.3;  % 1 nm
%k=0.98; % 2 nm
%k=3.9;  % 4 nm

if k==0
    x2=x1;
    NN0=1000;
    LMR=1e5;
    XX=linspace(0,1e-6*LMR,NN0);
    hh=zeros(1,NN0); 
else
    
    [M,~]=size(x1);
    dx=leng/M;
    [xx,yy]=meshgrid(((M-1)/2+1-(1:M)).*dx);
    xscale=((M-1)/2+1-(1:M))*dx;

    rng(seed)

    LMR=1e5;    %mirror length (mikrons)
    NN0=1000;   %number of points in mirror
    XX=linspace(0,1e-6*LMR,NN0);    %space along mirror (m)
    hh=zeros(1,NN0);    
    PP=linspace(1e-5,0.01,NN0);     %PSD argument
    Sg0=exp(-1.3*log(PP))*1.4e-8;    %PSD function

    % Sg0=Sg0.*0.3; % 1 nm
    % Sg0=Sg0.*0.98;% 2 nm
    % Sg0=Sg0.*3.9; % 4 nm
     Sg0=Sg0.*k;

    Sg=2*sqrt(2*Sg0/LMR)*NN0;

    Nc=NN0;
    ramp=randn(1,Nc);
    rphi=rand(1,Nc);

    hh=imag(ifft(Sg.*ramp.*exp(2*pi*1i*rphi)));
    hh=hh/1e6;
    XX=XX-max(XX)/2;


    %disp(['heigth rms=',num2str(std(hh)*1e9)]);
    % plot(XX,hh,'k')
    %  title(['FFT RMS ' num2str(std(hh)) '   Mean ' num2str(mean(hh))]);
    % grid on
    % drawnow
    disp(['field density/profile density',num2str(dx/tan(teta_d)/(1e-6*LMR/NN0))]);
    hh1=interp1(XX,hh,xscale./tan(teta_d),'spline');

    hh1=hh1([linspace(1,1,numel(hh1))],:);
    x2=abs(x1).*exp(1i*(angle(x1)+K*hh1.'*sin(teta)));
end

end