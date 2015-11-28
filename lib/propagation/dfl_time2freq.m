function [X2,lamdscale]=dfl_time2freq(X,sscale,xlamds)

% xlamds=d(1, 1).inp.xlamds;
% sscale=d(1, 1).outp.Sscale;

N=size(X,3);

dk=2*pi/(max(sscale)-min(sscale));
sc=-(N-1)/2:1:(N-1)/2;
k=2*pi/xlamds;
%dk=2*pi/(N*xlamds*zsep);
K=k+dk.*sc;
lamdscale=2*pi./K;

X2=fftshift(fft(X,[],3),3)./sqrt(N);

spectrum_i_avg=reshape(mean(mean(abs(X2).^2,1),2),1,[]);

%  figure
%  plot(lamdscale,spectrum_i_avg);

return