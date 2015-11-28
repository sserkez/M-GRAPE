function [k,filt,phs]=KKphase(mult,xlamds,resolshift,I_function)
% clear all
% close all



%datin = dlmread('profile.dat');
%TAKE AMPLITUDE

k0=2*pi/xlamds;
%mult = 5;

%lam = xlamds + datin(:,1)*xlamds*1.0E-5;
lam = xlamds + resolshift'*xlamds;
k   = 2*pi./lam(:);
dk = k(2)-k(1);
Ntot = length(k)*mult;
k = k0 + dk*linspace(-Ntot/2,Ntot/2,Ntot);

cont1 = zeros(1,Ntot/mult*(mult-1)/2)+I_function(1)';%*1e-10;
cont2 = zeros(1,Ntot/mult*(mult-1)/2)+I_function(end)';%*1e-10;
filt = cat(2,cont1,I_function,cont2);

ord=1;
nn=find(filt<100);
filt(1:nn(1))=filt(nn(1));
filt(nn(end):end)=filt(nn(end));
filt0=-log(filt)/2;



tic
phs=-ord*kkimbook2(k,filt0,0);

k(1:length(k))=k(length(k):-1:1);
phs(1:length(k))=phs(length(k):-1:1);
filt(1:length(k))=filt(length(k):-1:1);

% dataP = [k' phs' ];
% dataM = [k' filt'];
% save('Tpha_in.dat', 'dataP', '-ascii','-double')
% save('Tmod_in.dat', 'dataM', '-ascii','-double')
toc

figure(1001)
plot(2*pi./k,filt,'-r');
hold on
plot(2*pi./k,phs,'^-k');
title('1 C crystal 400')
grid on
drawnow
return


