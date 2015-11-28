function [emitx] =  beam_matching(filein, fileout, mbetax, malphax, mbetay, malphay)
% this is script applies a linear map to change the twiss parameters of an
% e-beam for matching purposes

%filein  = name of the file to be imported, should have no header
% columns: x, xp, y, yp, t, gamma

%fileout = name of the matched e-beam file

%datain=importdata(filein,' ',35); %if in ELEGANT format
%datain=importdata(filein,' ',5); %if in GENESIS format
datain=importdata(filein,' ',39); %%%% change this one depending on who gives the particle file (lanfa was 20)
text=datain.textdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=datain.data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size(data)

tmin=min(data(:,5))
tmax=max(data(:,5))
tavg = mean(data(:,5))
data(:,5)=data(:,5)-tavg;
deltat=(tmax-tmin)%./299792458

deltag=-8.377132360193900e+01;%-9.123109435339065e+02 %any gamma shift?

xavg=mean(data(:,1))
xpavg=mean(data(:,2))
yavg=mean(data(:,3))
ypavg=mean(data(:,4))


x=data(:,1)-xavg;
xp=data(:,2)-xpavg;
y=data(:,3)-yavg;
yp=data(:,4)-ypavg;
gamma=sqrt(data(:,6).^2+1)+deltag;  %%%%%this line depends on whether or not p or gamma is input column 6

gavg=mean(gamma)

emitx=sqrt(mean(x.^2).*mean(xp.^2)-mean(x.*xp).^2).*gavg
emity=sqrt(mean(y.^2).*mean(yp.^2)-mean(y.*yp).^2).*gavg
betax=mean(x.*x).*gavg./emitx
betay=mean(y.*y).*gavg./emity
alphax=-mean(x.*xp).*gavg./emitx
alphay=-mean(y.*yp).*gavg./emity

figure
plot(data(:,5),xp,'.')
title('(t,xp) before match')


figure
subplot(2,2,1)
plot(data(:,5),data(:,6),'.');
title('Longitudinal Phase Space (t,\gamma)');
xlabel('t');
ylabel('\gamma');
subplot(2,2,2)
plot(x,y,'.');
title('Transverse Distribution');
xlabel('x');
ylabel('y');
subplot(2,2,3)
plot(x,xp,'.');
title('X Phase Space');
xlabel('x');
ylabel('xp');
subplot(2,2,4)
plot(y,yp,'.');
title('Y Phase Space');
xlabel('y');
ylabel('yp');

%remove old correlation
xp=xp+alphax.*x./betax;
yp=yp+alphay.*y./betay;

%scale the beam
x=x.*sqrt(mbetax./betax);
y=y.*sqrt(mbetay./betay);
xp=xp.*sqrt(betax./mbetax);
yp=yp.*sqrt(betay./mbetay);

%add new correlation
xp=xp-malphax.*x./mbetax;
yp=yp-malphay.*y./mbetay;

emitxnew=sqrt(mean(x.^2).*mean(xp.^2)-mean(x.*xp).^2).*gavg
emitynew=sqrt(mean(y.^2).*mean(yp.^2)-mean(y.*yp).^2).*gavg
betaxnew=mean(x.*x).*gavg./emitxnew
betaynew=mean(y.*y).*gavg./emitynew
alphaxnew=-mean(x.*xp).*gavg./emitxnew
alphaynew=-mean(y.*yp).*gavg./emitynew

figure
plot(data(:,5),xp,'.')
title('(t,xp) after match')

figure
subplot(2,2,1)
plot(data(:,5),gamma,'.');
title('Matched Longitudinal Phase Space (t,\gamma)');
xlabel('t');
ylabel('\gamma');
subplot(2,2,2)
plot(x,y,'.');
title('Matched Transverse Distribution');
xlabel('x');
ylabel('y');
subplot(2,2,3)
plot(x,xp,'.');
title('Matched X Phase Space');
xlabel('x');
ylabel('xp');
subplot(2,2,4)
plot(y,yp,'.');
title('Matched Phase Space');
xlabel('y');
ylabel('yp');

figure
plot(data(:,5),gamma,'.');

%figure
%hist2d([data(:,5)';gamma'],500,500);
%title('Matched Longitudinal Phase Space (t,\gamma)');
%xlabel('t');
%ylabel('\gamma');

if(1==2)

figure
subplot(2,2,1)
hist2d([data(:,5)';gamma'],500,500);
title('Matched Longitudinal Phase Space (t,\gamma)');
xlabel('t');
ylabel('\gamma');
subplot(2,2,2)
hist2d([x';y'],500,500);
title('Matched Transverse Distribution');
xlabel('x');
ylabel('y');
subplot(2,2,3)
hist2d([x';xp'],500,500);
title('Matched X Phase Space');
xlabel('x');
ylabel('xp');
subplot(2,2,4)
hist2d([y';yp'],500,500);
title('Matched Phase Space');
xlabel('y');
ylabel('yp');

end

newdata=[x xp y yp data(:,5) gamma];
%charge=strcat('? charge = ',text(27))
%num=strcat('? size = ',text(28))
fid = fopen(fileout,'w');
%fprintf(fid, '%12s\n', '# S2E LCLS 145.6pC 1.5kA 3.5 GeV SF');
fprintf(fid, '%12s\n', '? version = 1.0');
fprintf(fid, '%12s\n', '? charge = write!e-12');
fprintf(fid, '%12s\n', '? size = 973051');
fprintf(fid, '%12s\n', '? COLUMNS X XPRIME Y YPRIME T GAMMA');

fprintf(fid, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',newdata');
fclose(fid);
size(newdata)

Qperpart = 145.6e-12/973051;
[n,xout]=hist(data(:,5),600);
%figure
%bar(xout,n*Qperpart/(xout(2)-xout(1)));
figure
plot(xout,n*Qperpart/(xout(2)-xout(1)))

temp = importdata(fileout,'\t',5);
figure
plot(temp.data(:,5),temp.data(:,6),'.')

NSLICE = round((max(temp.data(:,5))-min(temp.data(:,5))).*299792458/(2.339320754716981e-09*5))

[n,xout]=hist(temp.data(:,5),600);
%figure
%bar(xout,n*Qperpart/(xout(2)-xout(1)));
figure
plot(xout,n*Qperpart/(xout(2)-xout(1)))







