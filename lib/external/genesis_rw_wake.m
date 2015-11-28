function genesis_rw_wake(currentFile,outputFile)
%	to calculate rw_wake for LCLS undulator chamber, Al, rectangle, 5mm gap
%	are hardcoded. 
%   a current profile [z(m) current (A)] will be read, with bunch head on the left (from originally elegant2current used for Genwake by Sven).
%   output will be saved in outputFile.
% genesis_rw_wake('undcur.dat','LCLSwake.dat')

sig  = 3.5e7;  % Al: 'Conductivity (ohm-1*m-1)'
tau  = 8e-15;      % Al: relaxation time
rf   =  1;              % rf=1: rectangle chamber: rf=0: round chamber    
r    =2.5;             % mm, chamber radius 

c  = 2.99792458E8;
Z0 = 120*pi;
offset = 1.9e4; % DC offset for the energy loss due to the wakefield.  Makes so no linear taper necessary.
%offset=0;

[zs Ipk] = textread(currentFile,'%f %f','delimiter',' ');

% figure
% plot(zs)

Q=integrate(zs/c,Ipk);
r  = r*1E-3;
s0 = (2*r^2/(Z0*sig))^(1/3);

f = Ipk/integrate(zs,Ipk);

s = zs - zs(1);
w = rw_wakefield(s,r,s0,tau,rf);

n = length(s);
E = zeros(n,n);
for j = 1:n
  for i = 1:n
    if i==j
      break
    else
      E(i,j) = w(j-i)*f(i);
    end
  end
end

dz = mean(diff(zs));
Ez = Q*sum(E)*dz; % eV/m/

Ez_mean = integrate(zs,f'.*Ez);
Ez_rms  = sqrt(integrate(zs,f'.*(Ez-Ez_mean).^2));
%Ez_rmsg = 100*rw_esprd(E0/1E9,Ne/1E10,L,r,sigz*1E6,sig);

zs=max(zs)-zs;
zs=flipud(zs);
Ez=flipud(Ez')-Ez_mean+offset;
Ipk=flipud(Ipk);
size(Ez)
size(Ipk)
%\kostyl
Ez=Ez-sum(Ez.*Ipk)./sum(Ipk);
sum(Ez.*Ipk)
%end\kostyl

%tstr = ['AC Resistive-Wall Wake ({\it\tau} = ' sprintf('%4.1f',tau*1E15) ' fs, {\it\sigma_c} = ' sprintf('%4.2f',sig/1E7) '\times10^7' ...
%         ' /\Omega/m, {\itr} = ' sprintf('%4.1f',r*1E3) ' mm'];
% Ipk=flipud(Ipk);
% Ez=flipud(Ez);


out=[zs Ipk Ez];
 header1 = '? VERSION=1.0';
 header2 = ['? SIZE=' num2str(length(zs))];
 header3 = '? COLUMNS ZPOS CURPEAK ELOSS';
 
fid = fopen(outputFile,'wt');
fprintf(fid,'%s\n',header1);
fprintf(fid,'%s\n',header2);
fprintf(fid,'%s\n',header3);
fprintf(fid,'%14.10e %14.10e %14.10e\n',out');
fclose(fid);


% figure
% plot(zs*1E6,Ez)
% xlabel('{\itz} (\mum)'); ylabel('Eloss (eV/m)')
% figure
% plot(zs*1E6,Ipk/1E3)
% xlabel('\mum'); ylabel('kA');
% 
% 
% figure
% plotyy(zs*1E6,Ipk,zs*1E6,Ez)


figure;
[AX,H1,H2] = plotyy(zs*1E6,Ipk*1e-3,zs*1E6,Ez*144*1e-6);
%xlim([mint maxt])
set(gca,'FontSize',20,'FontName','Times New Roman')
title('\DeltaE over length of HXR undulator','FontSize',20,'FontName','Times New Roman')
xlabel('s [\mum]','FontSize',20,'FontName','Times New Roman')
axes(AX(1));
%axis square
ylabel('I [kA]','FontSize',20,'FontName','Times New Roman') 
axes(AX(2));
%axis square
%xlim([mint maxt])
%ylim([-20 20])
ylabel('\DeltaE [MeV]','FontSize',20,'FontName','Times New Roman') 
set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',2)
hold;
h1=area(zs*1e6,ones(1,length(zs)).*2*5.55,-2*5.55,'EdgeColor','g','FaceColor','g');
alpha(0.2);
h2=area(zs*1e6,ones(1,length(zs)).*2*4.44,-2*4.44,'EdgeColor','b','FaceColor','b');
alpha(0.2);
legend([H1,H2,h1,h2],'I','\DeltaE','2\rho_{1D}E','2\rho_{3D}E')

if(1==2)
figure
[AX,H1,H2] = plotyy(zs*1E6,Ipk*1e-3,zs*1E6,Ez*90.753*1e-6);
%xlim([mint maxt])
set(gca,'FontSize',20,'FontName','Times New Roman')
title('\DeltaE over length of SXR undulator','FontSize',20,'FontName','Times New Roman')
xlabel('s [\mum]','FontSize',20,'FontName','Times New Roman')
axes(AX(1))
%axis square
ylabel('I [kA]','FontSize',20,'FontName','Times New Roman') 
axes(AX(2))
%axis square
%xlim([mint maxt])
%ylim([-20 20])
ylabel('\DeltaE [MeV]','FontSize',20,'FontName','Times New Roman') 
set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',2)
hold;
h1=area(zs*1e6,ones(1,length(zs)).*2*5.8,-2*5.8,'EdgeColor','g','FaceColor','g')
alpha(0.2)
h2=area(zs*1e6,ones(1,length(zs)).*2*4.0,-2*4.0,'EdgeColor','b','FaceColor','b')
alpha(0.2)
legend([H1,H2,h1,h2],'I','\DeltaE','2\rho_{1D}','2\rho_{3D}')
end



