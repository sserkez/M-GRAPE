%filename1 = 'E:\WORK\New_matlab\1000eV_distfile.dat';
filename1 = 'C:\-D-\WORK\LCLS\tmp\par180pC_s.dat';
filename1 = 'C:\-D-\WORK\LCLS\tmp\3\ebeams\145pC1p5kA3p5GeVSF.dist';
%filename1 = 'unddist10M.dat';
%filename1 = 'dist_700.dat';

%nrg=6584; %500
%nrg=7785; %700
nrg=9293; %1000
ts=1.0; %transverse scale
%ts=1.0; %transverse scale

%filename1 = 'E:\WORK\New_matlab\150pCL1S25L2S362kA.dat';

delimiter = ' ';
startRow = 6;

formatSpec = '%f%f%f%f%f%f%*s%*s%[^\n\r]';

fileID = fopen(filename1,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);


X = dataArray{:, 1};
PX = dataArray{:, 2};
Y = dataArray{:, 3};
PY = dataArray{:, 4};
T = dataArray{:, 5};
Gamma = dataArray{:, 6};


clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%%
nT=300;
nX=200;
Xd = hist3([T*3e8 X],[nT nX]); 
%plot(T*3e8,E,'.');
% figure(10)
% plot(T*299792458,E,'.');
 T_=linspace(min(T),max(T),nT);
 X_=linspace(min(X),max(X),nX);
 T_=T_-min(T_);
figure(1001)
%pcolor(T_,E_,Ed')
%imagesc(T_.*1e6,fliplr(X_'),fliplr(Xd'))
imagesc(T_,fliplr(X_'),fliplr(Xd'))
xlabel('t [fs]');
ylabel('x [m]');
set(gca,'YDir','normal')

%%
nT=300;
nG=200;
Gd = hist3([T*3e8 Gamma],[nT nG]); 
%plot(T*3e8,E,'.');
% figure(10)
% plot(T*299792458,E,'.');
 T_=linspace(min(T.*3e8),max(T.*3e8),nT);
 G_=linspace(min(Gamma),max(Gamma),nG);
 T_=T_-min(T_);
figure(1002)
%pcolor(T_,E_,Ed')
%imagesc(T_.*1e6,fliplr(X_'),fliplr(Xd'))
imagesc(T_,fliplr(G_'),fliplr(Gd'))
xlabel('s [m]');
ylabel('\gamma');
set(gca,'YDir','normal')

%% Filter distribution file
%index=Gamma>9312|(Gamma>9302&T_val-min(T_val)>4.2e-14)|(Gamma<9282&T_val-min(T_val)<6.7e-14)|(Gamma<9270&T_val-min(T_val)<7.6e-14);
index=1:100:size(T,1);
%index=zeros(size(T,1),1);
%index(1:1:size(T_val,1))=1;
index=logical(index);

T(index)=[];
PX(index)=[];
X(index)=[];
Y(index)=[];
PY(index)=[];
Gamma(index)=[];
% % T(index)=[];
% % PX(index)=[];
% % X(index)=[];
% % Y(index)=[];
% % PY(index)=[];
% % Gamma(index)=[];
% % T=T-min(T);
% % Gamma=Gamma-mean(Gamma)+nrg;
X=X.*ts; PX=PX.*ts; Y=Y.*ts; PY=PY.*ts; 
%% Write distribution file
% nbins=500;
% t=linspace(min(VarName5),max(VarName5),nbins);
% Current=histc(VarName5,linspace(min(VarName5),max(VarName5),nbins));
% Current=flipud(Current/max(Current)*4400);
% z=t*299792458;
% z=(z-min(z))';
% plot(z,Current+1e-5);
% 
% filename2='beamfile.dat';
% dlmwrite(filename2, '? VERSION = 1.0 ', 'delimiter', '');
% dlmwrite(filename2, '? COLUMNS ZPOS CURPEAK  ', 'delimiter', '','-append');
% dlmwrite(filename2,[z Current+1e-5],'delimiter',' ','precision',10,'-append');

% figure
% plot(T-min(T),Gamma,'.');
% 
% %hist(T_val-min(T_val),100);
%% calculate current

I=hist(T,500);
Ti=linspace(0,max(T),500);
I=I./sum(I.*Ti(2)).*145.6e-12;
figure(6543)
plot(Ti,I);

%% Write

filename2='C:\-D-\WORK\LCLS\tmp\dist_500_1.3.dat';

Gamma=Gamma-mean(Gamma)+nrg;

dlmwrite(filename2, '? VERSION = 1.0 ', 'delimiter', '');
dlmwrite(filename2, ['? SIZE = ',num2str(size(Gamma,1))], 'delimiter', '','-append');
dlmwrite(filename2, '? CHARGE = 1.8E-10 ', 'delimiter', '','-append');
dlmwrite(filename2, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');

%dlmwrite(filename2,[X_val-4e-5 X_mom-mean(X_mom) Y_val Y_mom T_val+6.1e-14 Gamma],'delimiter',' ','precision',8,'-append');
%dlmwrite(filename2,[X PX Y PY T-min(T) Gamma],'delimiter',' ','precision',8,'-append');
dlmwrite(filename2,[X PX Y PY T Gamma],'delimiter',' ','precision',8,'-append');

clearvars -except I Ti
fclose all;