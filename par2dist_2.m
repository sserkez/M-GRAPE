clear all
%orig_distr = 'unddist10M.dat';
%orig_distr = 'dist_700.dat';

%  nm_p='C:\-D-\Work\LCLS\tmp\3\530_tdp\530_u1_tdp_7.out';
%   nm_p='C:\-D-\Work\LCLS\tmp\3\damage\1200_u1_1.out';
nm_p='C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_11.out';
%nm_p='C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u1_9.out';
%nm_p='C:\-D-\Work\LCLS\tmp\3\s-n study\1\1200_u1_3.out';

nm1=[nm_p,'.dpa'];
nm2=[nm_p,'.dat'];
% nm1='C:\-D-\Work\SASE3_chicane\run3\U1.1.out.dpa';
% nm2='C:\-D-\Work\SASE3_chicane\run3\U1.1.out.dat';

d=outread(nm_p,1);
d.outp.power=rmfield(d.outp.power,'v');

I=flipud(d.outp.current);
Ti=d.outp.Sscale/3e8;

 npart=d.inp.npart;
 xlamds=d.inp.xlamds;
 zsep=d.inp.zsep;
% 


    % %read current
    % fileID = fopen(orig_distr,'r');
    % delimiter = ' ';
    % startRow = 5;
    % formatSpec = '%f%f%f%f%f%f%*s%*s%[^\n\r]';
    % dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    % fclose(fileID);
    % T = dataArray{:, 5};
    % I=hist(T,100);
    % Ti=linspace(0,max(T),100);
    % figure(130)
     plot(Ti,I);
    % 
    % clear T dataArray

%%

fd1=fopen(nm1,'r');
if fd1==-1
    error('.par file not found')
end
% FileInfo = dir([pwd,'\',nm1]);
FileInfo = dir(nm1);
N_records=FileInfo.bytes/(8);
N_slice=floor(N_records/npart/6); %kostyl'
x=fread(fd1,N_records,'double');
x=single(x);


%x=x(N_records+1:end);

%N_records=N_records/2;
%N_slice=N_records/npart/6;

xx=zeros(6,npart,round(N_slice),'single');

% for n=1:N_slice
%     for j=1:6
%         for i=1:NPART
%             xx(j,i,n)=x(NPART*6*(n-1)+NPART*(j-1)+i);
%         end
%     end
% end

for n=1:N_slice
    for j=1:6
        for i=1:npart
            xx(j,i,n)=(x(npart*6*(n-1)+npart*(j-1)+i));
        end
    end
end

% figure
% for i=1:N_slice
% hold on
% gamma=xx(1,:,i);
% phase=xx(2,:,i)+i*25*2*pi;%+(rand(1,1024)-0.5)*2*pi*25;
% plot(phase/2/pi*1.24/1000,gamma,'.');
% hold off
% end
%%
% E=reshape(xx(1,:,:),[N_slice*npart 1]);
% Ph=reshape(xx(2,:,:),[N_slice npart]);
% for i=1:npart
%     Ph(:,i)=Ph(:,i)+(i-1)*zsep*2*pi+(rand(N_slice,1)-0.5)*2*pi*zsep;
% end
% T=Ph/2/pi*xlamds/299792458;
% T=reshape(T,[N_slice*npart 1]);
% 
% clear Ph;
% % Ph=reshape(xx(2,:,:),[N_slice*npart 1]);
% X=reshape(xx(3,:,:),[N_slice*npart 1]);
% Y=reshape(xx(4,:,:),[N_slice*npart 1]);
% XP=reshape(xx(5,:,:),[N_slice*npart 1]);
% YP=reshape(xx(6,:,:),[N_slice*npart 1]);
% 
% index=E<1000;
% T(index)=[];
% X(index)=[];
% Y(index)=[];
% XP(index)=[];
% YP(index)=[];
% E(index)=[];
% 
% npart=size(T,1);
% load('current.mat');
% curpeak=interp1(t,curpeak,linspace(0,max(max(T)),N_slice));
% curpeak=round(curpeak/max(curpeak)*npart);
% curpeak(isnan(curpeak))=0;
% curpeak(1)=1;
% 
% index=logical(zeros(N_slice*npart,1));
% for i=1:N_slice
%     index((i-1)*npart+1:(i-1)*npart+curpeak(i)+1)=1;
% end
% 
% T(index)=[];
% X(index)=[];
% Y(index)=[];
% XP(index)=[];
% YP(index)=[];
% E(index)=[];
% 
% T=(-T+max(T))/2;
% XP=XP./E;
% YP=YP./E;

E=reshape(xx(1,:,:),[npart N_slice]);

%tweak
Ph=reshape(xx(2,:,:),[npart N_slice]);

for i=1:N_slice
    Ph(:,i)=Ph(:,i)+(i-1)*zsep*2*pi+(rand(npart,1)-0.5)*2*pi*zsep;
end
Ph=-Ph+max(max(Ph));
%endtweak

T=Ph/2/pi*xlamds/299792458;
%T=reshape(T,[N_slice*npart 1]);

clear Ph;

% Ph=reshape(xx(2,:,:),[N_slice*npart 1]);
X=reshape(xx(3,:,:),[npart N_slice]);
Y=reshape(xx(4,:,:),[npart N_slice]);
XP=reshape(xx(5,:,:),[npart N_slice]);
YP=reshape(xx(6,:,:),[npart N_slice]);
clear xx

%load('current.mat');
I=interp1(Ti,I,mean(T,1)-min(T(E>100)));
I(isnan(I))=0;
I=round(I/max(I)*npart);
%curpeak(1)=1;
%curpeak=fliplr(curpeak);
figure(130)
plot(mean(T,1),I);

% index=(X==0&Y==0);
% T(index)=[];
% X(index)=[];
% Y(index)=[];
% XP(index)=[];
% YP(index)=[];
% E(index)=[];


% index1=logical(zeros(N_slice*npart,1));
% for i=1:N_slice
%     index1((i-1)*npart+1:(i-1)*npart+curpeak(i)+1)=1;
% end
index2=T<0|(X==0&Y==0);
index2=T<0|E<100;
index2=E<100;


index3=index2;
k=0;
for i=1:N_slice
    for j=1:npart
        if(index2(j,i)==0)&k<npart-I(i)
         index3(j,i)=1;
         k=k+1;
%          if k==sum(index2(i,:)==1)|k==npart-curpeak(i);
%              break
%          end
        end
    end
    k=0;
end
index=index3;


%break
T(index)=[];
X(index)=[];
Y(index)=[];
XP(index)=[];
YP(index)=[];
E(index)=[];

XP=XP./E;
YP=YP./E;
%T=T-min(min(T));

%%
skipstep=ceil(numel(T)/1000000);
%skipstep=6;
index=ones(size(T,2),1);
index(1:skipstep:end)=round(0);
index=index==1;
T(index)=[]; T=rot90(T);
X(index)=[]; X=rot90(X);
Y(index)=[]; Y=rot90(Y);
XP(index)=[]; XP=rot90(XP);
YP(index)=[]; YP=rot90(YP);
E(index)=[]; E=rot90(E);
%figure(135)

%%
nT=300;
nE=200;
Ed = hist3([T*3e8 E],[nT nE]); 
%plot(T*3e8,E,'.');
% figure(10)
% plot(T*299792458,E,'.');
 T_=linspace(min(T.*3e8),max(T.*3e8),nT);
 E_=linspace(min(E),max(E),nE);
figure(1001)
%pcolor(T_,E_,Ed')
imagesc(T_.*1e6,E_,fliplr(Ed)')
xlabel('s [\mum]');
ylabel('\gamma');
drawnow

%break
%% Write distribution

%nrg=6584; %500
%nrg=7785; %700
%nrg=9293; %1000

dlmwrite(nm2, '? VERSION = 1.0 ', 'delimiter', '');
dlmwrite(nm2, ['? SIZE = ',num2str(size(T,1))], 'delimiter', '','-append');
dlmwrite(nm2, ['? CHARGE = ' num2str(d.inp.charge)], 'delimiter', '','-append');
dlmwrite(nm2, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');

%dlmwrite(nm2,[X_val-4e-5 X_mom-mean(X_mom) Y_val Y_mom T_val+6.1e-14 Gamma],'delimiter',' ','precision',8,'-append');
dlmwrite(nm2,[X XP Y YP T E],'delimiter',' ','precision',8,'-append');

fclose all