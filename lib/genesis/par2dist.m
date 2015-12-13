function par2dist(nm_p)

nm1=[nm_p,'.dpa'];
nm2=[nm_p,'.dat'];

d=outread(nm_p,1,0,2);
d.outp.power=rmfield(d.outp.power,'v');

I=flipud(d.outp.current);
Ti=d.outp.Sscale/3e8;

 npart=d.inp.npart;
 xlamds=d.inp.xlamds;
 zsep=d.inp.zsep;

%%
disp(' -Processing .dpa file:');
disp(['  - opening ',nm1]);
fd1=fopen(nm1,'r');
% FileInfo = dir([pwd,'\',nm1]);
FileInfo = dir(nm1);
N_records=FileInfo.bytes/(8);
N_slice=floor(N_records/npart/6);
x=fread(fd1,N_records,'double');
x=single(x);
%disp('  +done');

disp('  - organizing .dpa contents');

xx=zeros(6,npart,round(N_slice),'single');

for n=1:N_slice
    for j=1:6
        for i=1:npart
            xx(j,i,n)=(x(npart*6*(n-1)+npart*(j-1)+i));
        end
    end
end
E=reshape(xx(1,:,:),[npart N_slice]);
Ph=reshape(xx(2,:,:),[npart N_slice]);
X=reshape(xx(3,:,:),[npart N_slice]);
Y=reshape(xx(4,:,:),[npart N_slice]);
XP=reshape(xx(5,:,:),[npart N_slice]);
YP=reshape(xx(6,:,:),[npart N_slice]);
clear xx
%disp('  +done');
%%
disp('  - smearing out microbunching');
for i=1:N_slice
    Ph(:,i)=Ph(:,i)+(i-1)*zsep*2*pi+(rand(npart,1)-0.5)*2*pi*zsep;
end
Ph=-Ph+max(max(Ph));
%endtweak
T=Ph/2/pi*xlamds/299792458;
%T=reshape(T,[N_slice*npart 1]);
clear Ph;
%disp('  +done');

% Ph=reshape(xx(2,:,:),[N_slice*npart 1]);

disp('  - calculating current');
%load('current.mat');
I=interp1(Ti,I,mean(T,1)-min(T(E>100)));
I(isnan(I))=0;
I=round(I/max(I)*npart);
% figure(130)
% plot(mean(T,1),I);

%index2=T<0|(X==0&Y==0);
index2=T<0|E<100;
disp('  - filtering particles');
index3=index2;
k=0;
for i=1:N_slice
    for j=1:npart
        if(index2(j,i)==0)&k<npart-I(i)
         index3(j,i)=1;
         k=k+1;
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
skipstep=round(numel(T)/200000);
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
% figure(135)
% plot(T*3e8,E,'.');
%plot(T*299792458,E,'.');

%% Write distribution
disp('  - writing .dat file');
dlmwrite(nm2, '? VERSION = 1.0 ', 'delimiter', '');
dlmwrite(nm2, ['? SIZE = ',num2str(size(T,1))], 'delimiter', '','-append');
dlmwrite(nm2, ['? CHARGE = ' num2str(d.inp.charge)], 'delimiter', '','-append');
dlmwrite(nm2, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');
dlmwrite(nm2,[X XP Y YP T E],'delimiter',' ','precision',8,'-append');
fclose all;
disp(' +done');
return