clear all
nm1='1000_u1_tdp.out.dpa';

N_part=2048;
xlamds=1.243e-9;
zsep=25+2; %костиль
fd1=fopen(nm1,'r');
FileInfo = dir([pwd,'\',nm1]);
N_records=FileInfo.bytes/(8);
N_slice=N_records/N_part/6;
x=fread(fd1,N_records,'double');
x=single(x);


%x=x(N_records+1:end);

%N_records=N_records/2;
%N_slice=N_records/N_part/6;

xx=zeros(6,N_part,round(N_slice));

% for n=1:N_slice
%     for j=1:6
%         for i=1:NPART
%             xx(j,i,n)=x(NPART*6*(n-1)+NPART*(j-1)+i);
%         end
%     end
% end

for n=1:N_slice
    for j=1:6
        for i=1:N_part
            xx(j,i,n)=x(N_part*6*(n-1)+N_part*(j-1)+i);
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
E=reshape(xx(1,:,:),[N_slice*N_part 1]);
Ph=reshape(xx(2,:,:),[N_slice N_part]);
for i=1:N_part
    Ph(:,i)=Ph(:,i)+(i-1)*zsep*2*pi+(rand(N_slice,1)-0.5)*2*pi*zsep;
end
Ph=reshape(Ph,[N_slice*N_part 1]);
T=Ph/2/pi*xlamds/299792458;
clear Ph;
% Ph=reshape(xx(2,:,:),[N_slice*N_part 1]);
X=reshape(xx(3,:,:),[N_slice*N_part 1]);
Y=reshape(xx(4,:,:),[N_slice*N_part 1]);
XP=reshape(xx(5,:,:),[N_slice*N_part 1]);
YP=reshape(xx(6,:,:),[N_slice*N_part 1]);

index=E<1000;
T(index)=[];
X(index)=[];
Y(index)=[];
XP(index)=[];
YP(index)=[];
E(index)=[];

T=(-T+max(T))/2;
XP=XP./E;
YP=YP./E;
%%
figure
plot(T,E,'.');

%plot(T*299792458,E,'.');
load('current.mat');
curpeak=interp1(t,curpeak,linspace(0,max(T),N_slice));
break
%% Write distribution

index=ones(size(T,1),1);
index(1:4:end)=round(0);
index=index==1;
T(index)=[];
X(index)=[];
Y(index)=[];
XP(index)=[];
YP(index)=[];
E(index)=[];


filename2='1000eV_u2_distfile.dat';
%nrg=6584; %500
%nrg=7785; %700
%nrg=9293; %1000

dlmwrite(filename2, '? VERSION = 1.0 ', 'delimiter', '');
dlmwrite(filename2, ['? SIZE = ',num2str(size(T,1))], 'delimiter', '','-append');
dlmwrite(filename2, '? CHARGE = 1.50E-10 ', 'delimiter', '','-append');
dlmwrite(filename2, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');

%dlmwrite(filename2,[X_val-4e-5 X_mom-mean(X_mom) Y_val Y_mom T_val+6.1e-14 Gamma],'delimiter',' ','precision',8,'-append');
dlmwrite(filename2,[X XP Y YP T E],'delimiter',' ','precision',8,'-append');
