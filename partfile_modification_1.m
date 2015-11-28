clear all
nm1='1000_u1_tdp.out.dpa';
nm2='1000_u1_tdp_new.out.dpa';
NPART=1024;
fd1=fopen(nm1,'r');
FileInfo = dir([pwd,'\',nm1]);
N_records=FileInfo.bytes/(8);
N_slice=N_records/NPART/6;
x=fread(fd1,N_records,'double');
x=single(x);


%x=x(N_records+1:end);

%N_records=N_records/2;
N_slice=N_records/NPART/6;

xx=zeros(6,NPART,round(N_slice));

% for n=1:N_slice
%     for j=1:6
%         for i=1:NPART
%             xx(j,i,n)=x(NPART*6*(n-1)+NPART*(j-1)+i);
%         end
%     end
% end

for n=1:N_slice
    for j=1:6
        for i=1:NPART
            xx(j,i,n)=x(NPART*6*(2*n-1)+NPART*(j-1)+i);
        end
    end
end

break

    clear x

% Variable    E         Ph      X       Y       X'      Y'
filter_value=[7785      0       0       0       0       0];
do_normalize=[0         0       0       0       0       0];

corr_value  =[0         0       0       0       0       0];
do_correct  =[0         0       0       0       0       0];

xx1=xx;
for i=1:6
index{i}=reshape(xx(i,:,:)~=filter_value(i),[NPART,N_slice]);
V=reshape(xx(i,:,:),[NPART,N_slice]);
Av_V(i)=mean(V(index{i}));
    if do_normalize(i)
        xx1(i,:,:)=V-index{i}.*Av_V(i);
    elseif do_correct(i)~=0
        xx1(i,:,:)=V-index{i}.*do_correct(i);
    end
end

X_mean=Av_V(3);
Y_mean=Av_V(4);
Xm_mean=Av_V(5);
Ym_mean=Av_V(6);


for n=1:N_slice
    for j=1:6
        for i=1:NPART
            x1(NPART*6*(n-1)+NPART*(j-1)+i)=xx1(j,i,n);
        end
    end
end

 fd2=fopen(nm2,'wb');
 fwrite(fd2,x1','double');

fclose all;
clear all;