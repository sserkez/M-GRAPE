function par=parread(nm_p)

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
par.E=reshape(xx(1,:,:),[npart N_slice]);
par.Ph=reshape(xx(2,:,:),[npart N_slice]);
par.X=reshape(xx(3,:,:),[npart N_slice]);
par.Y=reshape(xx(4,:,:),[npart N_slice]);
par.XP=reshape(xx(5,:,:),[npart N_slice]);
par.YP=reshape(xx(6,:,:),[npart N_slice]);
clear xx

return