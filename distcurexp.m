function distcurexp(dist,filename)

I=hist(dist.T,1000);
Ti=linspace(min(dist.T),max(dist.T),1000);
I=I./sum(I.*(Ti(2)-Ti(1))).*dist.charge;
figure
plot(Ti,I);

beamfile=[Ti'.*3e8 I'];
size(beamfile)
if ~isempty(filename)
    fid=fopen(filename,'w');
    fprintf(fid, '%12.10E\t %12.10E\t\n',beamfile');
end
fclose(fid);