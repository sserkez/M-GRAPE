% find optimal xlamds for the steady-state calculations
clear all
nm='c:\-D-\Work\LCLS\Thesis_source\geninp_700_u1';
nm_out='c:\-D-\Work\LCLS\Thesis_source\700_u1_1.out';
%copyfile('300.in','geninp')
xlamds_arr=1.77e-9:0.002e-9:1.78e-9;

P(1:numel(xlamds_arr))=0;
for k=1:numel(xlamds_arr)
fclose all;
delete([nm_out,'.dfl']);
fid = fopen(['c:\-D-\Work\LCLS\Thesis_source\geninp'],'w');
fd = fopen(nm,'r');

	%fd=fopen(nm0,'r');
	while 1
    		tline = fgetl(fd);
    		if ~ischar(tline), break, end
    		if findstr('xlamds',tline)
			fprintf(fid,' xlamds = %d\n',xlamds_arr(k));
            else
                fprintf(fid,[tline,'\n']);
			continue
            end
    end
fclose all;

dos('genesis301');
%             Xn=fieldimport('prop.out.dfl',M,1);
%     %        [H{13}]=fieldplot(13,Xn,leng_u2,'amplified field');
%             P(k)=sum(sum(abs(Xn).^2));
d=outread(nm_out,1);
P(k)=d.outp.power.v(end);

[X,N]=fieldimport_all([nm_out,'.dfl'],d.inp.ncar,1);
[H{15}]=fieldplot3d(15,X,1,0,d.inp.xlamds,'field',1);

figure(58874)
plot(xlamds_arr,P);

end

