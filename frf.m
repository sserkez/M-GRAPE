nm0='geninp_U1';
for ii=1:3
	fid = fopen(['s_2.' int2str(ii)], 'w');
	fd=fopen(nm0,'r');
    nm='U2';
	while 1
    		tline = fgetl(fd);
    		if ~ischar(tline), break, end
    		if findstr('ipseed',tline)
			fprintf(fid,' ipseed = %d\n',33*ii+1);
			continue
		end
    		if findstr('outputfile',tline)
			fprintf(fid,' outputfile =''U2.%d.out''\n',nm,ii);
			continue
		end
    		if findstr('fieldfile',tline)
%			fprintf(fid,' fieldfile =''%s.%d.6um.dat''\n',nm,ii);
%			fprintf(fid,' fieldfile =''%s.usr2.%d.out.dfl2''\n',nm,ii);
%			fprintf(fid,' fieldfile =''%s.soft.%d.3.105um.dat''\n',nm,ii);
			fprintf(fid,' fieldfile =''U1.%d.dat''\n',nm,ii);
			continue
		end
    		if findstr('partfile',tline)
			fprintf(fid,' partfile =''%s.%d.out.dpa''\n',nm,ii);
%			fprintf(fid,' partfile =''%s.15AQ8.%d.out.dpa''\n',nm,ii);
			continue
        end
            if findstr('distfile',tline)
			fprintf(fid,' distfile =''%s.%d.out.dpa''\n',nm,ii);
%			fprintf(fid,' partfile =''%s.15AQ8.%d.out.dpa''\n',nm,ii);
			continue
        end
            if findstr('beamfile',tline)
			fprintf(fid,' beamfile =''sase2.20pC.bds''\n');
%			fprintf(fid,' beamfile =''%s.mono.%d.out.bds''\n',nm,ii);
			continue
		end
		fprintf(fid,'%s\n',tline);
	end
	fclose(fid);
	fclose(fd);
end
