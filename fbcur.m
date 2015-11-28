x=dlmread('sase2.20pc.bds','',2,0);
NN=800;
y=zeros(NN,9);
y(:,1)=linspace(x(1,1),x(end,1),NN);
y(:,2)=interp1(x(:,1),x(:,2),y(:,1));
y(:,3)=interp1(x(:,1),x(:,3),y(:,1));
y(:,4)=interp1(x(:,1),x(:,4),y(:,1));
y(:,5)=interp1(x(:,1),x(:,5),y(:,1));
y(:,6)=interp1(x(:,1),x(:,6),y(:,1));
y(:,7)=interp1(x(:,1),x(:,7),y(:,1));
y(:,8)=interp1(x(:,1),x(:,8),y(:,1));
y(:,9)=interp1(x(:,1),x(:,9),y(:,1));
for tt=1:1
%inf=['sase3.mono1p5A.' int2str(tt) '.out']
inf=['U1.1.out']
%inf1=['sase3.mono.' int2str(tt) '.out'];
outf=[inf '.bds'];
fdr=fopen(inf,'r')
fdw=fopen(outf,'w')
while 1; 
    tline = fgetl(fdr);
    if ~ischar(tline), break, end
    if findstr('$end',tline), break, end
    if findstr('delgam',tline),k=findstr('=',tline); tline(k)=' '; delgam = sscanf(tline,'%*s%g'); end
    if findstr('rxbeam',tline),k=findstr('=',tline); tline(k)=' '; rxbeam = sscanf(tline,'%*s%g'); end
    if findstr('rybeam',tline),k=findstr('=',tline); tline(k)=' '; rybeam = sscanf(tline,'%*s%g'); end
    if findstr('alphax',tline),k=findstr('=',tline); tline(k)=' '; alphax = sscanf(tline,'%*s%g'); end
    if findstr('alphay',tline),k=findstr('=',tline); tline(k)=' '; alphay = sscanf(tline,'%*s%g'); end
    if findstr('emitx',tline),k=findstr('=',tline); tline(k)=' '; emitx = sscanf(tline,'%*s%g'); end
    if findstr('emity',tline),k=findstr('=',tline); tline(k)=' '; emity = sscanf(tline,'%*s%g'); end
    if findstr('xbeam',tline),k=findstr('=',tline); tline(k)=' '; xbeam = sscanf(tline,'%*s%g'); end
    if findstr('ybeam',tline),k=findstr('=',tline); tline(k)=' '; ybeam = sscanf(tline,'%*s%g'); end
    if findstr('pxbeam',tline),k=findstr('=',tline); tline(k)=' '; pxbeam = sscanf(tline,'%*s%g'); end
    if findstr('pybeam',tline),k=findstr('=',tline); tline(k)=' '; pybeam = sscanf(tline,'%*s%g'); end
    if findstr('curlen',tline),k=findstr('=',tline); tline(k)=' '; curlen = sscanf(tline,'%*s%g'); end
    if findstr('nslice',tline),k=findstr('=',tline); tline(k)=' '; nslice = sscanf(tline,'%*s%g'); end
    if findstr('ntail',tline),k=findstr('=',tline); tline(k)=' '; ntail = sscanf(tline,'%*s%g'); end
    if findstr('gamma0',tline)
		k=findstr('=',tline); 
		tline(k)=' ';
		k=findstr('D',tline); 
		if ~isempty(k)
			tline(k)='E';
		end
		gamma0 = sscanf(tline,'%*s%g');
     end
   	if findstr('zsep',tline)
		k=findstr('=',tline); 
		tline(k)=' ';
		k=findstr('D',tline); 
		if ~isempty(k)
			tline(k)='E';
		end
		zsep = sscanf(tline,'%*s%g');
     	end
   	if findstr('xlamds',tline)
		k=findstr('=',tline); 
		tline(k)=' ';
		k=findstr('D',tline); 
		if ~isempty(k)
			tline(k)='E';
		end
		xlamds = sscanf(tline,'%*s%g');
     	end
%    disp(tline)
end
while 1; 
    tline = fgetl(fdr);
    if ~ischar(tline), break, end
    if findstr('=================',tline), break, end
    if findstr('entries per record',tline), NR = sscanf(tline,'%g'); end
    if findstr('history records',tline), NS = sscanf(tline,'%g'); end
end
fprintf(fdw,' ? SIZE=%d\n',nslice);
%fprintf(fdw,' ? COLUMNS ZPOS CURPEAK ELOSS\n');
%fprintf(fdw,' ? COLUMNS ZPOS BETAX BETAY CURPEAK\n');
fprintf(fdw,' ? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY CURPEAK ELOSS\n');
%fprintf(fdw,' ? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY CURPEAK\n');
%delta=5.9976*curlen/nslice
delta=zsep*nslice*xlamds/nslice;
for n=1:NS
    k=1;
    while 1; 
        tline = fgetl(fdr);
       if findstr('power',tline), k = 0; end
       if k==NR
%           xx=sscanf(tline,'%*g%*g%*g%*g%*g%*g%g%*g%*g%*g%*g%*g%*g%g');
           xx=sscanf(tline,'%*g%*g%*g%*g%*g%g%*g%*g%*g%*g%g');
%             fprintf(fdw,'%g %f %g %g\n',delta*(n+ntail-1),gamma0+xx(1),xx(2),curpeak);
%             fprintf(fdw,'%g %g\n',delta*(n+ntail-1),3e3);
 % fprintf(fdw,'%g %g %g\n', y(n,1),y(n,8), y(n,9));
%  fprintf(fdw,'%g %f %g %g %g %g  %g %g %g\n', delta*(n+ntail-1), gamma0+xx(1),xx(2),y(n,4), y(n,5), 15.6983,27.9257,y(n,14), y(n,15));
%  fprintf(fdw,'%g %g  %g %g %g\n', delta*(n+ntail-1), 15.4571,26.0768,y(n,14), y(n,15));
  %fprintf(fdw,'%g %f %g %g %g %g  %g %g %g\n', delta*(n+ntail-1), gamma0+xx(1),xx(2),y(n,4), y(n,5), 15.4571,27.5315,y(n,14), y(n,15));
  fprintf(fdw,'%g %f %g %g %g %g  %g %g %g\n', delta*(n+ntail-1), gamma0+xx(1),xx(2),y(n,4), y(n,5), y(n,6), y(n,7), y(n,8), y(n,9));
%  fprintf(fdw,'%g %f %g %g %g %g  %g %g %g\n', delta*(n+ntail-1), gamma0+xx(1),xx(2),y(n,4), y(n,5), 19.7821,30.6236,y(n,14), 0);
%  fprintf(fdw,'%g %g  %g %g \n', delta*(n+ntail-1), 19.7212, 29.59130,y(n,14));
%  fprintf(fdw,'%g %g \n', delta*(n+ntail-1),y(n,14));
%  fprintf(fdw,'%g %f %g %g %g %g  %g %g\n', delta*(n+ntail-1), gamma0+xx(1),xx(2),y(n,4), y(n,5), 19.7821,30.6236,y(n,14));
        end
        if ~ischar(tline), break, end
        if findstr('current',tline), curpeak = sscanf(tline,'%g'); end
        if findstr('=================',tline), break, end
        k=k+1;
    end

end
fprintf(fdw,'\n');
fclose(fdr);
fclose(fdw);
end
