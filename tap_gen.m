%generate taper profile
AW=5.2563;
QF=8.5;
nsec=15;

Nwig=74;
NwigD=16;
minus=1;
QFwig=6;
QFwigD=84;

% z_start=30;
% z_end=(74+16)*0.068*nsec;
% a=0.0005;
% b=1.3;

z_start=25;
z_end=(74+16)*0.068*nsec;
a=0.001;
b=1.3;

% z_start=29;
% z_end=(74+16)*0.068*nsec;
% a=0.0007;
% b=1.3;

% z_start=25;
% z_end=(74+16)*0.068*nsec;
% a=0.0007;
% b=1.4;

z_start=25;
z_end=(74+16)*0.068*nsec;
a=0.0004;
b=1.3;

Z=0:(74+16)*0.068:z_end;
K=AW-a*(Z-z_start).^b;
K(Z<=z_start)=AW;
figure(58586)
hold all
plot(Z,K,'linewidth',2)

figure(58587)

plot(Z,K,'linewidth',2)

disp(z_end);



 fid = fopen('C:\-D-\Work\SASE3_SXRSS\tdp_3_tap\taper_stst.lat','w');
 fprintf(fid,'# header is included \n');
 fprintf(fid,'? VERSION= 1.00  including new format \n');
 fprintf(fid,'? UNITLENGTH= %g :unit length in header \n',0.068);
 % 1:a sections of wiggler
fprintf(fid,'AW      %g   %g  %g \n', K(1), Nwig, 0);
for i =2 :1:numel(Z)
    fprintf(fid,'AW      %g   %g  %g \n', K(i), Nwig, NwigD);
%    fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, 0);
end

for i =1 :1:ceil(numel(Z)/2)+1
    fprintf(fid,'QF       %g   %g  %g \n', minus*QF, QFwig, QFwigD);
    fprintf(fid,'QF      %g   %g  %g \n', -minus*QF, QFwig, QFwigD);
%    fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, 0);
end

fprintf(fid,'QX      0   3600  0 \n');
fprintf(fid,'QY      0   3600  0 \n');

for i =2 :1:numel(Z)
    fprintf(fid,'AD      %g   %g  %g \n', K(i), NwigD, Nwig);
%    fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, 0);
end


fprintf(fid,'SL      0.0000E+00  3600    0\n');
fprintf(fid,'CX      0.0000E+00  3600    0\n');
fprintf(fid,'CY      0.0000E+00  3600    0\n');
fprintf(fid,'AX      0.0000E+00  3600    0\n');
fprintf(fid,'AY      0.0000E+00  3600    0\n');
fclose all;

%%
cd('c:\-D-\Work\SASE3_SXRSS\tdp_3_tap\');
dos('c:\-D-\Work\SASE3_SXRSS\tdp_3_tap\genesis301.exe');
%%
nm_p='c:\-D-\Work\SASE3_SXRSS\tdp_3_tap\Taper_stst.out';
d(2)=outread(nm_p,1);
H{3}=outpot_e(3,d(2));
H{4}=outpot_ph(4,d(2));

Pz_s=d(2).outp.power.v(:,:);
Zscale=d(2).outp.Zscale;

figure(6875)
hold all
plot(Zscale,Pz_s,'linewidth',2);