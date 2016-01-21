%curbeamexp
% reads current profile from given *.out file and writes it to beam file
% *.out.beam with ZPOS and CURPEAK columns
function make_curr_beam(nm_p)

d=outread(nm_p,1,0,2);

disp(' -Reading the current profile')
CURPEAK=d.outp.current;
ZPOS=d.outp.Sscale;
disp('  +done');

nm_beam=[nm_p,'.beam'];

disp([' -Writing the current profile to ',nm_beam])
dlmwrite(nm_beam, '? VERSION = 1.0 ', 'delimiter', '');
% dlmwrite(nm_beam, ['? SIZE = ',num2str(size(ZPOS,1))], 'delimiter', '','-append');
% dlmwrite(nm_beam, ['? CHARGE = ' num2str(d.inp.charge)], 'delimiter', '','-append');
dlmwrite(nm_beam, '? COLUMNS ZPOS CURPEAK  ', 'delimiter', '','-append');
dlmwrite(nm_beam,[ZPOS CURPEAK],'delimiter',' ','precision',8,'-append');
disp('  +done');