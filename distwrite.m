function distwrite(dist,filename);

%filename='C:\-D-\WORK\LCLS\tmp\dist_500_1.3.dat';

dlmwrite(filename, '? VERSION = 1.0 ', 'delimiter', '');
dlmwrite(filename, ['? SIZE = ',num2str(size(dist.T,1))], 'delimiter', '','-append');
dlmwrite(filename, ['? CHARGE = ',num2str(dist.charge)], 'delimiter', '','-append');
dlmwrite(filename, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');

%dlmwrite(filename2,[X_val-4e-5 X_mom-mean(X_mom) Y_val Y_mom T_val+6.1e-14 Gamma],'delimiter',' ','precision',8,'-append');
%dlmwrite(filename2,[X PX Y PY T-min(T) Gamma],'delimiter',' ','precision',8,'-append');
dlmwrite(filename,[dist.X dist.PX dist.Y dist.PY dist.T dist.G],'delimiter',' ','precision',8,'-append');

clearvars -except I Ti
fclose all;