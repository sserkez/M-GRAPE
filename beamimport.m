%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\-D-\Work\SASE3_chicane\tmp\sase2.20pC.bds
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2014/07/11 14:34:10

%% Initialize variables.
filename = 'C:\-D-\Work\SASE3_SXRSS\tdp_3\beam_12GeV_0.5nC.bds';
%filename = 'C:\-D-\Work\SASE3_SXRSS\tdp_3\beam_0.5nC.txt';
delimiter = ' ';
startRow = 5;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% %% Allocate imported array to column variable names
% ZPOS = dataArray{:, 1};
% GAMMA0 = dataArray{:, 2};
% DELGAM = dataArray{:, 3};
% EMITX = dataArray{:, 4};
% EMITY = dataArray{:, 5};
% BETAX = dataArray{:, 6};
% BETAY = dataArray{:, 7};
% CURPEAK = dataArray{:, 8};
% ELOSS = dataArray{:, 9};

%% Allocate imported array to column variable names
ZPOS = dataArray{:, 1};
GAMMA0 = dataArray{:, 2};
DELGAM = dataArray{:, 3};
EMITX = dataArray{:, 4};
EMITY = dataArray{:, 5};
BETAX = dataArray{:, 6};
BETAY = dataArray{:, 7};
XBEAM  = dataArray{:, 8};
YBEAM  = dataArray{:, 9};
PXBEAM  = dataArray{:, 10};
PYBEAM  = dataArray{:, 11};
ALPHAX  = dataArray{:, 12};
ALPHAY = dataArray{:, 13};
CURPEAK = dataArray{:, 14};
ELOSS = dataArray{:, 15};


ZPOS=(ZPOS-min(ZPOS))*1e6-3; %DELETE WHEN EXPORT!!!!!!!!!!!!!!!!!!!!!!!!!
%% Plot stuff
figure(5861) 

plot(ZPOS,CURPEAK,'linewidth',2,'color','k');
    ylabel('I [A]');
    xlabel('s [\mum]');
    xlim([0 65]);
    
figure(5862); plot(ZPOS,GAMMA0,'linewidth',2)
    ylabel('\gamma');
    xlabel('s [\mum]');
    xlim([0 65]);
    
figure(5863); 
plot(ZPOS,EMITX,'linewidth',2,'color','b','linestyle','-')
hold on
plot(ZPOS,EMITY,'linewidth',2,'color','r','linestyle','--')
hold off
    ylabel('\epsilon_{x,y} [\mum]');
    xlabel('s [\mum]');
    legend('x','y');
    xlim([0 65]);
    

    
figure(5864); plot(ZPOS,DELGAM,'linewidth',2);
    ylabel('\sigma_\gamma');
    xlabel('s [\mum]');
    xlim([0 65]);
    
    
figure(5865); plot(ZPOS,ELOSS,'linewidth',2,'color',[0 0.5 0]);
    ylabel('\gamma');
    xlabel('E [V/m]');
    xlim([0 65]);
    
    
for i=1:5
    set(figure(i+5860), 'Position', [-250+i*300, -140+i*150, 450, 350]);
end

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
break
%% write distribution
nm2='C:\-D-\Work\SASE3_SXRSS\tdp_3\beam_12GeV_0.5nC.bds'
nm2='C:\-D-\Work\SASE3_SXRSS\tdp_3\beam_12GeV_0.5nC_wake1.bds'
 %dlmwrite(nm2, '? VERSION = 1.0 ', 'delimiter', '');
% dlmwrite(nm2, ['? SIZE = ',num2str(size(T,1))], 'delimiter', '','-append');
% dlmwrite(nm2, ['? CHARGE = ' num2str(d.inp.charge)], 'delimiter', '','-append');
% dlmwrite(nm2, '? COLUMNS X XPRIME Y YPRIME T P  ', 'delimiter', '','-append');
% dlmwrite(nm2,[ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS],'delimiter',' ','precision',8,'-append');
dlmwrite(nm2,[ZPOS CURPEAK ELOSS],'delimiter',' ','precision',8,'-append');