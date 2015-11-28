% get the current vs position info for wake calculation with 'wakecalc'

%%%%%%%%%%%%% current file for J. Qiang calculation %%%%%%%%%%%%% %%%%%%%%%%%%%%%%

%readOutput('C:\-D-\Work\LCLS\tmp\3\530_u1_tdp.out')
%%[current position] = genplot('current','profile');

% global info output
% figure
% plot(output.cur)
% current=fliplr(output.cur);
% position=output.t;%-info.zsep.*info.lambda;   %%% this is for J. Qiang calc
d=outread('C:\-D-\Work\LCLS\tmp\3\530_u1_tdp.out');
current=flipud(d.outp.current);
position=d.outp.Sscale;
size(current)
size(position)

beamfile=[position current];
size(beamfile)

fid2=fopen('C:\-D-\Work\LCLS\tmp\3\ebeams\150pC1p5kA3p5GeV.curr','w');
fprintf(fid2, '%12.10E\t %12.10E\t\n',beamfile');
fclose(fid2);
%clear current position beamfile



%%%%%%%%%%%% current file for S. Reiche calculation %%%%%%%%%%%%

%readOutput('C:\Users\gmarcus\Documents\GENESIS\LCLS-II\SXR\S2E\lcls2_sxr_S2E.out')
%[current position] = genplot('current','profile');
%global info output
%current=fliplr(output.cur);
%position=fliplr(output.t);%-info.zsep.*info.lambda/299792458;   %%% this is for J. Qiang calc
%size(current)
%size(position)

%beamfile=[position' current'];
%size(beamfile)

%fid2=fopen('current','w');
%fprintf(fid2, '%12.10E\t %12.10E\t\n',beamfile');
%fclose(fid2);
%clear current position beamfile
