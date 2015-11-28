% Optimal_call_script.m
% Claudio Emma 
% UCLA and SLAC
% October 2013

% This script runs genesis according to the optimal parameters
% obtained from J. Wu and J.Yiao's optimisation package

pathname0= pwd;
pathname1= '~/bin';
pathname2=[pathname0, '/genesis'];

% Define the input structure following optimization_STF_standard.m

% Now manually input the optimal values for the following variables
% they can be found in the results text file of the previous simulation

tapstart_op= 14.1061;
order_op= 1.95291;
ratio_op= 0.119707;
quad_ini_op = 26.2581;
quad_fin1_op = 15.2634;
quad_fin2_op = 39.9788;
quadchgstart1_op = 62.7211;
quadchgstart2_op = 195.211;

% Electron rest mass energy in eV

E00 = 0.51099906E6;

in = struct('E0', 13.64E9, 'DELE', 2.54404*E00, 'EMITX', 0.3E-6, ...
    'EMITY', 0.3E-6, 'CURPEAK', 4000, 'PRAD0', 5.0E6, 'XLAMD', ...
    0.032, 'AW0',  2.3832, 'NWIG', 106, 'NGAP', 32, 'NSEC', 50, ...
	    'LU', 200.0, 'IWITYP', 0, 'NPro', 9,'ITGAUS', 3);
    
        E0 = in.E0;
        DELE = in.DELE;
        EMITX = in.EMITX;
        EMITY = in.EMITY;
        CURPEAK = in.CURPEAK;
        PRAD0 = in.PRAD0;
        XLAMD = in.XLAMD;
        AW0 = in.AW0;
        NWIG = in.NWIG;
        NGAP = in.NGAP;
        NSEC = in.NSEC;
        LU = in.LU;
        IWITYP = in.IWITYP;
        NPro = in.NPro;
        ITGAUS = in.ITGAUS;
        disp(in);
        
        
        
        XLAMD= round(XLAMD*1000)/1000;
        
        DELZ= ceil((NWIG+NGAP)*NSEC/4000);
        
        if DELZ < 8 
    DL = ceil(8/DELZ)*DELZ;
    FL = DL;
else
    DL = DELZ;
    FL = DELZ;
end
NWIG = round(NWIG/DELZ)*DELZ;
NGAP = round(NGAP/DELZ)*DELZ;

disp(['DELZ is set to be ', num2str(DELZ)]);
disp(['DL = ', num2str(DL),' FL = ',num2str(FL)]);
disp(['NWIG = ',num2str(NWIG), ' NGAP = ',num2str(NGAP)]);





    AVERBETA = 30.0;  %average beta function along the undulator, unit: M
    XLAMDS = XLAMD*(1+AW0^2)/2/(E0/E00)^2; % radiation wavelength, unit : m
    ZRAYL = 10;         %Rayleigth length of the radiaton, unit: M
    ZWAIST = 0;        %Postion of the waist of the radiation, unit: M
    % belows are the magnetic parameters
%     FL = 8.0; % Focusing quadrupole length (undulator periods)
%     DL = 8.0; % Defocusing quadrupole length (undulator periods) 
    DRL = (2*(NWIG+NGAP)-FL-DL)/2; % Drifle length between adjacent quads (undulator periods)
    AWD = AW0; % Virtual undulator parameters of the gap
    F1ST = 0;
%     DELZ = 1; % the integration step for the genesis simulation

    % path name
    pathname0 = pwd;                % working path name
    % generate the directory to contain genesis files
    if exist('genesis','file') ~= 7
        mkdir('genesis');
    end
    pathname1 = '~/bin';  % where the 'genesis2.exe' is % jw
    pathname2 =  [pathname0,'/genesis'];   % where the 'genesis.cmd' and 'run.sh' are % jw
    
    % scan parameters
    scanpoint1 = 9; % Scan point number for each parameter scanning % jw
    iteration = 1 ; % 0, no iteration scan, 1, iteration scan
    scanpoint2 = 9; % scan point number for the iteration scan % jw
    %belows are initial set of the scan range for each parameter scan, they
    %will be changed automatically for different input parameters.
    ratio_min = 0.01;
    ratio_max  = 0.2; % jw
    tapstart_min = 0; 
    tapstart_max = 22; 
    order_min = 1.2; 
    order_max = 3; 

    quad_ini_min = 5; 
    quad_ini_max = 50; 
    quad_fin1_min = 1; 
    quad_fin1_max = 30; 
    quad_fin2_min = 10; 
    quad_fin2_max = 100; 
    quadchgstart1_min = 10; 
    quadchgstart1_max = 100;
    quadchgstart2_min = 50; 
    quadchgstart2_max = 200;
    
    % control parameters
    draw = 1 ; % to plot or not plot the parameter scanning curves


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set input parameters
    %Claudio Edit add ITGAUS as an input parameter and 

    input.ITGAUS = ITGAUS;

    % Claudio Edit manually tell GENESIS to dump field dist at the 
    % IPRADIth integration step
    
    input.IPRADI= 0;
    
    % Claudio Edit: Time dependence parameters
    
    ZSEP = floor((1*10^(-8))/XLAMDS);
    NSLICE = 10;
    
    input.ITDP = 1;
    input.NSLICE = NSLICE;
    input.ZSEP= ZSEP;
    input.NTAIL=-64;
    
 disp(['ZSEP is set to be ', num2str(ZSEP)]);
 disp(['XLAMDS is set to be ', num2str(XLAMDS)]);
 disp(['NSLICE is set to be ', num2str(NSLICE)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    input.pathname0 = pathname0;
    input.pathname1 = pathname1;
    input.pathname2 = pathname2;
    
    input.NGAP = NGAP;
    input.E0 =E0;
    input.DELE = DELE;
    input.LU = LU;
    input.XLAMD = XLAMD;
    input.NWIG  = NWIG ;
    input.NSEC = NSEC;
    input.FL = FL;
    input.DL = DL;
    input.DRL = DRL;
    input.DELZ = DELZ;
    input.F1ST = F1ST;
    input.AW0 = AW0;
    input.AWD = AWD;
    input.IWITYP = IWITYP;
    input.GAMMA0 = E0/E00;
    input.DELGAM = DELE/E00;
    input.EMITX = EMITX;
    input.EMITY = EMITY;
    input.AVERBETA = AVERBETA;
    input.CURPEAK = CURPEAK;
    input.XLAMDS = XLAMDS;
    input.ZRAYL = ZRAYL;
    input.ZWAIST = ZWAIST;
    input.PRAD0 = PRAD0;
    input.NPro = NPro;


     input.scanpoint1 = scanpoint1;
     input.iteration = iteration;
     input.scanpoint2 = scanpoint2;
     input.ratio_min = ratio_min;
     input.ratio_max = ratio_max;
     input.tapstart_min = tapstart_min;
     input.tapstart_max = tapstart_max;
     input.order_min = order_min;
     input.order_max = order_max;
     input.quad_fin1_min = quad_fin1_min;
     input.quad_fin1_max = quad_fin1_max;
     input.quad_ini_min = quad_ini_min;
     input.quad_ini_max = quad_ini_max;
     input.quad_fin2_min = quad_fin2_min;
     input.quad_fin2_max = quad_fin2_max;
     input.quadchgstart1_min = quadchgstart1_min;
     input.quadchgstart1_max = quadchgstart1_max;
     input.quadchgstart2_min = quadchgstart2_min;
     input.quadchgstart2_max = quadchgstart2_max;
     
     input.draw = draw;
 


  
% check the parameters
disp('Before scanning, first check parameters ...');
lfodo = 2*input.DRL + input.FL + input.DL;

% here we only consider one case, fodo cell length larger than a wiggler
% section length

if lfodo/input.NWIG < 2
    error('For this scanning, it requires FODO cell longer than two wiggler section!');
end

if mod(input.NWIG,1) ~=0 || mod(input.NGAP,1)~=0
    error('The NWIG and GAP is not integer, please change them.');
end

% 2, check the undulator length Lt > LU
Lt = input.NSEC*(input.NGAP+input.NWIG)*input.XLAMD;
input.Lt = Lt;
disp(['Total undulator length is ', num2str(Lt),'m']);
if Lt<input.LU
    input.LU = Lt;
    disp('The scanned undulator length is too long, set to total undulator length');
end


% 3, check the undulator type
if input.IWITYP ~=0 && input.IWITYP ~=1
    input.IWITYP = 0;
    disp('Wrong input for the undulator type, 0 for planar, 1 for helical');
end


% 4, check the up and low limit for the scanning parameters
if input.ratio_min > input.ratio_max
    input.ratio_min = 0.01;
    input.ratio_max = 0.5;
    warning('Wrong input for ratio scanning range, set to default values');
end

if input.tapstart_min > input.tapstart_max
    input.tapstart_min = 0;
    input.tapstart_max = 22;
    warning('Wrong input for tapering start point scanning range, set to default values');
end

if input.order_min > input.order_max
    input.order_min = 1.2;
    input.order_max = 3;
    warning('Wrong input for tapering profile order scanning range, set to default values');
end

if input.quad_fin1_min > input.quad_fin1_max
    input.quad_fin1_min = 1;
    input.quad_fin1_max = 40;
    warning('Wrong input for final quad gradient scanning range, set to default values');
end

if input.quad_fin2_min > input.quad_fin2_max
    input.quad_fin2_min = 10;
    input.quad_fin2_max = 100;
    warning('Wrong input for final quad gradient scanning range, set to default values');
end

if input.quad_ini_min > input.quad_ini_max
    input.quad_ini_min = 5;
    input.quad_ini_max = 50;
    warning('Wrong input for initial quad gradient scanning range, set to default values');
end

if input.quadchgstart1_min > input.quadchgstart1_max
    input.quadchgstart1_min = 10;
    input.quadchgstart1_max = 100;
    warning('Wrong input for quad gradient changing start point scanning range, set to default values');
end

if input.quadchgstart2_min > input.quadchgstart2_max
    input.quadchgstart2_min = 50;
    input.quadchgstart2_max = 200;
    warning('Wrong input for quad gradient changing start point scanning range, set to default values');
end

if input.quadchgstart2_max > input.LU 
    input.quadchgstart2_max = input.LU*0.9;
    warning('The quad gradient changing start point scanning maximum is too large, set to default values');
end



% main code
disp(input);


generate_taper_input_GUI_2(input,tapstart_op,order_op,ratio_op,...
quad_ini_op,quad_fin1_op,quad_fin2_op,quadchgstart1_op, quadchgstart2_op);
 
 
% cd(input.pathname2);
% 
% [a,b]=unix('./genesis.cmd');
% pause(1);
% 
% % [out]=read_genesis_mod_standard(input);
% % Since the above doesn't work read the gneesis input the pedestrian way so 
% % just copy the code from that function
% 
% 
% % to read the mod file
% % global input;
% cd(input.pathname2);
% fidin = fopen('mod.out','r');
% 
% 
% n=0;
% while ~feof(fidin)
%     tline = fgetl(fidin);
%     n=n+1;
%     if numel(tline)>6 && tline(5)=='z' && tline(6)=='[' && tline(7)=='m'
%         disp('hi')
%         n1 = n;
%     end
%     if numel(tline)>3 && tline(1)=='*' && tline(2)=='*' && tline(3)=='*' && tline(4)=='*'
%         n2 = n;
%     end
% end
% nl = n;
% 
% 
% data1 = [];
% data2 = [];
% fclose(fidin);
% fidin1 = fopen('mod.out','r');
% for i = 1:1:nl
%     tline1 = fgetl(fidin1);
%     if i>n1 & i<n2-1
%         tline2 = str2num(tline1);
%         data1 = [data1; tline2];
%     end
%     if i>n2+5
%         tline2 = str2num(tline1);
%         data2 = [data2; tline2];
%     end
% end
% fclose(fidin1);
% 
% 
% 
% out = [data1 data2];
% %temppower = out(:,4); % keep fundamental power jw
% %out(:,4) = out(:,24); % put 3rd harmonic power in the fundamental power slot jw
% %out(:,24) = temppower; % put fundamental power in the 3rd harmonic power slot jw
% %clear temppower; % jw
% clear data1;
% clear data2;
% 
% 
% cd(input.pathname0);
%     
%    [pmax0,I0] = max((out(:,4)));
%     L2 = out(I0,1);
% disp('******************************************************************');
% disp(' ');
% disp(' ');
% disp('After both the tapering and quad gradient variation optimization,');
% disp(['at L = ',num2str(L2),' m, the maximum power is ',num2str(pmax0/1e12),' TW']);
% disp(' ');
% disp(' ');
% disp('******************************************************************');
% 
% 
% cd(input.pathname0);
% 
% h = figure(9);
% 
% subplot(2,2,1); % power
% plot(out(:,1),out(:,4)/1e12);
% ylabel('Power (TW)');
% xlabel('Undulator length (m)');
% subplot(2,2,2); % transverse radiation beam size
% plot(out(:,1),out(:,8)*1e6);
% hold on;
% plot(out(:,1),sqrt(out(:,12).^2/2+out(:,13).^2/2)*1e6,'g')
% hold off;
% ylabel('\sigma_r_,_r (\mum) & \sigma_r_,_e ');
% xlabel('Undulator length (m)');
% subplot(2,2,3);
% plot(out(:,1),(out(:,4)*376.7/pi./out(:,8).^2).^(0.5));
% ylabel('Electric field E (V/m)');
% xlabel('Undulator length (m)');
% subplot(2,2,4);
% plot(out(:,1),out(:,11));
% ylabel('Bunching factor');
% xlabel('Undulator length (m)');
% 

return;





        