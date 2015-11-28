function generate_taper_input_GUI_2(input,tapstart, order, ratio, quad_ini, quad_fin1, quad_fin2, quadchgstart1, quadchgstart2)

% a general code to generate the genesis main input file and lattice file
% with tapering and varying beta functions along the undulator

% to generate the main input file, it is almost the same with the generate 
% SASE simulation, the only difference is that one should change the RXBEAM
% RYBEAM, ALPHAX and ALPHAY according to difference quad_ini
global E00 ;
E00 = 0.51099906E6; % unit: eV


QUADF = quad_ini;
QUADD = input.FL * quad_ini/input.DL;
NWIG = input.NWIG  ;
FL = input.FL ;
DL = input.DL ;
DRL = input.DRL ;
F1ST = input.F1ST ;
XLAMD = input.XLAMD ;
GAMMA0 = input.GAMMA0;
AW0 = input.AW0 ;
AWD = input.AWD ;
NSEC = input.NSEC ;

% calculate the beta function with QUADF and QUADD
% 1, calculate the fodo cell length, note that lfodo is still in units of
% XLAMD.
lfodo = 2*DRL + FL + DL;

% here we only consider one case, fodo cell length larger than a wiggler
% section length

if lfodo/NWIG < 1
    error('For LCLS, it should has FODO cell longer than a wiggler section!');
end
a1 =  floor(lfodo/NWIG); a1 =2;
ldrift = (lfodo - a1*NWIG)/a1;
if mod(ldrift,1)~=0
    error(['The drift length should be multiple of the unit length ',num2str(XLAMD)]);
end


% 2, calculate f, focal length of quad with thin lens approximation
brho = E00 * GAMMA0/1e9 / 0.29979;
f = 1/(QUADF/brho)/FL/XLAMD;
sinmu2 = lfodo*XLAMD/4/f;
sinmu = 2*sinmu2*sqrt(1-sinmu2^2);
cosmu = 1- 2*sinmu2^2;

if F1ST/FL == 1/2
    betax = lfodo*XLAMD*(1+sinmu2)/sinmu;
    betay = lfodo*XLAMD*(1-sinmu2)/sinmu;
    RXBEAM = sqrt(input.EMITX*betax/GAMMA0);
    RYBEAM = sqrt(input.EMITY*betay/GAMMA0);
    ALPHAX = 0;
    ALPHAY = 0;
elseif F1ST/FL == 0
    L1 = F1ST *XLAMD;
    L = lfodo*XLAMD/2;
    betax0 = L*(2*f^2 + f*(L-2*L1)+L1*(-L+L1))/f^2/sinmu;
    alfax0 = ((f^2-f*L+L*(-L+L1))/f^2-cosmu)/sinmu;
    betay0 = L*(2*f^2 - f*(L-2*L1)+L1*(-L+L1))/f^2/sinmu;
    alfay0 = ((f^2+f*L+L*(-L+L1))/f^2-cosmu)/sinmu;
    
    RXBEAM = sqrt(input.EMITX*betax0/GAMMA0);
    RYBEAM = sqrt(input.EMITY*betay0/GAMMA0);
    ALPHAX = alfax0;
    ALPHAY = alfay0;
else
    error('The value of F1ST should be set to half of FL or zero!');
end

RXBEAM, RYBEAM, ALPHAX, ALPHAY

[RXBEAM, RYBEAM, ALPHAX, ALPHAY]= cal_twiss(F1ST,brho,QUADF,XLAMD, input.FL,input.DL, input.DRL,input.EMITX,input.EMITY,GAMMA0)


%cd(input.pathname1);
fid = fopen('2wbaglineardetuned1tdp.in','w');

 fprintf(fid,'  $NEWRUN \n');
 fprintf(fid,'  VERSION= %g, \n', 1.0);
 fprintf(fid,'  NPART= 4096,\n ');
 fprintf(fid,' NBINS= 8,\n ');
 fprintf(fid,' XLAMD=  %g,\n ', input.XLAMD);
 fprintf(fid,' IPHSTY=1,\n ');
 fprintf(fid,' AW0=  %g,\n ', input.AW0);
 fprintf(fid,' AWD=  %g,\n ', input.AWD);
 if input.IWITYP == 1
     fprintf(fid,' XKX= 0.5,\n ');
    fprintf(fid,' XKY= 0.5,\n ');
     fprintf(fid,' IWITYP= 1,\n ');
 elseif input.IWITYP == 0
    fprintf(fid,' XKX= 0.0,\n ');
    fprintf(fid,' XKY= 1,\n ');
     fprintf(fid,' IWITYP= 0,\n ');
 else
     fprintf(fid,' XKX= 0.0,\n ');
     fprintf(fid,' XKY= 1,\n ');
      fprintf(fid,' IWITYP= 0,\n ');
     warning('Wrong input wiggler type, use planar undulator');
 end
 fprintf(fid,' wcoefz=  0.0,   0.0,   0.0,\n ');
 fprintf(fid,' CURPEAK=  %g,\n ', input.CURPEAK);
 fprintf(fid,' CURLEN= 1.27310E-6,\n ');
 fprintf(fid,' PRAD0=   %g,\n ', input.PRAD0);
 fprintf(fid,' GAMMA0=  %g,\n ', input.GAMMA0);
 fprintf(fid,' DELGAM=  %g,\n ', input.DELGAM);
 fprintf(fid,' IGAMGAUS= 1,\n ');
 fprintf(fid,' EMITX=  %g,\n ', input.EMITX);
 fprintf(fid,' EMITY=  %g,\n ', input.EMITY);
 fprintf(fid,' XLAMDS=  %g,\n ', input.XLAMDS);
 fprintf(fid,' ZRAYL=  %g,\n ', input.ZRAYL);
 fprintf(fid,' ZWAIST=  %g,\n ', input.ZWAIST);
 fprintf(fid,' NHARM= 5,\n '); % jw
 fprintf(fid,' iallharm = 1,\n '); % jw
 fprintf(fid,' LOUT= 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1,\n '); % jw
 fprintf(fid,' IDMPPAR= 0,\n ');
 fprintf(fid,' IDMPFLD= 0,\n ');
 fprintf(fid,' MAGINFILE = ''taperwbaglineartdp.lat'' ,\n');
 fprintf(fid,'  OUTPUTFILE= ''modwbaglineardetuned1.out'' ,\n ');
 fprintf(fid,' NWIG= %g,\n ', input.NWIG);
 fprintf(fid,' ZSTOP=  %g,\n ',input.LU);
 fprintf(fid,' NCAR= 191,\n ');
 fprintf(fid,' ISRSIG= 1,\n ');
 fprintf(fid,' ISRAVG= 1,\n ');
 fprintf(fid,' DELZ= %g,\n ', input.DELZ);
 fprintf(fid,' ZSEP= %g,\n ', input.DELZ*16);
 fprintf(fid,' NSEC= %g,\n ', input.NSEC);
 fprintf(fid,' RXBEAM=  %g ,\n ', RXBEAM);
 fprintf(fid,' RYBEAM=  %g ,\n ', RYBEAM);
 fprintf(fid,' ALPHAX=  %g,\n ', ALPHAX);
 fprintf(fid,' ALPHAY=  %g,\n ', ALPHAY);
 fprintf(fid,' XBEAM=  0.,\n ');
 fprintf(fid,' PXBEAM=  0.,\n ');
 fprintf(fid,' YBEAM=  0.,\n ');
 fprintf(fid,' PYBEAM=  0.,\n ');
 fprintf(fid,' QUADF= %g ,\n ',QUADF);
 fprintf(fid,' QUADD= %g ,\n ',QUADD);
 fprintf(fid,' FL= %g,\n ', input.FL);
 fprintf(fid,' DL= %g,\n ', input.DL);
 fprintf(fid,' DRL= %g,\n ',input.DRL);
 fprintf(fid,' F1ST= %g,\n ', input.F1ST);
 fprintf(fid,' ITDP= 1,\n ');
 %Claudio Edit
 fprintf(fid,' ITGAUS = %g ,\n',input.ITGAUS);
 fprintf(fid,' SHOTNOISE= 1 ,\n ');
 fprintf(fid,' NSLICE= 400 ,\n ');
 fprintf(fid,' NTAIL= 0 ,\n ');
 % Claudio Edit 
% fprintf(fid,' maginfile=''3.lat'',\n');
 fprintf(fid,' $end\n');
 fclose(fid);
 
 % to generate the lattice input file
 % 1, calculate the tapering start point with regard to the wiggler section
 a = floor(tapstart/XLAMD/(ldrift+NWIG)); % apply tapering from (a+1)th section
 b = floor((tapstart-a*XLAMD*(ldrift+NWIG))/XLAMD/input.DELZ); % apply tapering from (b+1)th period in (a+1)th section
 
 % calculate the coefficient of the taper profile
 L = input.Lt;
 if L <= tapstart
    % no taper
    disp('The start point of the taper is longer than the undulator length');
    disp('no taper');
    c = 0;
 end
 c = ratio/(L-tapstart)^order;
 
 
 % 2, calculate the coeficient of the quadrupole gradient decreasing
d1 = floor((quadchgstart1+F1ST*XLAMD)/((2*DRL + FL + DL)/2*XLAMD)); % from (d+1) quadrupole, decrease the gradient
d2 = floor((quadchgstart2+F1ST*XLAMD)/((2*DRL + FL + DL)/2*XLAMD)); % from (d+1) quadrupole, decrease the gradient
if quad_fin1  == quad_ini && quad_fin2 == quad_ini
    disp('no decreasing in quad gradients');
    g1 = 0;
    g2 = 0;
elseif L <= quadchgstart2
    disp('Start point to decrease quad gradient > undulator length');
    g1 = 0;
    g2 = 0;
else
    g1 = (1-quad_fin1/quad_ini)/(quadchgstart2-quadchgstart1);
    g2 = (1-quad_fin2/quad_fin1)/(L-quadchgstart2);
end

 XLAMD = XLAMD*input.DELZ;
 NWIG = NWIG/input.DELZ;
 F1ST = round(F1ST/input.DELZ);
 DL = DL/input.DELZ;
 FL = FL/input.DELZ;
 DRL = DRL/input.DELZ;
 ldrift = ldrift/input.DELZ;

 
 


% cd(input.pathname2);
 fid = fopen('taperwbaglineartdp.lat','w');
 fprintf(fid,'# header is included \n');
 fprintf(fid,'? VERSION= 1.00  including new format \n');
 fprintf(fid,'? UNITLENGTH= %g :unit length in header \n',XLAMD);
 % 1:a sections of wiggler
fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, 0);
for i =2 :1:a
    fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, ldrift);
%    fprintf(fid,'AW      %g   %g  %g \n', AW0, NWIG, 0);
end
% from (a+1)th section, starts to taper
b1= b;
if b1 > NWIG
    b1 = NWIG;
end

if b1 >1
    fprintf(fid,'AW      %g   %g  %g \n', AW0, 1, ldrift);
    for i =2 :1: b1
        fprintf(fid,'AW      %g   %g  %g \n', AW0, 1, 0);
    end
    for t=(b1+1):1:NWIG
        z = a*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
        if z <0
            z =0;
        end
        fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, 0);
    end
elseif b1 == 1
    fprintf(fid,'AW      %g   %g  %g \n', AW0, 1, ldrift);
    for t = 2:1:NWIG
        z = a*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
        if z < 0
            z = 0;
        end
        fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, 0);
    end
else
    t = 1; 
    z = a*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
    if z < 0
        z = 0;
    end
    fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, ldrift);
    for t = 2:1:NWIG
        z = a*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
        if z<0
            z = 0;
        end
        fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, 0);
    end
end
% from (a+2)th section to NSECth section
for i = (a+2):1:NSEC
    t = 1; 
    z = (i-1)*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
    fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, ldrift);
    for t = 2:1:NWIG
        z = (i-1)*(ldrift+NWIG)*XLAMD + t*XLAMD - tapstart;
        fprintf(fid,'AW      %g   %g  %g \n', AW0*(1-c*z^order), 1, 0);
    end
end

% then to define the Quads
if F1ST/FL == 1/2
    fprintf(fid,'QF      %g   %g  %g \n', quad_ini, F1ST, 0);
    for i = 2:1:d2
        temp = 1;
        if i>d1
            z = (2*DRL + FL + DL)/2*XLAMD*i - F1ST*XLAMD - quadchgstart1;
            temp = 1 - g1*z;
        end        
        
        if mod(i,2) == 0
            fprintf(fid,'QF      %g   %g  %g \n', quad_ini*(-1)*FL/DL*temp, DL, DRL);
        else
            fprintf(fid,'QF      %g   %g  %g \n', quad_ini*temp, FL, DRL);
        end
    end

    for i = d2+1:1:ceil(NSEC/a1*2)

            z = (2*DRL + FL + DL)/2*XLAMD*i - F1ST*XLAMD - quadchgstart2;
            temp = 1 - g2*z;
        
        
        if mod(i,2) == 0
            fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*(-1)*FL/DL*temp, DL, DRL);
        else
            fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*temp, FL, DRL);
        end
    end
    
    
    
    z = (2*DRL + FL + DL)/2*XLAMD*(ceil(NSEC/a1*2)+1) - F1ST*XLAMD - quadchgstart2;
    if mod(ceil(NSEC/a1*2)+1,2) == 0
        fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*FL/DL*(-1)*(1-g2*z), F1ST, DRL);
    else
        fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*(1-g2*z), F1ST, DRL);
    end
elseif F1ST == 0
    fprintf(fid,'QF      %g   %g  %g \n', quad_ini, FL, 0);
        for i = 2:1:d2
        temp = 1;
        if i>d1
            z = (2*DRL + FL + DL)/2*XLAMD*i - F1ST*XLAMD - quadchgstart1;
            temp = 1 - g1*z;
        end        
        
        if mod(i,2) == 0
            fprintf(fid,'QF      %g   %g  %g \n', quad_ini*(-1)*FL/DL*temp, DL, DRL);
        else
            fprintf(fid,'QF      %g   %g  %g \n', quad_ini*temp, FL, DRL);
        end
        end

        for i = d2+1:1:ceil(NSEC/a1*2)

            z = (2*DRL + FL + DL)/2*XLAMD*i - F1ST*XLAMD - quadchgstart2;
            temp = 1 - g2*z;
     
        
        if mod(i,2) == 0
            fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*(-1)*FL/DL*temp, DL, DRL);
        else
            fprintf(fid,'QF      %g   %g  %g \n', quad_fin1*temp, FL, DRL);
        end
        end        
end


 
 for i = 1:1:a

     % jw add - start
     if ldrift ==0
         AWDeq = AW0;
     else
      % NP is the radiation sllipage periods number, used to calucate the AWDeq
     NP = ceil(input.NGAP/(1+AW0^2));
     if NP == input.NGAP/(1+AW0^2)
         NP = NP +1;
     end
     % FLAG CLAUDIO EDIT OUT
     % NP = NP + 2/3; % jw
     AWDeq = sqrt((1+AW0^2-input.NGAP/NP)/(input.NGAP/NP));
     end
     fprintf(fid,'AD      %g   %g  %g \n', AWDeq, ldrift, NWIG);
     % jw add - end   
%     fprintf(fid,'AD      %g   %g  %g \n', AW0, ldrift, NWIG); % jw
 end
 
 for i = a+1:input.NSEC
     z = (i-1)*(ldrift+NWIG)*XLAMD + NWIG*XLAMD + ldrift/2*XLAMD - tapstart;
     if z < 0
         z =0;
     end
     AWD = AW0*(1-c*z^order);
     if ldrift ==0
         AWDeq = AWD;
     else
      % NP is the radiation sllipage periods number, used to calucate the AWDeq
     NP = ceil(input.NGAP/(1+AWD^2)); 
     if NP == input.NGAP/(1+AWD^2)
         NP = NP +1;
     end
     % FLAG CLAUDIO EDIT OUT
     % NP = NP + 2/3; % jw
     AWDeq = sqrt((1+AWD^2-input.NGAP/NP)/(input.NGAP/NP));
     end
     fprintf(fid,'AD      %g   %g  %g \n', AWDeq, ldrift, NWIG);
 end
 
 
%  
%  for i=a+1 : 1: input.NSEC
%      if ldrift ~= 0
%     t = 1; 
%     z = (i-1)*(ldrift+NWIG)*XLAMD + NWIG*XLAMD + t*XLAMD - tapstart;
%     if z < 0
%         z = 0;
%     end
%     fprintf(fid,'AD      %g   %g  %g \n', AW0*(1-c*z^order), 1, NWIG);
%     for t = 2:1:ldrift
%         z = (i-1)*(ldrift+NWIG)*XLAMD + NWIG*XLAMD + t*XLAMD - tapstart;
%         if z<0
%             z = 0;
%         end
%         fprintf(fid,'AD      %g   %g  %g \n', AW0*(1-c*z^order), 1, 0);
%     end
%      else
%              t = 0; 
%     z = (i-1)*(ldrift+NWIG)*XLAMD + NWIG*XLAMD + t*XLAMD - tapstart;
%     if z < 0
%         z = 0;
%     end
%     fprintf(fid,'AD      %g   %g  %g \n', AW0*(1-c*z^order), 0, NWIG);
%      end
%          
% end    
     
 
 fclose(fid);

%cd(input.pathname0);

function [RXBEAM, RYBEAM, ALPHAX, ALPHAY]= cal_twiss(F1ST,brho,QUADF,XLAMD, FL,DL, DRL,EMITX,EMITY,GAMMA0)
K = QUADF/brho;
Lqf = FL*XLAMD;
Lqd = DL * XLAMD;
Ld = DRL *XLAMD;

if F1ST == 0
    M1 = [cos(sqrt(K)*Lqf),1/sqrt(K)*sin(sqrt(K)*Lqf);-sqrt(K)*sin(sqrt(K)*Lqf),cos(sqrt(K)*Lqf)];
    M2 = [cosh(sqrt(K)*Lqd),1/sqrt(K)*sinh(sqrt(K)*Lqd);sqrt(K)*sinh(sqrt(K)*Lqd),cosh(sqrt(K)*Lqd)];
    M3 = [1, Ld; 0, 1];
    M = M3*M2*M3*M1;
    cosmu = (M(1,1)+M(2,2))/2
    sinmu = sqrt(1-cosmu^2);
    sinmu2 = sqrt(1/2-sinmu/2);
    mu2 = asin(sinmu2);
    
    betax0 = M(1,2)/sinmu;
    alfax0 = (M(1,1)-cosmu)/sinmu;
    
    M1y = [cos(sqrt(K)*Lqd),1/sqrt(K)*sin(sqrt(K)*Lqd);-sqrt(K)*sin(sqrt(K)*Lqd),cos(sqrt(K)*Lqd)];
    M2y = [cosh(sqrt(K)*Lqf),1/sqrt(K)*sinh(sqrt(K)*Lqf);sqrt(K)*sinh(sqrt(K)*Lqf),cosh(sqrt(K)*Lqf)];
    M3y = [1, Ld; 0, 1];
    My = M3y*M1y*M3y*M2y;
    cosmuy = (My(1,1)+My(2,2))/2;
    sinmuy = sqrt(1-cosmuy^2);
    betay0 = My(1,2)/sinmuy;
    alfay0 = (My(1,1)-cosmuy)/sinmuy;
    

    
    RXBEAM = sqrt(EMITX*betax0/GAMMA0);
    RYBEAM = sqrt(EMITY*betay0/GAMMA0);
    ALPHAX = alfax0;
    ALPHAY = alfay0;
end
    
return;
 
 
 
 
  
 
 
 