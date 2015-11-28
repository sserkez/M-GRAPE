%Xfel raprms 1

nnnn=-0.5;
% z3_t_offset=1.156;
% slitpos=0.944;
D0=1120; %[l/mm]
D1=2.12; %[l/mm2]
D1=1.43+0.01*nnnn-0.00; %[l/mm2]1.41
D2=0.002; %[l/mm3]

% gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)         !!!!!!
% gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) !!
gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)         !!!!!!
gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k) !!

%R_g_tang=195;           %tangential radius of grating     !!!!!!
R_g_tang=240+10*nnnn;%225           %tangential radius of grating new according to PRD !!!!!!
R_g_sag=0.4;           %sagital radius of grating     !!!!!
R_g_sag=0.42;           %sagital radius of grating     !!!!!
%R_g_sag=0.47;

U1_to_G=2.77;           %U1 to Grating distance in m !!!!!
%U1_to_G=0;

% grating roll effect calculations
% R_g_tang1=R_g_tang-R_g_sag;
% R_g_sag1=R_g_sag;
% 
% roll_angle=2*pi/360*0;
% R_g_tang=1/((2*R_g_sag1-R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
% R_g_sag=1/((2*R_g_sag1+R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
% clear R_g_sag1 R_g_tang1

%Theta_g_i=1*0.01745;     %incidence angle         !!!!!!!
Theta_g_i=deg2rad(1.00);     %incidence angle new according to PRD        !!!!!!!
nnn=-10;%3
R_m_tang=31-0.6*nnn;          %tangential radius of mirror M2 !!!
Theta_m_i=0.015;         %M2 incidence angle original!!!

% R_m_tang=31-0.6*nnn-5;          %tangential radius of mirror M2 !!!
% Theta_m_i=0.0175;         %M2 incidence angle original!!!

% %Theta_g_i=1*0.01745;     %incidence angle         !!!!!!!
% Theta_g_i=deg2rad(1.00);     %incidence angle new according to PRD        !!!!!!!
% nnn=-3;%3
% R_m_tang=31-0.6*nnn+2.5;          %tangential radius of mirror M2 !!!
% Theta_m_i=0.014;         %M2 incidence angle original!!!
% %Theta_m_i=0.02;

z3_t_offset=1.22;       %Grating to refocusing mirror distance original !!!
z3_t_offset=1.74-0.005*nnn+0.002;       %Grating to refocusing mirror distance 

slitpos=1.000+0.5;          %slit position !!!
sourcesize_x=1.0;
deltaU=6.1;
% deltaU=0;

G_to_U2=7.2-U1_to_G;




% %Xfel raprms 1
% 
% nnnn=-0.5;
% % z3_t_offset=1.156;
% % slitpos=0.944;
% D0=1120; %[l/mm]
% D1=2.12; %[l/mm2]
% D1=1.43+0.01*nnnn; %[l/mm2]1.41
% D2=0.002; %[l/mm3]
% 
% % gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)         !!!!!!
% % gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) !!
% gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)         !!!!!!
% gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k) !!
% 
% %R_g_tang=195;           %tangential radius of grating     !!!!!!
% R_g_tang=240+10*nnnn;%225           %tangential radius of grating new according to PRD !!!!!!
% R_g_sag=0.4;           %sagital radius of grating     !!!!!
% R_g_sag=0.42;           %sagital radius of grating     !!!!!
% %R_g_sag=0.47;
% 
% U1_to_G=2.77;           %U1 to Grating distance in m !!!!!
% %U1_to_G=0;
% 
% % grating roll effect calculations
% % R_g_tang1=R_g_tang-R_g_sag;
% % R_g_sag1=R_g_sag;
% % 
% % roll_angle=2*pi/360*0;
% % R_g_tang=1/((2*R_g_sag1-R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
% % R_g_sag=1/((2*R_g_sag1+R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
% % clear R_g_sag1 R_g_tang1
% 
% %Theta_g_i=1*0.01745;     %incidence angle         !!!!!!!
% Theta_g_i=deg2rad(1.00);     %incidence angle new according to PRD        !!!!!!!
% nnn=-10;%3
% R_m_tang=31-0.6*nnn;          %tangential radius of mirror M2 !!!
% Theta_m_i=0.015;         %M2 incidence angle original!!!
% %Theta_m_i=0.02;        
% 
% 
% % %Theta_g_i=1*0.01745;     %incidence angle         !!!!!!!
% % Theta_g_i=deg2rad(1.00);     %incidence angle new according to PRD        !!!!!!!
% % nnn=-3;%3
% % R_m_tang=31-0.6*nnn+2.5;          %tangential radius of mirror M2 !!!
% % Theta_m_i=0.014;         %M2 incidence angle original!!!
% % %Theta_m_i=0.02;
% 
% z3_t_offset=1.22;       %Grating to refocusing mirror distance original !!!
% z3_t_offset=1.74-0.005*nnn+0.004;       %Grating to refocusing mirror distance 
% 
% slitpos=1.000+0.5;          %slit position !!!
% sourcesize_x=0.9;
% deltaU=6.1;
% %deltaU=0;
% 
% G_to_U2=7.2-U1_to_G;