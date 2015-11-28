D0=1123; %[l/mm]
D1=1.6; %[l/mm2]
D2=0.002; %[l/mm3]

% gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)         !!!!!!
% gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) !!
gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)         !!!!!!
gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k) !!

%R_g_tang=195;           %tangential radius of grating     !!!!!!
R_g_tang=185;           %tangential radius of grating new according to PRD !!!!!!
R_g_sag=0.18;           %sagital radius of grating     !!!!!

U1_to_G=1.26;           %U1 to Grating distance in m !!!!!

R_g_tang1=R_g_tang-R_g_sag;
R_g_sag1=R_g_sag;

roll_angle=2*pi/360*0;
R_g_tang=1/((2*R_g_sag1-R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
R_g_sag=1/((2*R_g_sag1+R_g_tang1*cos(2*roll_angle)+R_g_tang1)/(2*R_g_sag1*(R_g_sag1+R_g_tang1)));
clear R_g_sag1 R_g_tang1

%Theta_g_i=1*0.01745;     %incidence angle         !!!!!!!
Theta_g_i=1.04*0.01745;     %incidence angle new according to PRD        !!!!!!!

R_m_tang=23.2;          %tangential radius of mirror M2 !!!
Theta_m_i=0.015;         %M2 incidence angle original!!!
%Theta_m_i=0.02;         

z3_t_offset=1.53;       %Grating to refocusing mirror distance original !!!
%z3_t_offset=1.59;       %Grating to refocusing mirror distance !!!

slitpos=1.350;          %slit position !!!

deltaU=3.87;
%deltaU=0;

G_to_U2=3.353;