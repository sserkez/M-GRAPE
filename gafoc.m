%% calculation of gaussian beam parameters after propagation through the lens            .  
%                     .
%                    /|\  
%    *-------z1------|||------------------d2_s-----------------*
%    s1              \|/                  d2_t              s2_s
%                     ' f_g_s, f_g_t                        s2_t

xlamds=1239.8411/energy*1e-9;  

s1=(36.5-energy*0.0118)*1e-6/2.35; %sourse sigma
z1=(1.787+2.37*energy/1000); %sourse distance
z1=(0.52+1.267+2.37*energy/1000); %sourse distance FIXED DUE TO PPT MONO_GAUSS DATA
%z1=2.2;
%z1=z1-1; %sourse distance corrected by 1 m
z1=z1+deltaU; %sourse distance with U8 removed
Theta_g_d=acos(cos(Theta_g_i)-xlamds/gr_sigma);%4.39 for 500 eV, 3.18 for 1keV !!!
%Theta_g_d=Theta_g_i; %!!!!!! remove !!!!!
w1=s1*2;

zr=pi*w1^2/xlamds;  %Raileygh distance


disp(' ');
disp('SLAC ');
disp(['E=',num2str(energy),'eV,   wav=',num2str(xlamds*1e9),'nm']);
disp(['z1=',num2str(z1),'m   s1=',num2str(s1*1e6),'mikr   zr=',num2str(zr),'m']);
disp(['w1=',num2str(w1*1e6)]);
%% sagital focal distance calculation
f_g_s=R_g_sag/(sin(Theta_g_i)+sin(Theta_g_d));
disp(['f_g_s=',num2str(f_g_s),'m']);
%% tangential focal distance calculation

% vls contribution
f_g_t_vls=gr_sigma^2*sin(Theta_g_d)^2/xlamds/gr_alpha; %!!! original!!!!
%f_g_t_vls=inf;  %!!!!!! remove !!!!!
%f_g_t_vls=gr_sigma^2/xlamds/gr_alpha;

% curvature contribution
f_g_t_c=R_g_tang/(1/Theta_g_d+Theta_g_i/Theta_g_d^2); %tangential focal distance due to curvature

% final tangential focal distance
f_g_t=1/(1/f_g_t_vls+1/f_g_t_c);
%f_g_t=1.07;

zr_prime=zr*(Theta_g_d/Theta_g_i)^2;
z1_prime=z1*(Theta_g_d/Theta_g_i)^2;
w1_prime=w1*(Theta_g_d/Theta_g_i);

disp(['f_g_t=',num2str(f_g_t),'m,  f_g_t_c=',num2str(f_g_t_c),'m,  f_g_t_vls=',num2str(f_g_t_vls),'m']);
%disp(['z1''=',num2str(z1_prime),'m,  s1''=',num2str(w1_prime/2*1e6),'mikr,   zr''=',num2str(zr_prime),'m']);
%% propagation 1

%sagital
z4_s=(f_g_s^2*(z1-f_g_s))/((z1-f_g_s)^2+zr^2)+f_g_s;%original!!!
z4_s=1/(1/f_g_s-1/(z1+zr^2/(z1-f_g_s)));
w4_s=w1/sqrt((1-z1/f_g_s)^2+1/f_g_s^2*zr^2);
%tangential
z2_t=(f_g_t^2*(z1_prime-f_g_t))/((z1_prime-f_g_t)^2+zr_prime^2)+f_g_t;
w2_t=w1_prime/sqrt((1-z1_prime/f_g_t)^2+1/f_g_t^2*zr_prime^2);
% w2_t_simplified=w1*Theta_g_i/Theta_g_d*f_g_t/sqrt(z1^2+zr^2);
% disp(num2str(w2_t/w2_t_simplified));
z2_t_geom=1/(1/f_g_t-1/z1_prime);

s4_s=w4_s/2;
s2_t=w2_t/2;

zr_slit=pi*w2_t^2/xlamds;

disp('-res-');
%disp(['z2_s=',num2str(z2_s),'m   w2_s=',num2str(w2_s*1e6),'  s2_s=',num2str(s2_s*1e6),'mikr   magn_s=',num2str(s2_s/s1)]);
disp(['z2_t_geom=',num2str(z2_t_geom),'m']);
disp(['z2_t=',num2str(z2_t),'m   w2_t=',num2str(w2_t*1e6),'  s2_t=',num2str(s2_t*1e6),'mikr   magn_s=',num2str(s2_t/s1)]);
disp(['zr_slit=',num2str(zr_slit),'m']);
%% tangential propagation after the slit

%M3pos=1.53;

%w3_t=sqrt(xlamds*6e-3/pi);
w3_t=w2_t;
%w3_t=3e-6;

zr_3=pi*w3_t^2/xlamds;
%z3_t=M3pos-z2_t;
%R_m_tang=60.0023;%25;

f_m_t=R_m_tang*sin(Theta_m_i)/2;
z3_t=z3_t_offset-z2_t;
%
% z3_t=0.18;
% f_m_t=0.174;

%z3_t=0.5;
%f_m_t=0.45;


z4_t=(f_m_t^2*(z3_t-f_m_t))/((z3_t-f_m_t)^2+zr_3^2)+f_m_t;
w4_t=w3_t/sqrt((1-z3_t/f_m_t)^2+1/f_m_t^2*zr_3^2);
s4_t=w4_t/2;

z4_t_geom=1/(1/f_m_t-1/(z3_t_offset-z2_t_geom));
z4_s_geom=1/(1/f_g_s-1/(z1));

disp(['z4_s=',num2str(z4_s-4.43),'m   w4_s=',num2str(w4_s*1e6),'  s4_s=',num2str(s4_s*1e6),'mikr  z4_s_r=',num2str(pi*w4_s^2/xlamds),'m']);
disp(['z4_t=',num2str(z2_t+z3_t+z4_t-4.43),'m   w4_t=',num2str(w4_t*1e6),'  s4_t=',num2str(s4_t*1e6),'mikr  z4_t_r=',num2str(pi*w4_t^2/xlamds),'m']);

%z3_t-f_m_t

%% number of illuminated grooves

% N_gr=w1*sqrt(1+(z1/zr)^2)/Theta_g_i/gr_sigma;
w_gr=w1*sqrt(1+(z1/zr)^2)/Theta_g_i;
N_gr=w_gr/gr_sigma;
disp(['Grating footprint (intensity fwhm)=',num2str(w_gr*1e6/(sqrt(2*log(2))))]); %to remove
%%
w_slit=w2_t*sqrt((z2_t-slitpos)^2/zr_slit^2+1);

sdx=8e-4; %resolution scale limit
snx=1000; %resolution scale points number
ds_arr=linspace(0,6e-6,30); %slit size array
res=linspace(0,0,size(ds_arr,1));
salpha=linspace(0,0,size(ds_arr,1));

%
%ds_arr=w_slit;
%ds_arr=2e-6;
ds_arr=0;
%

for index1=1:size(ds_arr,2)
   
%ds=10e-6;
ds=ds_arr(index1); %slit size
dlpl=linspace(-sdx,sdx,snx); %delta lambda per lambda
%dlpl=-2e-4:1e-6:2e-4;

distr_g=exp(-2.*(dlpl.*f_g_t.*xlamds./Theta_g_d./gr_sigma./w_slit).^2); %gauss distribution
distr_ft=double(abs(dlpl)<ds*Theta_g_d*gr_sigma/f_g_t/xlamds/2);         %mono distribution

if ds==0
    distr_eff=distr_g; % effective distribution (convolution)
else
    distr_eff=conv2(distr_g,distr_ft)/max(conv2(distr_g,distr_ft));
end

res_fwhm_g=findFWHM(linspace(-sdx,sdx,snx),distr_g);
res_fwhm_ft=ds*Theta_g_d*gr_sigma/f_g_t/xlamds;
res_fwhm_eff=findFWHM(linspace(-sdx*2,sdx*2,snx*2-1),distr_eff);

%%% tmp
%res_fwhm_eff=sqrt(2*log(2))*w_slit*Theta_g_d*gr_sigma/f_g_t/xlamds;
%%%


% if showpictures
%     figure(30);
%     clf
%     plot(dlpl,distr_g);
%     hold all
%     plot(dlpl,distr_ft);
%     if ds==0
%         plot(dlpl,distr_eff);
%     else
%         plot(linspace(-sdx*2,sdx*2,snx*2-1),distr_eff);
%     end
%     hold off
% end

%res(index1)=res_fwhm_g/res_fwhm_eff;
res(index1)=res_fwhm_eff;
%salpha(index1)=ds/w_slit;
salpha(index1)=ds;
%pause(0.5);
end

%!!! tmp
resolution_slit=res_fwhm_eff;

zzres_arg=salpha';
zzresfunc=(res.^-1)';
%resolution_slit=1/(f_g_t.*xlamds./Theta_g_d./gr_sigma./w_slit)
%resolution_slit=Theta_g_i*gr_sigma/pi/(w1*sqrt(1+(z1/zr)^2));
%!!! tmp
%1/resolution_slit
%resolution_slit=w_slit*Theta_g_d*gr_sigma./xlamds/f_g_t;

% figure(31);
% clf
% plot(salpha,res);
% V=axis;
% text(double(V(2)),double(V(4)),sprintf('\n resolution = %.3e  ', 1/resolution_slit),...
%        'HorizontalAlignment','right','VerticalAlignment','top','FontSize',10);
resolution_1=pi*w_gr/(1.18*gr_sigma);
resolution_slit=1/resolution_1;
%% Abberation

U2_res=2*Theta_g_d*gr_sigma*f_m_t/(pi*w4_t)/(z3_t_offset-f_m_t);
U2_res=1/U2_res;

ang_disp_U2=(z3_t_offset-f_m_t)/(Theta_g_d*gr_sigma*f_m_t);

C_30=(Theta_m_i^2/z3_t-Theta_m_i/R_m_tang)/2/z3_t+(Theta_m_i^2/z4_t-Theta_m_i/R_m_tang)/2/z4_t;
C_30*2*pi/xlamds*(8e-5)^3/Theta_m_i^3;