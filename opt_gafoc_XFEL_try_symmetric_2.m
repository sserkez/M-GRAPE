%optimization script
showpictures=1;
showpublicationpictures=0;
ienergy=300:100:1200;
%ienergy=900;
xlamdsout=1239.8411./ienergy.*1e-9;

% gr_alpha=1.79e-6;    %alpha parameter of grating (old=1,8)
% R_g_tang=140;         %tangential radius of grating
% R_g_sag=0.23;         %sagital radius of grating
% R_m_tang=26.7;
% 
% z3_t_offset=1.156;
% slitpos=0.944;
D0=1120; %[l/mm]
D1=2.12; %[l/mm2]
D1=1.43; %[l/mm2]1.41
D2=0.003; %[l/mm3]

% gr_sigma=8.9047e-7;   %sigma parameter of grating (1/k)         !!!!!!
% gr_alpha=1.269e-6;    %alpha parameter of grating (old=1,8) (n1/k/k) !!
gr_sigma=1/(D0*1e3);          %sigma parameter of grating (1/k)         !!!!!!
gr_alpha=D1*1e6/(D0*1e3)^2;       %alpha parameter of grating (old=1,8) (n1/k/k) !!

%R_g_tang=195;           %tangential radius of grating     !!!!!!
R_g_tang=240;%225           %tangential radius of grating new according to PRD !!!!!!
 R_g_sag=0.4;           %sagital radius of grating     !!!!!
R_g_sag=0.4;           %sagital radius of grating     !!!!!
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

R_m_tang=41.5-5.1;          %tangential radius of mirror M2 !!!
Theta_m_i=0.015;         %M2 incidence angle original!!!
%Theta_m_i=0.02;         

z3_t_offset=1.22;       %Grating to refocusing mirror distance original !!!
z3_t_offset=1.22+0.62-0.055;       %Grating to refocusing mirror distance 

slitpos=1.000+0.5;          %slit position !!!

deltaU=6.1;
%deltaU=0;

G_to_U2=7.2-U1_to_G;

z1_out=linspace(1,1,size(ienergy,2));
zr_slit_out=linspace(1,1,size(ienergy,2));
z2_t_out=linspace(1,1,size(ienergy,2));
z2_t_geom_out=linspace(1,1,size(ienergy,2));
w2_t_out=linspace(1,1,size(ienergy,2));
z4_s_out=linspace(1,1,size(ienergy,2));
z4_t_out=linspace(1,1,size(ienergy,2));
w4_s_out=linspace(1,1,size(ienergy,2));
w4_t_out=linspace(1,1,size(ienergy,2));
w1_out=linspace(1,1,size(ienergy,2));
N_gr_out=linspace(1,1,size(ienergy,2));
Theta_g_d_out=linspace(1,1,size(ienergy,2));
f_g_t_out=linspace(1,1,size(ienergy,2));
w_slit_out=linspace(1,1,size(ienergy,2));
resolution_slit_out=linspace(1,1,size(ienergy,2));

for index=1:size(ienergy,2)
    energy=ienergy(index);
    gafoc_XFEL;
    z1_out(index)=z1;
    z2_t_out(index)=z2_t;
    z2_t_geom_out(index)=z2_t_geom;
    w2_t_out(index)=w2_t;
    zr_slit_out(index)=zr_slit;
    z4_s_out(index)=z4_s-G_to_U2;
    z4_s_geom_out(index)=z4_s_geom-G_to_U2;
    z4_t_out(index)=z3_t_offset+z4_t-G_to_U2;
    z4_t_geom_out(index)=z3_t_offset+z4_t_geom-G_to_U2;
    w4_s_out(index)=w4_s;
    w4_t_out(index)=w4_t;
    w1_out(index)=w1;
    N_gr_out(index)=N_gr;

    w_slit_out(index)=w_slit;
    resolution_slit_out(index)=resolution_slit;
end

if showpictures
    
    figure(20);
    clf
    subplot(3,1,1);
    %title(['R=',num2str(R_g_tang),'m  D0=',num2str(D0),' l/mm  D1=',num2str(D1),' l/mm2']);
    title('Slit position');
    hold all
    plot(ienergy,z2_t_out,'linewidth',2,'color','b');
    plot(ienergy,z2_t_out-zr_slit_out,'color','r');
    plot(ienergy,z2_t_out+zr_slit_out,'color','g');
    plot(ienergy,linspace(slitpos,slitpos,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    plot(ienergy,z2_t_geom_out,'marker','*','color','m','MarkerSize',7);

    hold off
    xlabel('energy, eV');
    ylabel('Waist position [m]');
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    subplot(3,1,2);
    hold all
    plot(ienergy,1e6.*w2_t_out./2.*2.35);
    hold off
    xlabel('energy, eV');
    ylabel('Waist size FWHM [um]');
    subplot(3,1,3);
    hold all
    plot(ienergy,1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1));
    hold off
    xlabel('energy, eV');
    ylabel('Spot at slit size FWHM [um]');


    %w_slit=w2_t_out.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1);
    figure(21);
    plot(ienergy,1./resolution_slit_out);
    hold off
    xlabel('energy, eV');
    ylabel('Resolution on slit');

    figure(22);
    clf
    hold all
    plot(ienergy,N_gr_out);
    hold off
    xlabel('energy, eV');
    ylabel('Illuminated grooves');

    figure(24);
    clf
    subplot(3,1,1);
    title('Reimage (undulator entrance)');
    hold all
    plot(ienergy,z4_s_out,'linestyle','-','linewidth',2,'color','b');
    plot(ienergy,z4_s_out-pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','r');
    plot(ienergy,z4_s_out+pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','g');
%     plot(ienergy,z4_s_geom_out,'marker','*','linestyle','-','color','m','MarkerSize',7);

    plot(ienergy,z4_t_out,'linestyle','--','linewidth',2,'color','b');
    plot(ienergy,z4_t_out-pi*w4_t_out.^2./xlamdsout,'linestyle','--','color','r');
    plot(ienergy,z4_t_out+pi*w4_t_out.^2./xlamdsout,'linestyle','--','color','g');
%     plot(ienergy,z4_t_geom_out,'marker','*','linestyle','--','color','m','MarkerSize',7);

    plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    %plot(ienergy,z4_t_out,ienergy,z4_t_out-pi*w4_t_out.^2./xlamdsout,ienergy,z4_t_out+pi*w4_t_out.^2./xlamdsout,'linestyle','--');
    hold off
    xlabel('energy, eV');
    ylabel('Sagital (-) & tangential (--) waist position [m]');
    % subplot(2,2,2);
    % hold all
    % plot(ienergy,z4_t_out-4.43,ienergy,z4_t_out-4.43-pi*w4_t_out.^2./xlamdsout,ienergy,z4_t_out-4.43+pi*w4_t_out.^2./xlamdsout,'linestyle','--');
    % hold off
    % xlabel('energy, eV');
    % ylabel('Tangential waist position');
    subplot(3,1,2);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35,'linestyle','-','color','r');
    plot(ienergy,1e6.*w4_t_out./2.*2.35,'linestyle','--','color','r');
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    xlabel('energy, eV');
    ylabel('Waist size FWHM [um] s(-), t(--), orig(*)');

    subplot(3,1,3);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-','color','r');
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','--','color','r');
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-1)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-','color','b');
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-1)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','--','color','b');
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-2)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-','color','g');
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-2)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','--','color','g');
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    xlabel('energy, eV');
    ylabel('Spot FWHM [um] at 0m-r 1m-b 2m-g');

end

if showpublicationpictures

    
    figure(19);
    clf
    subplot(2,1,1);
    %title('Source position');
    hold all
    plot(ienergy,-z1_out+U1_to_G,'linewidth',2,'color','b');
    plot(ienergy,-z1_out+U1_to_G-pi*w1_out.^2./xlamdsout,'color','r');
    plot(ienergy,-z1_out+U1_to_G+pi*w1_out.^2./xlamdsout,'color','g');
    plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    hold off
    xlabel('Energy [eV]');
    ylabel('Waist position [m]');
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w1_out./2.*2.35,'linestyle','-','linewidth',2,'color','r');
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Energy [eV]');
    ylabel('Waist size [\mum]');
    
    
    figure(20);
    clf
    subplot(2,1,1);
    %title(['R=',num2str(R_g_tang),'m  D0=',num2str(D0),' l/mm  D1=',num2str(D1),' l/mm2']);
    %title('Slit position');
    hold all
    plot(ienergy,z2_t_out,'linewidth',2,'color','b');
    plot(ienergy,z2_t_out-zr_slit_out,'color','r');
    plot(ienergy,z2_t_out+zr_slit_out,'color','g');
    plot(ienergy,linspace(slitpos,slitpos,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    plot(ienergy,z2_t_geom_out,'marker','*','color','m','MarkerSize',7);

    hold off
    %xlabel('Energy [eV]');
    ylabel('Waist position [m]');
    xlabel('Energy [eV]');
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w2_t_out./2.*2.35,'linestyle','-','linewidth',2,'color','r');
    plot(ienergy,1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1),'linestyle','--','linewidth',1.5);
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Energy [eV]');
    ylabel('Size [\mum]');
    legend('Waist size','Spot at slit position width')
%     subplot(3,1,3);
%     hold all
%     plot(ienergy,1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1));
%     hold off
%     xlabel('Energy [eV]');
%     ylabel('Spot at slit size _{FWHM} [\mum]');


    %w_slit=w2_t_out.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1);
%     figure(21);
%     plot(ienergy,1./resolution_slit_out);
%     hold off
%     xlabel('Energy [eV]');
%     ylabel('Resolution on slit');
% 
%     figure(22);
%     clf
%     hold all
%     plot(ienergy,N_gr_out);
%     hold off
%     xlabel('Energy [eV]');
%     ylabel('Illuminated grooves');

    
    figure(241);
    clf
    subplot(2,1,1);
    hold all
    plot(ienergy,z4_s_out,'linestyle','-','linewidth',2,'color','b');
    plot(ienergy,z4_s_out-pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','r');
    plot(ienergy,z4_s_out+pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','g');
    plot(ienergy,z4_s_geom_out,'marker','*','linestyle','-','color','m','MarkerSize',7);
    plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([-1 0.5]);
    ylabel('Waist position [m]');
    xlabel('Energy [eV]');
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35,'linestyle','-','color','r','linewidth',2);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([5 35]);
    ylabel('Waist size [\mum]');
    xlabel('Energy [eV]');
    
    figure(242);
    clf
    subplot(2,1,1);
    hold all
    plot(ienergy,z4_t_out,'linestyle','-','linewidth',2,'color','b');
    plot(ienergy,z4_t_out-pi*w4_t_out.^2./xlamdsout,'linestyle','-','color','r');
    plot(ienergy,z4_t_out+pi*w4_t_out.^2./xlamdsout,'linestyle','-','color','g');
    plot(ienergy,z4_t_geom_out,'marker','*','linestyle','-','color','m','MarkerSize',7);
    plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([-1 0.5]);
    ylabel('Waist position [m]');
    xlabel('Energy [eV]');
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w4_t_out./2.*2.35,'linestyle','-','color','r','linewidth',2);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([5 35]);
    ylabel('Waist size [\mum]');
    xlabel('Energy [eV]');
    
%     figure(24);
%     clf
%     subplot(3,1,1);
%     %title('Reimage (undulator entrance)');
%     hold all
%     plot(ienergy,z4_s_out,'linestyle','-','linewidth',2,'color','b');
%     plot(ienergy,z4_s_out-pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','r');
%     plot(ienergy,z4_s_out+pi*w4_s_out.^2./xlamdsout,'linestyle','-','color','g');
%     plot(ienergy,z4_s_geom_out,'marker','*','linestyle','-','color','m','MarkerSize',7);
% 
%     plot(ienergy,z4_t_out,'linestyle','--','linewidth',2,'color','b');
%     plot(ienergy,z4_t_out-pi*w4_t_out.^2./xlamdsout,'linestyle','--','color','r');
%     plot(ienergy,z4_t_out+pi*w4_t_out.^2./xlamdsout,'linestyle','--','color','g');
%     plot(ienergy,z4_t_geom_out,'marker','*','linestyle','--','color','m','MarkerSize',7);
% 
%     plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
%     %plot(ienergy,z4_t_out,ienergy,z4_t_out-pi*w4_t_out.^2./xlamdsout,ienergy,z4_t_out+pi*w4_t_out.^2./xlamdsout,'linestyle','--');
%     hold off
%     xlabel('Energy [eV]');
%     ylabel('Sagital (-) & tangential (--) waist position [m]');
%     % subplot(2,2,2);
%     % hold all
%     % plot(ienergy,z4_t_out-4.43,ienergy,z4_t_out-4.43-pi*w4_t_out.^2./xlamdsout,ienergy,z4_t_out-4.43+pi*w4_t_out.^2./xlamdsout,'linestyle','--');
%     % hold off
%     % xlabel('energy, eV');
%     % ylabel('Tangential waist position');
%     subplot(3,1,2);
%     hold all
%     plot(ienergy,1e6.*w4_s_out./2.*2.35,'linestyle','-','color','r');
%     plot(ienergy,1e6.*w4_t_out./2.*2.35,'linestyle','--','color','r');
%     plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
%     hold off
%     xlabel('Energy [eV]');
%     ylabel('Waist size _{FWHM} [\mum] s(-), t(--), orig(*)');
figure(243);
clf
    subplot(1,2,1);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-','color','r','linewidth',1.5);
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-1)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','--','color','b','linewidth',1.5);
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-2)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-.','color','g','linewidth',1.5);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    ylim([0 300]);
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Energy [eV]');
    ylabel('Spot size [\mum]');
    legend('entrance','1m','2m');
    
    subplot(1,2,2);
    hold all
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','-','color','r','linewidth',1.5);
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-1)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','--','color','b','linewidth',1.5);
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-2)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','-.','color','g','linewidth',1.5);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    ylim([0 300]);
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Energy [eV]');
    ylabel('Spot size [\mum]');
    
    set(figure(19), 'Position', [300, 300, 450, 650]);  
    set(figure(20), 'Position', [310, 310, 450, 650]);  
    set(figure(241), 'Position', [320, 320, 450, 650]);  
    set(figure(242), 'Position', [330, 330, 450, 650]); 
    set(figure(243), 'Position', [340, 340, 1026, 290]); 

    
end
z3_t_offset-slitpos
f_m_t
%1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1)