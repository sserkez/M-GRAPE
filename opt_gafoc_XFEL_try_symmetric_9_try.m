%optimization script

%nowI try to move grating upstream the undulator to get larger waist at the
%undulator 2 and greater resolution at the slit


showpictures=1;
showpublicationpictures=0;
ienergy=300:20:1200;
%ienergy=1200;
xlamdsout=1239.8411./ienergy.*1e-9;

opt_parms_5;

U2_res_out=linspace(1,1,size(ienergy,2));
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
    
    xlamds=1239.8411/energy*1e-9;  
%old and probably obsolete values that give strange result. Recalculated
%based on the origin files and late report plots
% s1=(36.5-energy*0.0118)*1e-6/2.35; %sourse sigma
% z1=(1.787+2.37*energy/1000)+U1_to_G; %sourse distance

% s1=(32.86-0.00863*energy)*1e-6; %sourse sigma
% z1=(0.09226+0.00101*energy)+U1_to_G; %sourse distance

% s1=(32.86-0.00863*energy)*1e-6/1.2*sourcesize_x; %sourse sigma
% z1=((energy-300)/700+1)+U1_to_G; %sourse distance
% 
% %last recalculated
% s1=(25.7-(energy-300)/700*8.7)*1e-6;
% z1=0.85+(energy-300)/700*1.95+U1_to_G;
% 
% %last recalculated poly
% s1=(33.1-0.0283*energy+1.223e-5*energy.^2)*1e-6*sourcesize_x;
% z1=0.00327*energy-1.4506e-6*energy.^2+U1_to_G;
    
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
    U2_res_out(index)=U2_res;
end

if showpictures
    
    figure(19);
    clf
    subplot(2,1,1);
    %title('Source position');
    hold all
    plot(ienergy,-z1_out+U1_to_G,'linewidth',2,'color','b');
    plot(ienergy,-z1_out+U1_to_G+pi*w1_out.^2./xlamdsout,'color','r');
    plot(ienergy,-z1_out+U1_to_G-pi*w1_out.^2./xlamdsout,'color','g');
    plot(ienergy,linspace(0,0,numel(ienergy)),'LineWidth',1,'linestyle','--','color','k')%!!!!!!!!!!!!!!!!!!!!!!!!
    hold off
    xlabel('Energy [eV]');
    ylabel('Waist position [m]');
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w1_out./2.*2.35);
    hold off
    xlabel('Energy [eV]');
    ylabel('Waist size _{FWHM} [\mum]');
    
    figure(18);
    clf
    plot(ienergy,U2_res_out,'linewidth',2,'color','b');
    xlabel('Energy [eV]');
    ylabel('effective resolving power in undulator');
    xlim([min(ienergy) max(ienergy)]);
    ylim([min(U2_res_out) max(U2_res_out)]);
    
    
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
    ylim([-2 10]);
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
    xlabel('Photon energy [eV]');
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
    xlabel('Photon energy [eV]');
    ylabel('Waist FWHM size [\mum]');
    
    
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
    %xlabel('Photon energy [eV]');
    ylabel('Waist position [m]');
    xlabel('Photon energy [eV]');
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
    xlabel('Photon energy [eV]');
    ylabel('Size [\mum]');
    legend('Waist size','Spot at slit position')
%     subplot(3,1,3);
%     hold all
%     plot(ienergy,1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1));
%     hold off
%     xlabel('Photon energy [eV]');
%     ylabel('Spot at slit size _{FWHM} [\mum]');


    %w_slit=w2_t_out.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1);
%     figure(21);
%     plot(ienergy,1./resolution_slit_out);
%     hold off
%     xlabel('Photon energy [eV]');
%     ylabel('Resolution on slit');
% 
%     figure(22);
%     clf
%     hold all
%     plot(ienergy,N_gr_out);
%     hold off
%     xlabel('Photon energy [eV]');
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
    ylim([-0.5 12]);
    ylabel('Waist position [m]');
    xlabel('Photon energy [eV]');
    
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35,'linestyle','-','color','r','linewidth',2);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([20 70]);
    ylabel('Waist FWHM size [\mum]');
    xlabel('Photon energy [eV]');
    
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
    ylim([-1 12]);
    ylabel('Waist position [m]');
    xlabel('Photon energy [eV]');
    subplot(2,1,2);
    hold all
    plot(ienergy,1e6.*w4_t_out./2.*2.35,'linestyle','-','color','r','linewidth',2);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    ylim([20 70]);
    ylabel('Waist FWHM size [\mum]');
    xlabel('Photon energy [eV]');
    
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
%     xlabel('Photon energy [eV]');
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
%     xlabel('Photon energy [eV]');
%     ylabel('Waist size _{FWHM} [\mum] s(-), t(--), orig(*)');
figure(243);
clf
    subplot(1,2,1);
    hold all
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-','color','r','linewidth',1.5);
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-2)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','--','color','b','linewidth',1.5);
    plot(ienergy,1e6.*w4_s_out./2.*2.35.*sqrt(1+((z4_s_out-5)./(pi.*w4_s_out.^2./xlamdsout)).^2),'linestyle','-.','color','g','linewidth',1.5);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    ylim([0 150]);
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Photon energy [eV]');
    ylabel('Spot size [\mum]');
    
    
    subplot(1,2,2);
    hold all
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','-','color','r','linewidth',1.5);
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-1)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','--','color','b','linewidth',1.5);
    plot(ienergy,1e6.*w4_t_out./2.*2.35.*sqrt(1+((z4_t_out-2)./(pi.*w4_t_out.^2./xlamdsout)).^2),'linestyle','-.','color','g','linewidth',1.5);
    plot(ienergy,1e6.*w1_out./2.*2.35,'marker','*','color','k');
    hold off
    ylim([0 150]);
    if numel(ienergy)>1
        xlim([min(ienergy) max(ienergy)]);
    end
    xlabel('Photon energy [eV]');
    ylabel('Spot size [\mum]');
    legend('entrance','2m','5m');
    
    set(figure(19), 'Position', [300, 300, 450, 650]);  
    set(figure(20), 'Position', [310, 310, 450, 650]);  
    set(figure(241), 'Position', [320, 320, 450, 650]);  
    set(figure(242), 'Position', [330, 330, 450, 650]); 
    set(figure(243), 'Position', [340, 340, 1026, 290]); 

    
end
% z3_t_offset-slitpos
% f_m_t

disp(['Source size=',num2str(s1*1e6),'[um]']);
disp(['Source ang=',num2str(xlamds/4/pi/s1),'[rad]']);
disp(['Source distance=',num2str(z1),'[m]']);

%1e6.*w2_t_out./2.*2.35.*sqrt((z2_t_out-slitpos).^2./zr_slit_out.^2+1)