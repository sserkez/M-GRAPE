%compare amplitudes and phases of the radiation obtained from .out files
clear all

nm_p{1}='D:\Work\!PROJECTS\Phase_controlled_harmonics\SASE3_v4\stage_1.out';
nm_p{2}='D:\Work\!PROJECTS\Phase_controlled_harmonics\SASE3_v4\stage_2.out';
 
nm_p{1}='D:\Work\!PROJECTS\ocelot_test\test_13\run_0\run.0.s1.gout';
nm_p{2}='D:\Work\!PROJECTS\ocelot_test\test_13\run_0\run.0.s3.gout';

 %% Import from .out
 
 for Di=1:2 %Data index
    if ~isempty(nm_p{Di})
        d(Di)=outread(nm_p{Di},1,0,2); 
    else
        d(Di)=struct;
    end
end
disp ('------------end-of-import-------------');

%%
Sn=min(d(1).outp.Sn,d(2).outp.Sn);


 for Di=1:2 %Data index
        d(Di).outp.Sscale=d(Di).outp.Sscale(1:Sn);
        d(Di).outp.power.v=d(Di).outp.power.v(1:Sn,:);
        d(Di).outp.p_mid.v=d(Di).outp.p_mid.v(1:Sn,:);
        d(Di).outp.phi_mid.v=d(Di).outp.phi_mid.v(1:Sn,:);
end



%% General plots

phasezoom_value=0.5; %if 0 -> no zoom

Z1=100; %[m]
Z2=12;


% H{1}=outpot_e(1,d(1));
% H{2}=outpot_ph(2,d(1));
% H{3}=outpot_e(11,d(2));
% H{4}=outpot_ph(12,d(2));
% 
% outplot_z(4,d(1),Z1);
% outplot_z(5,d(2),Z2);

% Compare phases/intensities

Sscale1=d(1).outp.Sscale;
Zscale1=d(1).outp.Zscale;
Sscale2=d(2).outp.Sscale;
Zscale2=d(2).outp.Zscale;
xlamds1=d(1).inp.xlamds*1e9;
xlamds2=d(2).inp.xlamds*1e9;


if Z1>max(Zscale1)
    disp('Z1 parameter exceeds data limits');
    Z1=max(Zscale1);
elseif Z1<min(Zscale1)
    Z1=min(Zscale1);
    disp('Z1 parameter exceeds data limits');
end
Zi1=find(Zscale1>=Z1,1,'first');
Z1=Zscale1(Zi1);
minsc1=min(Sscale1);
maxsc1=max(Sscale1);

if Z2>max(Zscale2)
    disp('Z2 parameter exceeds data limits');
    Z2=max(Zscale2);
elseif Z2<min(Zscale2)
    Z2=min(Zscale2);
    disp('Z2 parameter exceeds data limits');
end
Zi2=find(Zscale2>=Z2,1,'first');
Z2=Zscale2(Zi2);
minsc2=min(Sscale2);
maxsc2=max(Sscale2);

shift_points=0;
d(1).outp.power.v=circshift(d(1).outp.power.v,[shift_points 0]);
d(1).outp.phi_mid.v=circshift(d(1).outp.phi_mid.v,[shift_points 0]);
amplitude1=sqrt(d(1).outp.power.v(:,Zi1));
amplitude2=sqrt(d(2).outp.power.v(:,Zi2));



A1=max(amplitude1);
A2=max(amplitude2);

% figure(6)
% clf
% hline(1)=plot(Sscale1,power1);
% hold all
% hline(2)=plot(Sscale2,power2);
% set(hline(1),'LineWidth',2,'color','r','linestyle','-');
% set(hline(2),'LineWidth',2,'color','b','linestyle','-');
% legend(num2str(xlamds1),num2str(xlamds2))

figure(7)
clf;
ax1=subplot(2,2,1);
%hline1(1)=plot(Sscale1,amplitude1/A1,'LineWidth',1,'color','r','linestyle','-');
hline1(1)=plot(Sscale1,amplitude1,'LineWidth',1,'color','r','linestyle','-');
hold on
%hline1(2)=plot(Sscale2,amplitude2/A2,'LineWidth',1,'color','b','linestyle','-');
hline1(2)=plot(Sscale2,amplitude2,'LineWidth',1,'color','b','linestyle','-');
hold off
legend(['\lambda_1=',num2str(xlamds1),'nm'],['\lambda_2=',num2str(xlamds2),'nm \newlineZ_2=',num2str(Z2),'m'])
title(['radiation amplitude envelopes normalized by peak values ratio ',num2str(A1/min(A1,A2)),':',num2str(A2/min(A1,A2))]);

ax2=subplot(2,2,2);

if A2>=A1
plot(Sscale1,amplitude2./amplitude1,'LineWidth',1,'color','k','linestyle','-');
else
plot(Sscale1,amplitude1./amplitude2,'LineWidth',1,'color','k','linestyle','-');
end
title('amplitude envelopes ratio')


phase1=[0; diff(unwrap((xlamds1/xlamds2).*d(1).outp.phi_mid.v(:,Zi1)))];
phase2=[0; diff(unwrap(d(2).outp.phi_mid.v(:,Zi2)))];


% figure(8)
% clf
ax3=subplot(2,2,3);
hline2(1)=plot(Sscale1,phase1,'LineWidth',1,'color','r','linestyle','-');
hold on
hline2(2)=plot(Sscale2,phase2,'LineWidth',1,'color','b','linestyle','-');
line([minsc1,maxsc1],[0,0],'color','k','linestyle','--');
hold off
%legend(['\lambda_1=',num2str(xlamds1),'nm'],['\lambda_2=',num2str(xlamds2),'nm'])
title('radiation phase (1-st stage values are scaled by \lambda_1/\lambda_2)')
xlabel('[m]');
ylabel('[rad]');
if phasezoom_value==0
    ylim([-pi pi]);
else
    ylim([-phasezoom_value phasezoom_value]);
end

ax4=subplot(2,2,4);
hline3=plot(Sscale1,phase2-phase1,'LineWidth',1,'color','k','linestyle','-'); %would not work with different Sscales!!!!!!!!!
hold on
line([minsc1,maxsc1],[0,0],'color','k','linestyle','--');
hold off
title('radiation phase difference')
xlabel('[m]');
ylabel('[rad]');
if phasezoom_value==0
    ylim([-pi pi]);
else
    ylim([-phasezoom_value phasezoom_value]);
end
linkaxes([ax1 ax2 ax3 ax4],'x');
xlim([max([minsc1 minsc2]) min([maxsc1 maxsc2])])
%Handle=linkprop([hline1(1) hline1(2) hline2(1) hline2(2) hline3], 'XLim');
