%calculating the required resolving power to obtain a transform-limited
%pulse for a given energy and ebeam duration

% clear all
% 
% energy=300:50:1000;
% xlamds=1239.8411./energy*1e-9;
% c=3e8;
% omega=2*pi*c./xlamds;
% 
% for t_fwhm=[10e-15 20e-15 30e-15];
% omega_fwhm=2.7726./t_fwhm;
% res_power=omega./omega_fwhm;
% 
% figure(40496)
% plot(energy,res_power,'linewidth',2)
% hold all
% end
% hold off
% xlabel('Photon energy [eV]');
% ylabel('Resolving power');
% legend('10 fs','20 fs','30 fs');

clear all
energy=300:50:1200;
xlamds=1239.8411./energy*1e-9;
c=3e8;
%omega=2*pi*c./xlamds;

for t_fwhm=[5e-15 10e-15 20e-15 40e-15];
    l_fwhm=t_fwhm.*c;
%omega_fwhm=2.7726./t_fwhm;
res_power=l_fwhm./xlamds./(2*log(2)./pi);

figure(40496)
plot(energy,res_power,'linewidth',2)
hold all
end
hold off
xlabel('Photon energy [eV]');
ylabel('Resolving power');
legend('5 fs (1.5\mum)','10 fs (3\mum)','20 fs (6\mum)','40 fs (12\mum)');
xlim([energy(1) energy(end)]);
ylim([0 2e4]);
clear all
%t_fwhm=0:1e-15:20e-15;
energy=300:50:1200;
c=3e8;
%omega=2*pi*c./xlamds;

for res_power=[3000 5000 7000 10000];
    xlamds=1239.8411./energy*1e-9;

%omega_fwhm=2.7726./t_fwhm;
    l_fwhm=res_power.*(xlamds.*2*log(2)./pi);
    t_fwhm=l_fwhm./c*1e15;
figure(40495)
plot(energy,t_fwhm,'linewidth',2)
hold all
end
hold off
xlabel('Photon energy [eV]');
ylabel('Pulse duration [fs]');
legend('3000','5000','7000','10000');
ylim([0 45]);
xlim([energy(1) energy(end)]);

set(figure(40496), 'Position', [100, 100, 500, 400]);
set(figure(40495), 'Position', [100, 100, 500, 400]);

%%
c=3e8;
energy =300;
res_power=9400;
xlamds=1239.8411./energy*1e-9;
l_fwhm=res_power.*(xlamds.*2*log(2)./pi);
t_fwhm=l_fwhm./c*1e15;
disp(l_fwhm*1e6);