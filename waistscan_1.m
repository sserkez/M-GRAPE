function [z_opt,fwhm_min,fwhm_0]=waistscan_1(figure_number,X,leng,xlamds,Z_array)

%cross-check x/y directions

[M,~,N]=size(X);
dx=leng/M;
xscale=((M-1)/2+1-(1:M))*dx;

index=size(Z_array,2);

sigma=1.5e-5;
z=0;
for i=1:index

X1=prop_TF(X, leng, xlamds,Z_array(i));
%fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z_array(i)),'m']);
fieldplot(315,X1,leng,['propagated field, Z=',num2str(Z_array(i)),'m'],1);
            
    Ix=sum(sum(abs(X1).^2,3),1)'./sum(sum(sum(abs(X1).^2)));
    Iy=sum(sum(abs(X1).^2,3),2)./sum(sum(sum(abs(X1).^2)));

    %trying a variance approach
    mean_x(i)=sum(xscale.*Ix');
    mean_y(i)=sum(xscale.*Iy');
    std_x(i)=sqrt(sum((xscale-mean_x(i)).^2.*Ix'))*1e6;
    std_y(i)=sqrt(sum((xscale-mean_y(i)).^2.*Iy'))*1e6;
%     rms_x(i)=sqrt(sum(xscale.^2.*Ix'));
%     rms_y(i)=sqrt(sum(xscale.^2.*Iy'));

        try
    fwhm_x(i)=findFWHM(xscale,Ix);
    fwhm_y(i)=findFWHM(xscale,Iy);
        catch
        fwhm_x(i)=0;
        fwhm_y(i)=0;
        end
    maxI(i)=max(max(max(abs(X1).^2)));

end

% [fwhm_min,n]=min(fwhm);
% z_opt=Z_array(n);
% fwhm_0=findFWHM(xscale,sum(abs(X).^2,dim)');


figure(figure_number)
axes1=subplot(3,1,1);
plot(Z_array,std_x,'linewidth',2,'color','b','linestyle','--');
hold on
plot(Z_array,std_y,'linewidth',2,'color','b','linestyle','-');
hold off
ylabel('\sigma _{x(--) y(-)} [m]');
grid on

waistx_pos=Z_array(std_x==min(std_x));
waistx_val=min(std_x);
waisty_pos=Z_array(std_y==min(std_y));
waisty_val=min(std_y);


text(0,1,sprintf(' \\sigma_x = %.1f \\mum at %.1fm\n \\sigma_y = %.1f \\mum at %.1fm', waistx_val, waistx_pos, waisty_val, waisty_pos),...
                    'HorizontalAlignment','left','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','k');

axes2=subplot(3,1,2);
plot(Z_array,fwhm_x*1e6,'linewidth',2,'color','r','linestyle','--');
hold on
plot(Z_array,fwhm_y*1e6,'linewidth',2,'color','r','linestyle','-');
hold off
ylabel('FWHM _{x(--) y(-)} [m]');
grid on

axes3=subplot(3,1,3);
plot(Z_array,maxI/max(maxI),'linewidth',2,'color','k');
ylabel('Max. itensity [arb. units]');
xlabel('distance [m]');
grid on

%Handle=linkprop([axes1 axes2 axes3], 'XLim');
xlim([Z_array(1) Z_array(end)]);