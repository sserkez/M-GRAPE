function [z_opt,fwhm_min,fwhm_0]=waistscan(X,leng,xlamds,Z_array,dim)

[M,~]=size(X);
dx=leng/M;
xscale=((M-1)/2+1-(1:M))*dx;



index=size(Z_array,2);

parfor i=1:index

X1=prop_TF(X, leng, xlamds,Z_array(i));
fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z_array(i)),'m']);
    try
        fwhm(i)=findFWHM(xscale,sum(abs(X1).^2,dim)');
        %fwhm(i)=mygaussfit(xscale*1e6,sum(abs(X1).^2,dim)')/1e6*2.355;
    catch
        fwhm(i)=0;
    end
maxI(i)=max(max(abs(X1).^2));

end

[fwhm_min,n]=min(fwhm);
z_opt=Z_array(n);
fwhm_0=findFWHM(xscale,sum(abs(X).^2,dim)');


figure(68)
subplot(2,1,1)
plot(Z_array,fwhm,'linewidth',2,'color','k');
ylabel('FWHM [m]');

subplot(2,1,2)
plot(Z_array,maxI/max(maxI),'linewidth',2,'color','k');
ylabel('Maximum Itenssity [a.u.]');
xlabel('distance [m]');
