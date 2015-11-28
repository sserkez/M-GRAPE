function [z_opt,fwhm_min,fwhm_0]=waistscan_g(figure_number,X,leng,xlamds,Z_array)

[M,~,N]=size(X);
dx=leng/M;
xscale=((M-1)/2+1-(1:M))*dx;

index=size(Z_array,2);



%X=fieldgaussian(M,leng,sigma,sigma,z,z,xlamds,1);
    for i=1:index

    X1=prop_TF(X, leng, xlamds,Z_array(i));
    %fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z_array(i)),'m']);
    fieldplot(315,X1,leng,['propagated field, Z=',num2str(Z_array(i)),'m'],1);

    
Xx=X1(round((M+1)/2),:,:)';
Xy=X1(:,round((M+1)/2),:);

Xx=abs(Xx).^2;
Xy=abs(Xy).^2;

Xx=Xx./sum(sum(Xx));
Xy=Xy./sum(sum(Xy));

 
    
 
    
        Ix=sum(sum(abs(X1).^2,3),1)'./sum(sum(sum(abs(X1).^2)));
        Iy=sum(sum(abs(X1).^2,3),2)./sum(sum(sum(abs(X1).^2)));
        
        
    Ix=Xx;
    Iy=Xy;
        
        mean_x(i)=sum(xscale.*Ix');
        mean_y(i)=sum(xscale.*Ix');
        std_x(i)=sqrt(sum((xscale-mean_x(i)).^2.*Ix'))*1e6;
        std_y(i)=sqrt(sum((xscale-mean_y(i)).^2.*Iy'))*1e6;
        maxI(i)=max(max(max(abs(X1).^2)));
    end

figure(figure_number)
subplot(3,1,1)
hold all
%delete(hline);
hline=plot(Z_array,std_x,'linewidth',2,'color','k','linestyle',':');
hold off

waistx_pos=Z_array(std_x==min(std_x));
if numel(waistx_pos)==2
    waistx_pos=mean(waistx_pos); %case if waist in between grid points
end
waistx_val=min(std_x);

text(1,1,sprintf(' \\sigma_{gauss} = %.1f \\mum at %.1fm\n', waistx_val, waistx_pos),...
                    'HorizontalAlignment','right','VerticalAlignment',...
                    'top','FontSize',10,'units','normalized','color','k');
                
subplot(3,1,2)
hold on
%delete(hline);
hline=plot(Z_array,std_x*2.35,'linewidth',2,'color','k','linestyle',':');
hold off


subplot(3,1,3);
hold on
plot(Z_array,maxI/max(maxI),'linewidth',2,'color','k','linestyle',':');
grid on