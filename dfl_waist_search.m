clear all

nm_p='1000.out';
inread;
%break
energy=1239.8/xlamds/1e9;
X=fieldimport('1000.out.dfl',M,1);

[H{1}]=fieldplot(1,X,leng,'initial field');
%[h2]=fieldplot(2,X1.*2,leng,leng,'field');
P0=sum(sum(abs(X).^2));
%disp(P0);
S_size=1.0e-5;
S_pos=3.1;
U_to_Gr=1.265;
X=prop_TF(X,leng,xlamds,U_to_Gr);
%%%X=fieldgaussian(301,leng,S_size/2.35,S_size/2.35,-(S_pos-U_to_Gr),-(S_pos-U_to_Gr),xlamds);
%X=fieldgaussian(301,leng,S_size,S_size,-S_pos+1,-S_pos+1,xlamds,1000);
X=fieldgaussian(301,leng,S_size,S_size,-S_pos,-S_pos,xlamds,1000);

%[H{1}]=fieldplot(1,X,leng,'initial field');

[M,~]=size(X);
dx=leng/M;
xscale=((M-1)/2+1-(1:M))*dx;

Z0=6:-0.2:-8;
%Z0=0;
index=size(Z0,2);

for i=1:index

X1=prop_TF(X, leng, xlamds,Z0(i));
%[H{2}]=fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z0(i)),'m']);
    try
        fwhm(i)=findFWHM(xscale,sum(abs(X1).^2,2)');
    catch
        fwhm(i)=0;
    end
maxI(i)=max(max(abs(X1).^2));

end

figure(69)
subplot(2,1,1)
hold all
plot(Z0,fwhm, 'linewidth',2);
ylabel('FWHM [m]');
hold off
subplot(2,1,2)
hold all
plot(Z0,maxI/max(maxI), 'linewidth',2);
ylabel('Maximum Itenssity [a.u.]');
xlabel('distance [m]');
hold off


% %clear all
% nm='mod.out.dfl';
% 
% fd=fopen(nm,'r');
% 
% datarun=1;
% %1-imported
% %2-template
% 
%  showpictures=1;    
%  trldZ=2;
% 
%  vertical=1;
% 
%  if datarun==1
%      
%      inread;
% 
%      
% %      M=301;
% %      rxbeam=  1.266636E-05;
% %      rybeam=  1.564884E-05;
% %      xlamds=  4.141000E-09;
% %      zrayl =  3.015540E-01;
% %      rmax0 =  20.000000E+00;
% 
%      rbeam=sqrt(rxbeam^2+rybeam^2);
%      ray=sqrt(zrayl*xlamds/pi);
%      if dgrid==0
%         leng=rmax0*(rbeam+ray); %%%! *2
%      else
%          leng=dgrid*2;
%      end
%      dx=leng/M;     dy=dx;
%      K=2*pi/xlamds;
%      
%      if fd==-1
%          error('no such file');
%      end
% 
%      x=fread(fd,2*M*M,'double');
%      b=x(1:2:end)+1i*x(2:2:end);
%      x1=reshape(b,M,M);
%      clear x b
% 
%  elseif datarun==2
% 
% %    [xx,yy]=meshgrid(exp(-((((M-1)/2+1-(1:M))/5).^2)));
% %    x1=(xx.*yy+0*rand(M));
% 
% % % % % % %gauss
% % % % % % M=301;
% % % % % % leng=5e-4;
% % % % % % dx=leng/M; dy=dx;
% % % % % % 
% % % % % % inpsigm=2e-5;
% % % % % % xlamds=4e-9;
% % % % % % K=2*pi/xlamds;
% 
% M=301;
%      rxbeam=  1.266636E-05;
%      rybeam=  1.564884E-05;
%      xlamds=  4.141000E-09;
%      zrayl =  3.015540E-01;
%      rmax0 =  20.000000E+00;
%      rbeam=sqrt(rxbeam^2+rybeam^2);
%      ray=sqrt(zrayl*xlamds/pi);
%      leng=rmax0*(rbeam+ray); %%%! *2
%      dx=leng/M;     dy=dx;
%      K=2*pi/xlamds;
% 
% 
% q=1j*pi*(2*inpsigm^2)/xlamds-0;
% [xx,yy]=meshgrid((M-1)/2+1-(1:M));
% x1=exp(-1j.*K.*((xx.*dx).^2+(yy*dx).^2)./2./q);
% 
% % %sinc
% % [xx,yy]=meshgrid((M-1)/2+1-(1:M));
% % x1=20*(abs((yy).*dx)<inpsigm).*(abs((xx).*dx)<inpsigm);%.*(abs(xx.*dx)<inpsigm);
%    
% 
% 
%  else break
%  end
%  
%  
%  
%  linefitw=1+100:M-100;
%  
% I1=abs(x1).^2;
% 
% xscale=((M-1)/2+1-(1:M))*dx; yscale=xscale;
% x=dx.*(1:M);
% y=x;
% 
% inppower=sum(sum(x1.*conj(x1)));
% 
% 
% 
% %if showpictures
% fig1=figure(1);
% clf;
% set(fig1,'name','Intensity_before','numbertitle','off')
% %imagesc(x1.*conj(x1));
% imagesc(xscale,yscale,I1);
% axis equal tight
% colorbar('SouthOutside');
% title(['Intensity, z= 0m']);
% V=axis;
% text(double(V(1)),double(V(3)),sprintf('  \\lambda = %.3e\n  P_{sum} = %.3e', xlamds, inppower),...
%        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10,'color','white');
% 
% fig2=figure(2);
% clf;
% %surface(imag(x1),'EdgeColor','none');
% imagesc(xscale,yscale,angle(x1));
% axis equal tight
% set(fig2,'name','Phase_before','numbertitle','off');
% colorbar('SouthOutside');
% title(['Phase, z= 0m']);
% %end
% 
% 
% 
% Zr=K*leng^2/2;
% 
% Z0=-5:0.5:1;
% %Z0=0;
% index=size(Z0,2);
% sigmaout=linspace(1,1,index);
% Ifwhmout=linspace(1,1,index);
% I2max=linspace(1,1,index);
% 
% % FN=leng^2/xlamds/Z0(1);
% 
% %x2=f_prop(0.005, 500, x1, xlamds);
% %x2=prop_IR(x1, 0.0005, xlamds, 2000);
% %x2=f_2D_prop_spectr(dx,dx,0.000001,x1);
% %x2=f_2D_prop_fresnel(dx,dx,0,x1);
% 
% 
% 
% 
% 
% 
% for index=1:size(Z0,2); %Z=Z0
%     Z=Z0(index);
% 
% x2=prop_TF(x1, leng, xlamds, Z);
% I2=abs(x2);
% 
% if vertical
%     [sigma,mu,A]=mygaussfit(xscale(linefitw), I2((M-1)/2+1,linefitw).*(I2((M-1)/2+1,linefitw)>max(I2((M-1)/2+1),linefitw)/trldZ)); %Vertical
% else
%     [sigma,mu,A]=mygaussfit(xscale(linefitw), I2(linefitw,(M-1)/2+1).*(I2(linefitw,(M-1)/2+1)>max(I2(linefitw,(M-1)/2+1))/trldZ)); %horizontal
% end
% 
% %Ifwhm=findFWHM(xscale(linefitw), (I2(linefitw,(M-1)/2+1).*(I2(linefitw,(M-1)/2+1)>max(I2(linefitw,(M-1)/2+1))/trldZ))');
% if showpictures
% fig4=figure(4);
% clf;
% set(fig4,'name','Intensity_after','numbertitle','off');
% %imagesc((x2).*conj((x2)));
% imagesc(xscale,yscale,I2);
% axis equal tight;
% title(['Intensity, z= ',num2str(Z),' m']);
% colorbar('SouthOutside');
% 
% fig5=figure(5);
% clf;
% set(fig5,'name','Phase_after','numbertitle','off');
% imagesc(xscale,yscale,angle(x2));
% axis equal tight;
% title(['Phase z= ',num2str(Z),' m']);
% colorbar('SouthOutside');
% 
% fig6=figure(6);
% clf;
% set(fig6,'name','I_fit','numbertitle','off');
% if vertical
%     plot(xscale, I2((M-1)/2+1,:),'LineWidth',2,'linestyle','*');
% else
%     plot(xscale, I2(:,(M-1)/2+1),'LineWidth',2,'linestyle','*');
% end
% hold all
% %plot(xscale, I21(:,(M-1)/2+1),'LineWidth',2,'linestyle','*');
% plot(xscale,A.*exp(-(xscale-mu).^2/(2*sigma^2)),'LineWidth',2,'color','r');
% hold off
% axis tight
% 
% % fig8=figure(8);
% % clf;
% % set(fig5,'name','Phasefront within 2\\sigma','numbertitle','off');
% % slls=(M-1)/2+1-round(sigma/dx):(M-1)/2+1+round(sigma/dx); %two sigma line length
% % plot(slls,angle(x2(slls,(M-1)/2+1)));
% % axis tight
% % V=axis;
% % text(double(V(1)),double(V(4)),sprintf('  std = %.3f \n', std(angle(x2(slls,(M-1)/2+1)))),...
% %        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10);
% 
% end
% 
% 
% 
% % % % sigmaxy=coeffvalues(gspotfit(xscale,yscale,I2));
% % % % sigmaxy=(sigmaxy(2)+sigmaxy(3))/2;
% 
% % pause(1);
% %waitforbuttonpress;
% %sigmaout=[sigmaout sigma];
% sigmaout(index)=sigma;
% % Ifwhmout(index)=Ifwhm;
% I2max(index)=max(max(I2));
% 
% end
% 
% if size(Z0,2)==1
%     Z=Z0;
% end
% 
%     Z=Z0(I2max==(max(I2max)));
%     %Z=Z0(sigmaout==(min(sigmaout)));
%     %Z=Z0(Ifwhmout==(min(Ifwhmout)));
%     x2=prop_TF(x1, leng, xlamds, Z);
%     I2=abs(x2);
% %     x21=prop_IR(x1, leng, xlamds, Z);
% %     I21=abs(x21);
%     %
%     if vertical
%         [pos,~] =find(I2==max(I2(:))); %vertical
%         [sigma,mu,A]=mygaussfit(xscale(linefitw), I2(pos,linefitw).*(I2(pos,linefitw)>max(I2(pos,linefitw),1)/trldZ));
%     else 
%         [~,pos] =find(I2==max(I2(:))); %horizontal
%         [sigma,mu,A]=mygaussfit(xscale(linefitw), I2(linefitw,pos).*(I2(linefitw,pos)>max(I2(linefitw,pos),1)/trldZ));
%     end
% %     [sigma,mu,A]=mygaussfit(xscale(linefitw),
% %     I2(linefitw,pos).*(I2(linefitw,pos)>max(I2(linefitw,pos),1)/trldZ));
% 
% 
% fig4=figure(4);
% clf;
% set(fig4,'name','Intensity_after','numbertitle','off');
% %imagesc((x2).*conj((x2)));
% imagesc(xscale,yscale,abs(x2));
% axis equal tight;
% title(['Intensity, z= ',num2str(Z),' m']);
% colorbar('SouthOutside');
% 
% 
% fig5=figure(5);
% clf;
% set(fig5,'name','Phase_after','numbertitle','off');
% imagesc(xscale,yscale,angle(x2));
% axis equal tight;
% title(['Phase, z= ',num2str(Z),' m']);
% colorbar('SouthOutside');
% 
% fig6=figure(6);
% set(fig6,'name','I_fit','numbertitle','off');
% %plot(xscale, I2(:,pos),'LineWidth',2,'linestyle','*'); %vertical
% if vertical
%     plot(xscale, I2(pos,:),'LineWidth',2,'linestyle','*'); %vertical 
% else
%     plot(xscale, I2(:,pos),'LineWidth',2,'linestyle','*'); %vertical 
% end
% hold all
% % plot(xscale, I21(:,(M-1)/2+1),'LineWidth',2,'linestyle','*');
% plot(xscale,A.*exp(-(xscale-mu).^2/(2*sigma^2)),'LineWidth',2,'color','r');
% hold off
% axis tight
% title(['Line z= ',num2str(Z),' m']);
% %end
% 
% fig7=figure(7);
% set(fig7,'name','Sigma & I vs length','numbertitle','off');
% %set(fig7,'name','I vs length','numbertitle','off');
% subplot(2,1,1);
% plot(Z0,sigmaout);
% %plot(Z0,I2max);
% %hold all
% %legend('sigma','FWHM');
% %hold off
% axis tight
% V=axis;
%    text(double(V(1)),double(V(4)),sprintf('  Z = %.3f \n  \\sigma=%.3e', Z, sigma),...
%        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10);
% subplot(2,1,2);
% plot(Z0,I2max);
% axis tight
%    
% % fig8=figure(8);
% % clf;
% % set(fig5,'name','Phasefront within 2\\sigma','numbertitle','off');
% % slls=(M-1)/2+1-round(sigma/dx):(M-1)/2+1+round(sigma/dx); %two sigma line length
% % plot(slls,angle(x2(slls,(M-1)/2+1)));
% % V=axis;
% % text(double(V(1)),double(V(4)),sprintf('  std = %.3f \n', std(angle(x2(slls,(M-1)/2+1)))),...
% %        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10);
% % axis tight
% 
% %text(Z0(sigmaout==(min(sigmaout))))
% 
% % if datarun==2
% % (sigmaout(end)/sigmaout(Z0==0))^2-(2*Z0(end)/K/(sigmaout(Z0==0))^2/2)^2;
% % end