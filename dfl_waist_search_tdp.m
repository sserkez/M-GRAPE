clear all

nm_p='500_tdp.out';
inread;
%break
energy=1239.8/xlamds/1e9;
[X,nslice]=fieldimport_all('500_tdp.out.dfl',M,1);

X=fftshift(fft(X,[],3),3);

P0=reshape(sum(sum(abs(X).^2)),[],1);

[H{1}]=fieldplot(1,X(:,:,500),leng,'initial field',1);

%%
%disp(P0);

% S_size=1.0e-5;
% S_pos=3.1;
% U_to_Gr=1.265;

%X=prop_TF(X,leng,xlamds,U_to_Gr);

%X=fieldgaussian(301,leng,S_size,S_size,-S_pos,-S_pos,xlamds,1000);

%[M,~]=size(X);
%dx=leng/M;
xscale=((M-1)/2+1-(1:M))*dx;

Z0=4:-0.3:-6;
%Z0=0;
index=size(Z0,2);


% old algorithm
% for i=1:index
%     
%     for s=1:nslice
% 
%     X1=prop_TF(X(:,:,s), leng, xlamds,Z0(i));
%     %[H{2}]=fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z0(i)),'m']);
%         try
%             fwhm(i,s)=(findFWHM(xscale,sum(abs(X1).^2,2)')+findFWHM(xscale,sum(abs(X1).^2,1)'))/2;
%         catch
%             fwhm(i,s)=0;
%         end
%     maxI(i,s)=max(max(abs(X1).^2));
% 
%     end
%     
% %     figure(2243)
% %     plot(Z0,fwhm(:,s));
% %     title(num2str(s));
%     
% end

fwhm_weighted=zeros(index,1);
maxI=zeros(index);

for i=1:index

    X1_a=zeros(M,M);
    for s=1:nslice
    X1=prop_TF(X(:,:,s), leng, xlamds,Z0(i));
    %[H{2}]=fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z0(i)),'m']);
    X1_a=X1_a+abs(X1).^2;
    end

    try
        fwhm_weighted(i)=(findFWHM(xscale,sum(X1_a,2)')+findFWHM(xscale,sum(X1_a,1)'))/2;
    catch
        fwhm_weighted(i)=0;
    end
    maxI(i)=max(max(abs(X1).^2));
%     figure(2243)
%     plot(Z0,fwhm(:,s));
%     title(num2str(s));

end
%%
 S_size=13e-6;
 S_pos=2.1;
 U_to_Gr=1.265;
 X=fieldgaussian(151,leng,S_size,S_size,-S_pos+U_to_Gr,-S_pos+U_to_Gr,xlamds,1000);
% U_to_Gr=1.265;

    for i=1:index
    X1=prop_TF(X, leng, xlamds,Z0(i));
    %[H{2}]=fieldplot(2,X1,leng,['propagated field, Z=',num2str(Z0(i)),'m']);
        try
            fwhm_g(i)=findFWHM(xscale,sum(abs(X1).^2,2)');
        catch
            fwhm_g(i)=0;
        end
    maxI_g(i)=max(max(abs(X1).^2));
    end



% fwhm_weighted=zeros(index,1);
% for i=1:nslice
%     fwhm_weighted=fwhm_weighted+fwhm(:,i).*P0(i)/sum(P0);
% end

figure(2244)
%hold all %to be removed
plot(Z0,sum(fwhm_weighted,s),'linewidth',2);
hold all
plot(Z0,fwhm_g,'linewidth',2);
xlabel('distance [m]');
ylabel('FWHM intensity [m]');
hold off

figure(2245)
plot(Z0,sum(maxI,2)/max(sum(maxI,2)),'linewidth',2);
hold all
plot(Z0,maxI_g/max(maxI_g),'linewidth',2);
xlabel('distance [m]');
ylabel('Central intensity [a.u.]');
hold off

figure(2249)
plot(Z0,fwhm_weighted,'linewidth',2);
xlabel('distance [m]');
ylabel('FWHM intensity [m]');