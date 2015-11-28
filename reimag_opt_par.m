% reimaging optimization script
M=151;
leng=20e-4;
xlamds=1.77e-9;

Ns=10;
Nz=10;
s=linspace(0.4e-5,3e-5,Ns);
z=linspace(-2,6,Nz);

P=cell(Ns,1);
iteration=0;
P_old=0;

parfor is=1:numel(s)
    t_iter=tic;
    %iteration=iteration+1;
    for iz=1:numel(z)
        
        

            X=fieldgaussian(M,leng,s(is),s(is),z(iz),z(iz),xlamds,P0);
    %        [H{12}]=fieldplot(12,X,leng,'amplified field');
            fieldexport(X,'prop_field.dfl');
            fieldexport(X,'dummy.dfl');
            dos('genesis301');
            Xn=fieldimport('700_prop.out.dfl',M,1);
    %        [H{13}]=fieldplot(13,Xn,leng_u2,'amplified field');
            P{is}(iz)=sum(sum(abs(Xn).^2));
            
%             if P(is,iz)~=P_old
%                 P_old=P(is,iz);
%                 break
%             end

        
        
        
        
        
        
        
    
    end

%     time=toc(t_iter);
%         disp(['time per iteration  ',num2str(time),' sec.']);
%         
%     if time*(Ns*Nz-iteration)/3600>1
%             disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)/3600),' hours.']);
%         elseif time*(Ns*Nz-iteration)/60>1
%             disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)/60),' min.']);
%         else
%             disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)),' sec.']);
%     end
    
end
    
P=cell2mat(P);

if Ns>1 && Nz>1
    figure(5986)
    imagesc(z,s,P);
    colorbar;
end

if Ns==1 && Nz~=1
    figure(5987)
    plot(z,P)
end

if Nz==1 && Ns~=1
    figure(5986)
    plot(s,P)
end