% reimaging optimization script
clear all
P0=1000;
M=151;
leng=10e-4;
xlamds=2.48E-09;


% Ns=30;
% Nz=30;
Ns=20;
Nz=20;
s=linspace(0.5e-5,4e-5,Ns); %sigma intensity
z=linspace(-2,8,Nz);

P=zeros(Ns,Nz);
iteration=0;
P_old=0;
for is=1:numel(s)
    for iz=1:numel(z)
        iteration=iteration+1;
        t_iter=tic;
%         while 1
            X=fieldgaussian(M,leng,s(is),s(is),z(iz),z(iz),xlamds,P0);
            [H{12}]=fieldplot(12,X,leng,'genegated field',1);
            fieldexport(X,'prop_field.dfl');
            %fieldexport(X,'dummy.dfl');
            fclose all;
            %pause(0.1);
            dos('genesis301');
            Xn=fieldimport('prop.out.dfl',M,1);
            [H{13}]=fieldplot(13,Xn,leng,'amplified field',1);
            P(is,iz)=sum(sum(abs(Xn).^2));
            
%             if P(is,iz)~=P_old
%                 P_old=P(is,iz);
%                 break
%             else
%                 delete('prop_field.dfl');
%             end
%          end
        
        if Ns>1 && Nz>1
            figure(5986)
            imagesc(z,s*1e6,P./max(max(P)));
            colorbar;
            xlabel('Waist position downstream undulator [m]');
            ylabel('Waist size \sigma_{intensity}, [\mum]');
        end

        if Ns==1 && Nz~=1
            figure(5987)
            plot(z,P)
        end

        if Nz==1 && Ns~=1
            figure(5986)
            plot(s,P)
        end
        
        time=toc(t_iter);
        disp(['time per iteration  ',num2str(time),' sec.']);
        
        if time*(Ns*Nz-iteration)/3600>1
            disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)/3600),' hours.']);
        elseif time*(Ns*Nz-iteration)/60>1
            disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)/60),' min.']);
        else
            disp(['estimated finish    ',num2str(time*(Ns*Nz-iteration)),' sec.']);
        end
        
    end
end