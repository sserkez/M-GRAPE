% interpolate the field.
% X - old field file
% if relative=1:
%   interpL - ratio of new_grid_size/old_grid_size
%   interpM - ratio of new_mesh_points/old_mesh_points
% else
%   interpL and interpM are new grid and mesh size
% method - argument can be any of the following strings that specify alternative interpolation methods: 'linear', 'nearest', 'cubic', or 'spline'
function X_interp=fieldinterpolate_a(X_orig,Lx_orig,Ly_orig,relative,interpMx,interpMy,interpLx,interpLy,method)

[My_orig,Mx_orig,~]=size(X_orig);

if relative

        if interpMx==0 || interpLx==0 ||interpMy==0 || interpLy==0
            error('Interpolation values should be >0');
        elseif interpMx==1 && interpLx==1 && interpMy==1 && interpLy==1
            X_interp=X_orig;
        else
            
            if mod(round(interpMx*Mx_orig),2) == 0
                Mx_interp=round(interpMx*Mx_orig)+1;
            else
                Mx_interp=round(interpMx*Mx_orig);
            end
            
            if mod(round(interpMy*My_orig),2) == 0
                My_interp=round(interpMy*My_orig)+1;
            else
                My_interp=round(interpMy*My_orig);
            end
            
            
            if interpMx==interpLx && interpMx>1 %&& mod(round(interpM),2) ~= 0
                nmatrix=zeros(My_orig,Mx_interp);                
                nmatrix(:,round((Mx_interp-Mx_orig+1)/2:(Mx_interp-Mx_orig)/2+Mx_orig))=X_orig;
                X_interp1=nmatrix;
            elseif interpMx==interpLx && interpMx<1 %&& mod(round(1/interpM),2) ~= 0
                X_interp1=X_orig(:,round((Mx_orig-Mx_interp)/2+1:(Mx_orig-Mx_interp)/2+Mx_interp));
            elseif size(X_orig,2)~=1
                Lx_interp=Lx_orig*interpLx;
                Ly_interp=Ly_orig;
                %Ly_interp=Ly_orig*interpLy;
                

                
                if interpLx<1
                    %Po=sum(sum(abs(X_orig(round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)),round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):    round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)))).^2));
                    Po=sum(sum(abs(X_orig(:,round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig)+1:    round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)))).^2));
                else
                    Po=sum(sum(abs(X_orig).^2));
                end
                
                [xxo,yyo]=meshgrid((Mx_orig+1)/2+1-(1:Mx_orig),(My_orig+1)/2+1-(1:My_orig));
                dxo=Lx_orig/Mx_orig;
                dyo=Ly_orig/My_orig;
%                 size(xxo)
%                 xxo(1,1)*dxo
%                 xxo(end,end)*dxo
%                 yyo(1,1)*dyo
%                 yyo(end,end)*dyo
                [xxi,yyi]=meshgrid((Mx_interp+1)/2+1-(1:Mx_interp),(My_orig+1)/2+1-(1:My_orig));
                dxi=Lx_interp/Mx_interp;
                dyi=Ly_orig/My_orig;
%                 size(xxi)
%                 xxi(1,1)*dxi
%                 xxi(end,end)*dxi
%                 yyi(1,1)*dyi
%                 yyi(end,end)*dyi
%                 size(X_orig)
                if size(X_orig,1)==1
                    X_interp1=interp1(xxo(1,:)*dxo,X_orig(1,:),xxi(1,:)*dxi,method,min(min(X_orig)));
                else
                    X_interp1=interp2(xxo*dxo,yyo*dyo,X_orig,xxi*dxi,yyi*dyi,method,min(min(X_orig)));
                end
%                 size(X_interp1);
                Pi=sum(sum(abs(X_interp1).^2));
                
                if Pi~=0
                X_interp1=X_interp1*sqrt(Po/Pi);
                end
            else
                X_interp1=X_orig;
            end
      %      disp('done x')
%             Po/Pi
%             size(X_interp1)
%             My_interp=My_interp
%             My_orig=My_orig
%             X_interp=X_interp1;
            if interpMy==interpLy && interpMy>1 %&& mod(round(interpM),2) ~= 0
                nmatrix=zeros(My_interp,Mx_interp);                
                nmatrix(round((My_interp-My_orig+1)/2:(My_interp-My_orig)/2+My_orig),:)=X_interp1;
                X_interp=nmatrix;
            elseif interpMy==interpLy && interpMy<1 %&& mod(round(1/interpM),2) ~= 0
                X_interp=X_interp1(round((My_orig-My_interp)/2+1:(My_orig-My_interp)/2+My_interp),:);
            elseif size(X_interp1,1)~=1
                Lx_interp=Lx_orig*interpLx;
                Ly_interp=Ly_orig*interpLy;
                
                %Ly_interp=Ly_orig*interpLy;
                
                if interpLy<1
                    %Po=sum(sum(abs(X_orig(round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)),round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):    round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)))).^2));
                    Po=sum(sum(abs(X_interp1(round(My_orig*(Ly_orig-Ly_interp)/2/Ly_orig)+1:round(My_orig*(1-(Ly_orig-Ly_interp)/2/Ly_orig)),:)).^2));
                else
                    Po=sum(sum(abs(X_interp1).^2));
                end
                
                [xxo,yyo]=meshgrid((Mx_interp+1)/2+1-(1:Mx_interp),(My_orig+1)/2+1-(1:My_orig));
                dxo=Lx_interp/Mx_interp;
                dyo=Ly_orig/My_orig;
%                 yyo(1,1)*dyo
%                 yyo(end,end)*dyo
                [xxi,yyi]=meshgrid((Mx_interp+1)/2+1-(1:Mx_interp),(My_interp+1)/2+1-(1:My_interp));
                dxi=Lx_interp/Mx_interp;
                dyi=Ly_interp/My_interp;
%                 yyi(1,1)*dyi
%                 yyi(end,end)*dyi

%                 if size(X_orig,1)==1
%                     X_interp1=interp1(xxo(1,:)*dxo,X_orig(1,:),xxi(1,:)*dxi,method,min(min(X_orig)));
%                 else
%                     X_interp1=interp2(xxo*dxo,yyo*dyo,X_orig,xxi*dxi,yyi*dyi,method,min(min(X_orig)));
%                 end
                X_interp=interp2(xxo*dxo,yyo*dyo,X_interp1,xxi*dxi,yyi*dyi,method,min(min(X_orig)));
                Pi=sum(sum(abs(X_interp).^2));
                
                if Pi~=0
                X_interp=X_interp*sqrt(Po/Pi);
                end
                
            else
                X_interp=X_interp1;
        %    disp('done y')
            end
%             Po/Pi
       end
else
    
    Lx_interp=interpLx;
    Mx_interp=round(interpMx);
    Ly_interp=interpLy;
    My_interp=round(interpMy);
    
    if Lx_interp<Lx_orig && Ly_interp<Ly_orig
        Po=sum(sum(abs(X_orig(round(My_orig*(Ly_orig-Ly_interp)/2/Ly_orig):    round(My_orig*(1-(Ly_orig-Ly_interp)/2/Ly_orig)),round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)))).^2));
    elseif Lx_interp<Lx_orig && Ly_interp>=Ly_orig
        Po=sum(sum(abs(X_orig(:,round(Mx_orig*(Lx_orig-Lx_interp)/2/Lx_orig):round(Mx_orig*(1-(Lx_orig-Lx_interp)/2/Lx_orig)))).^2));
    elseif Ly_interp<Ly_orig && Lx_interp>=Lx_orig
        Po=sum(sum(abs(X_orig(round(My_orig*(Ly_orig-Ly_interp)/2/Ly_orig):    round(My_orig*(1-(Ly_orig-Ly_interp)/2/Ly_orig)),:)).^2));
    else
        Po=sum(sum(abs(X_orig).^2));
    end
    
    [xxo,yyo]=meshgrid((Mx_orig-1)/2+1-(1:Mx_orig),(My_orig-1)/2+1-(1:My_orig));
%     xxo(1,1)
%     xxo(end,end)
%     yyo(1,1)
%     yyo(end,end)
    dxo=Lx_orig/Mx_orig;
    dyo=Ly_orig/My_orig;
    
    [xxi,yyi]=meshgrid((Mx_interp-1)/2+1-(1:Mx_interp),(My_interp-1)/2+1-(1:My_interp));
%     xxi(1,1)
%     xxi(end,end)
%     yyi(1,1)
%     yyi(end,end)
    dxi=Lx_interp/Mx_interp;
    dyi=Ly_interp/My_interp;
    
    X_interp=interp2(xxo*dxo,yyo*dyo,X_orig,xxi*dxi,yyi*dyi,method,min(min(X_orig)));
    Pi=sum(sum(abs(X_interp).^2));
    
    if Pi~=0
    X_interp=X_interp*sqrt(Po/Pi);
    end
    
end