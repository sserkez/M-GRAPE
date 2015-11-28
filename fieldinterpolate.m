% interpolate the field.
% X - old field file
% if relative=1:
%   interpL - ratio of new_grid_size/old_grid_size
%   interpM - ratio of new_mesh_points/old_mesh_points
% else
%   interpL and interpM are new grid and mesh size
% method - argument can be any of the following strings that specify alternative interpolation methods: 'linear', 'nearest', 'cubic', or 'spline'
function X_interp=fieldinterpolate(X_orig,L_orig,relative,interpM,interpL,method)

[M_orig,~]=size(X_orig);

if relative

        if interpM==0 || interpL==0
            error('Interpolation values should be >0');
        elseif interpM==1 && interpL==1
            X_interp=X_orig;
        elseif interpM==interpL && interpM>1 %&& mod(round(interpM),2) ~= 0


            if mod(round(interpM*M_orig),2) == 0
                M_interp=round(interpM*M_orig)+1;
            else
                M_interp=round(interpM*M_orig);
            end

            %X_interp=padarray(X_orig,[M_orig*(interpM-1),M_orig*(interpM-1)],0,'both');
        %     disp(M_orig);
        %     disp(interpM);
            nmatrix=zeros(M_interp);
        %     size(nmatrix)
        %     round(M_orig*(interpM-1)/2)
        %     round(M_orig*(interpM-1)/2+M_orig)

        %     size(nmatrix(round(M_orig*(interpM-1)/2:M_orig*(interpM-1)/2+M_orig),round(M_orig*(interpM-1)/2:M_orig*(interpM-1)/2+M_orig)))
        %     size(nmatrix(round((M_interp-M_orig)/2:(M_interp-M_orig)/2+M_orig),round((M_interp-M_orig)/2:(M_interp-M_orig)/2+M_orig)))
            nmatrix(round((M_interp-M_orig+1)/2:(M_interp-M_orig)/2+M_orig),round((M_interp-M_orig+1)/2:(M_interp-M_orig)/2+M_orig))=X_orig;
            X_interp=nmatrix;

        elseif interpM==interpL && interpM<1 %&& mod(round(1/interpM),2) ~= 0


             if mod(round(interpM*M_orig),2) == 0
                M_interp=round(interpM*M_orig)+1;
             else
                 M_interp=round(interpM*M_orig);
             end


            X_interp=X_orig(round((M_orig-M_interp)/2+1:(M_orig-M_interp)/2+M_interp),round((M_orig-M_interp)/2+1:(M_orig-M_interp)/2+M_interp));

        else

            L_interp=L_orig*interpL;

            if mod(round(interpM*M_orig),2) == 0
                M_interp=round(interpM*M_orig)+1;
            else
                M_interp=round(interpM*M_orig);
            end


            if interpL<1
                Po=sum(sum(abs(X_orig(round(M_orig*(L_orig-L_interp)/2/L_orig):round(M_orig*(1-(L_orig-L_interp)/2/L_orig)),round(M_orig*(L_orig-L_interp)/2/L_orig):    round(M_orig*(1-(L_orig-L_interp)/2/L_orig)))).^2));
            else
                Po=sum(sum(abs(X_orig).^2));
            end

            [xxo,yyo]=meshgrid((M_orig-1)/2+1-(1:M_orig));
            dxo=L_orig/M_orig;

            [xxi,yyi]=meshgrid((M_interp-1)/2+1-(1:M_interp));
            dxi=L_interp/M_interp;

            X_interp=interp2(xxo*dxo,yyo*dxo,X_orig,xxi*dxi,yyi*dxi,method,min(min(X_orig)));
            Pi=sum(sum(abs(X_interp).^2));
            
            if Pi~=0
            X_interp=X_interp*sqrt(Po/Pi);
            end

        end
    
else
    
    L_interp=interpL;
    M_interp=round(interpM);


        if L_interp<L_orig
            Po=sum(sum(abs(X_orig(round(M_orig*(L_orig-L_interp)/2/L_orig):round(M_orig*(1-(L_orig-L_interp)/2/L_orig)),round(M_orig*(L_orig-L_interp)/2/L_orig):    round(M_orig*(1-(L_orig-L_interp)/2/L_orig)))).^2));
        else
            Po=sum(sum(abs(X_orig).^2));
        end

        [xxo,yyo]=meshgrid((M_orig-1)/2+1-(1:M_orig));
        dxo=L_orig/M_orig;

        [xxi,yyi]=meshgrid((M_interp-1)/2+1-(1:M_interp));
        dxi=L_interp/M_interp;

        X_interp=interp2(xxo*dxo,yyo*dxo,X_orig,xxi*dxi,yyi*dxi,method,min(min(X_orig)));
        Pi=sum(sum(abs(X_interp).^2));
        
        if Pi~=0
        X_interp=X_interp*sqrt(Po/Pi);
        end
        
end