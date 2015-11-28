% creates gaussian field file
% M-number of points in mesh x or y (preferred odd number)
% sx, sy - desired intensity rms of waist [m]
% wx, wy - desired amplitude rms of waist [m]
% zx, zy - desired position of waist (negative - downstream, positive -
% upstream) [m]
% lambda - wavelength [m]

function X=fieldgaussian_a(Mx,My,Lx,Ly,sigmint_x,sigmint_y,zx,zy,lambda,P)
    
    if mod(round(Mx),2) == 0
        Mx=Mx-1;
    end
    
    if mod(round(My),2) == 0
        My=My-1;
    end

    dx=Lx/Mx;
    dy=Ly/My;
%     wx=sigmint_x*sqrt(2);
%     wy=sigmint_y*sqrt(2);
%     qx=(2j*pi*(wx^2)/lambda+zx);
%     qy=(2j*pi*(wy^2)/lambda+zy);
    wx=sigmint_x*2;
    wy=sigmint_y*2;
    qx=(1j*pi*(wx^2)/lambda+zx);
    qy=(1j*pi*(wy^2)/lambda+zy);
    %[yy,xx]=meshgrid((My-1)/2+1-(1:My),(Mx-1)/2+1-(1:Mx));
       [xx,yy]=meshgrid((Mx-1)/2+1-(1:Mx),(My-1)/2+1-(1:My));
    K=2*pi/lambda;
    %disp(xx.*dx);
    X=single(exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dy).^2./qy)./2));
     Pi=sum(sum(abs(X).^2));
     X=X*sqrt(P/Pi);