% creates gaussian field file
% M-number of points in mesh x or y (preferred odd number)
% sx, sy - desired intensity rms of waist [m]
% wx, wy - desired amplitude rms of waist [m]
% zx, zy - desired position of waist (negative - downstream, positive -
% upstream) [m]
% lambda - wavelength [m]

function X=fieldgaussian(M,leng,sigmint_x,sigmint_y,zx,zy,lambda,P)
    
    if mod(round(M),2) == 0
        M=M-1;
    else
        M=M;
    end

    dx=leng/M;
%     wx=sigmint_x*sqrt(2);
%     wy=sigmint_y*sqrt(2);
%     qx=(2j*pi*(wx^2)/lambda+zx);
%     qy=(2j*pi*(wy^2)/lambda+zy);
    wx=sigmint_x*2;
    wy=sigmint_y*2;
    qx=(1j*pi*(wx^2)/lambda+zx);
    qy=(1j*pi*(wy^2)/lambda+zy);
    [xx,yy]=meshgrid((M-1)/2+1-(1:M));
    K=2*pi/lambda;
    %disp(xx.*dx);
    X=single(exp(-1j.*K.*((xx.*dx).^2./qx+(yy*dx).^2./qy)./2));
     Pi=sum(sum(abs(X).^2));
     X=X*sqrt(P/Pi);