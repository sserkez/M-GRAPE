function[u1]=prop_TF_a(u1,Lx,Ly,lambda,z)

% by Svitozar Serkez

    if z==0
    else
        %Propagation-Transfer Function;

        % u1     - input field matrix A+iB
        % L      - trasverse length of field (m)
        % lambda - wavelength (m)
        % z      - propagation distance
        [My,Mx,N]=size(u1);
        dx=Lx/Mx;
        dy=Ly/My;
        %u1=gpuArray(u1);
        fx=linspace(-1/(2*dx)+1/Lx/2,1/(2*dx)-1/Lx/2,Mx);
        fy=linspace(-1/(2*dy)+1/Ly/2,1/(2*dy)-1/Ly/2,My);
        %fx= gpuArray.linspace(-1/(2*dx)+1/L/2,1/(2*dx)-1/L/2,M);
        [FX,FY]=meshgrid(fx,fy);

        H=exp(-1j*pi*lambda*z*(FX.^2+FY.^2));

        if N==1

             u1=fftshift(fft2(ifftshift(u1)));
             %fieldplot(301,u1,1,'field before',1)
             %fieldplot(302,H,1,'propagator',1)
 
             
             u1=H.*u1;                    
             %fieldplot(303,u1,1,'field after',1) %!original
             u1=fftshift(ifft2(ifftshift(u1)));

        else

        %     %u1=ifftshift(ifftshift(u1,1),2);
        %     u1=fft2(u1);
        %     u1=fftshift(fftshift(u1,1),2);
        % 
        %     for i=1:N
        %         u1(:,:,i)=H.*u1(:,:,i);
        %     end
        % 
        %     u1=ifftshift(ifftshift(u1,1),2);
        %     u1=ifft2(u1);
        %     %u1=fftshift(fftshift(u1,1),2);

            for i=1:N
                u=u1(:,:,i);
                u=fftshift(fft2(u));
                u=H.*u;
                u=ifft2(ifftshift(u));
                u1(:,:,i)=u;
            end

        end
    end
%u2=gather(u2);
end