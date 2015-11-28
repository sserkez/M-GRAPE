%analyze dfls at a distance
d=outread('/data/netapp/xfel/svitozar/xfel/sase3_SXRSS/tdp_4/U2_t.1.out',0,0,1);
ncar=201;
xlamds=d.inp.xlamds;
zsep=d.inp.zsep;
nslice=d.inp.nslice;
leng=d.inp.leng;
for iii=1:5
nm_dfl=['/data/netapp/xfel/svitozar/xfel/sase3_SXRSS/tdp_4/U1.',num2str(iii),'.out_f7.dfl'];
disp(iii);
        disp(' -import of .dfl field file,');
        tic
        [Xt,dfl.nslice]=fieldimport_all(nm_dfl,ncar,1);
        t=toc;
%       disp('  +done');
        dfl.Sscale=linspace(0,xlamds*zsep*nslice,nslice);
        dfl.E=sum(sum(sum(abs(Xt).^2)))/nslice*abs(dfl.Sscale(end)-dfl.Sscale(1))/3e8;
        dfl.power=reshape(sum(sum(abs(Xt).^2,1),2),1,[]);
        [Xf,dfl.Lamdscale]=dfl_time2freq(Xt,dfl.Sscale,xlamds);
        clear Xt
        dfl.spectrum=reshape(sum(sum(abs(Xf).^2,1),2),1,[]);
                disp(' -backpropagation of .dfl'); 
        Xf=prop_TF(Xf,leng,xlamds,-1.1);

        disp(' -calculation of .dfl field transverse distribution');        
        
        Xfz=fftshift(fft2(ifftshift(Xf)));

        dfl.X=sqrt(sum(abs(Xf).^2,3)).*exp(1i.*angle(mean(Xf,3))); %2D field
        clear Xf
        dfl.Xfz=sqrt(sum(abs(Xfz).^2,3)).*exp(1i.*angle(mean(Xfz,3)));
        clear Xfz       
        disp('   +done');

ext(iii)=dfl;
disp(' ');
end
%%
save('ext','ext');
break
%%
    figure(5959)
    hold on
for iii=1:4
    plot(1239.8./ext(iii).Lamdscale/1e9-698.45,ext(iii).spectrum);
    %plot(ext(iii).Lamdscale,ext(iii).spectrum)
end
hold off