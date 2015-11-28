% x = -1000:1:1000;
% s = 2;
% m = 3;
% y = 1/(sqrt(2*pi)* s ) * exp( - (x-m).^2 / (2*s^2)) + 0.02*randn( 1, 2001 );


x=(d(1).outp.Lamdscale-1.239e-9)*1e9*1e5;
x=(d(1).outp.Lamdscale)*1e9;
y=d(1).outp.spectrum_mid.v(:,500);
y=y./max(y);
% x(y<0.01)=[];
% y(y<0.01)=[];
tic
% [sigma,mu] = gaussfit( x, y, 1e-4, 1);
[sigma,mu] = gaussfit( x, y);
toc
xp = x;
yp =  exp( - (xp-mu).^2 / (2*sigma^2));
figure(787)
plot( x, y, '.', xp, yp, '-' ); 

%%
% clear all
% close all
distname='C:\-D-\Work\LCLS\tmp\3\ebeams\better\collim_matched';
distname='C:\-D-\Work\LCLS\tmp\3\ebeams\better\146pC1p3kA3p6GeVSF_lcut.dist';
distname='C:\-D-\Work\LCLS\tmp\3\ebeams\150pC1p5kA3p5GeV.dist';

distname='C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u1_1.out.dat';
distname='C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_1.out.dat';

distname='C:\-D-\Work\SASE3_SXRSS\beam_0.1nC.txt';
%distname='C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\slotted.dist_cut';
%distname='C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\slotted.dist_cut';
dist=distread(distname);
% dist.T=dist.T-min(dist.T);
distplot(1287,dist,:);
%%
dist=distcut(dist,[],[],[],[],[40e-15 70e-15],[]);
%%
distplot(dist);
% distcurexp(dist,[distname,'.curr'])
% drawnow
% 

% genesis_rw_wake([distname,'.curr'],[distname,'.wakexx']);
%%

% dist0=distread('C:\-D-\Work\LCLS\tmp\3\ebeams\better\0\slotted.dist');
dist0=distread('C:\-D-\Work\LCLS\tmp\3\ebeams\better\0\collim.dist');
dist0.G=dist0.G-mean(dist0.G)+8954;
distplot(dist0);
dist=dist0;
distcut;
%distwrite(dist,'C:\-D-\Work\LCLS\tmp\3\ebeams\better\0\slotted.dist_0');
dist0=dist;
distplot(dist0)
break
%%
dist=dist0;

emittance_scale_ar=[1.0 1.05 1.1 1.15];
gamma_scale_ar=[1.0 1.1 1.2 1.4];

for i=1:numel(emittance_root_scale_ar)
    for j=1:numel(gamma_scale_ar)

        emittance_scale=emittance_scale_ar(i);
        gamma_scale=gamma_scale_ar(j);

        nmnm=['C:\-D-\Work\LCLS\tmp\3\ebeams\better\0\dist_slotted_930_em',num2str(emittance_scale),'_ga',num2str(gamma_scale),'.dist'];

        dist.G=(dist.G-mean(dist.G))*gamma_scale+mean(dist.G);
        %distwrite(dist,'C:\-D-\Work\LCLS\tmp\3\ebeams\better\collim_930.dist');



        % dist=distread('C:\-D-\Work\LCLS\tmp\3\930_tdp\collim_930.dist');
        % dist=distread('C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u1_7.out.dat');

        %dist.G=dist.G-mean(dist.G)+8954;


        dist.X=dist.X.*emittance_scale.^2;
        dist.PX=dist.PX.*emittance_scale.^2;
        dist.Y=dist.Y.*emittance_scale.^2;
        dist.PY=dist.PY.*emittance_scale.^2;

        gavg=mean(dist.G);

        dist.emitx=sqrt(mean(dist.X.^2).*mean(dist.PX.^2)-mean(dist.X.*dist.PX).^2).*gavg;
        dist.emity=sqrt(mean(dist.Y.^2).*mean(dist.PY.^2)-mean(dist.Y.*dist.PY).^2).*gavg;
        dist.betax=mean(dist.X.*dist.X).*gavg./dist.emitx;
        dist.betay=mean(dist.Y.*dist.Y).*gavg./dist.emity;
        dist.alphax=-mean(dist.X.*dist.PX).*gavg./dist.emitx;
        dist.alphay=-mean(dist.Y.*dist.PY).*gavg./dist.emity;

        distwrite(dist,nmnm);
    end
end
break
%%
% [X,N]=fieldimport_all('C:\-D-\Work\LCLS\tmp\3\damage\500_u1_1.out.dfl',151,1);
% fieldplot3d(11111,X,12e-4,1,1:1:N,'name',1);

[X,N]=fieldimport_all('C:\-D-\Work\LCLS\tmp\3\damage\500_u1_1.out_f1.dfl',201,1);
fieldplot3d(11112,X,12e-4,1,1:1:N,'name',1);

[X,N]=fieldimport_all('C:\-D-\Work\LCLS\tmp\3\damage\500_u1_1.out_f2.dfl',201,1);
fieldplot3d(11113,X,12e-4,1,1:1:N,'name',1);
break
%%
% E=[1.9e-5 8.2e-6 6.7e-6 1.54e-6 1.4e-6;
%     1.24e-5 5.43e-6 4.26e-6 1.58e-6 6.7e-7;
%     7.68e-6 3.41e-6 2.62e-6 1.9e-7 1.7e-7]';
% g=[1.0 1.2 1.4];
% e=[1.0 1.05 1.1 1.2 1.4];

E=[1.54e-5 2.05e-6 6.8207e-08  4.3098e-08;
    2e-5 9.29e-7 5.5905e-08 4.1909e-08;
    1.2e-5 2.82e-7 4.8460e-08 4.0385e-08;
    2.74e-6 7.78e-8 4.3903e-08 3.8078e-08]';
g=[1.0 1.1 1.2 1.4];
e=[1.0 1.05 1.1 1.2];
figure(20314)
surf(g,e,E);
xlabel('\Delta\gamma_{scale}');
ylabel('\epsilon_{scale}');
zlabel('Power [a.u.]');
view([1,1,1])
break
%%

for i=1:9
      %nm_p{i}=['C:\-D-\Work\SASE3_chicane\run6\U1.',num2str(i),'.out'];
      %nm_p{i}=['C:\-D-\Work\SASE3_chicane\run6\U2.',num2str(i),'.out'];
      nm_p=['/nfs/slac/g/beamphysics/svitozar/950_tdp/Slotted_mm/930_sl_u1_',num2str(i),'.out'];
      %nm_p=['/data/netapp/xfel/svitozar/LCLS/run3/930_tdp/Slotted_mm/930_sl_u1_',num2str(i),'.out'];
      %nm_p=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u1_',num2str(i),'.out'];
      nm_f=[nm_p,'_f.dfl'];
        d=outread(nm_p,1); 
      [X,N]=fieldimport_all(nm_f,d.inp.ncar,1);
      
              %scale the power!  notice! notice! notice! notice! notice! notice! %notice!
        % (sum(sum(sum(abs(X).^2)))/N/3e8)/(sum(sum(sum(abs(X./5).^2)))/N/3e8)
        X=X./sqrt(5);
        %
        fieldexport(X,[nm_p,'_f5.dfl']);
      
end

break
%%
%dist=dist0;
clear all
close all

emittance_scale=1.0;

betax_new=1.498871e+01/930*1200;
betay_new=6.560118e+00/930*1200;
alphax_new=1.471605e+00;
alphay_new=-6.855406e-01;

%dist0=distread('C:\-D-\Work\LCLS\tmp\3\ebeams\better\0\slotted.dist');
%dist0=distread('c:\-D-\Work\LCLS\tmp\dist_1000.dat');
%dist0=distread('C:\-D-\Work\LCLS\tmp\3\s-n study\1\slotted_mm_1.03.dist');

dist0=distread('C:\-D-\Work\LCLS\tmp\3\ebeams\150pC1p5kA3p5GeV.dist');
%delete this line!
%dist0=distread('C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_1.out.dat');
%dist0=distread('C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\slotted.dist_cut');
betax_new=9.419507787376498;
betay_new=23.003246071749620;
alphax_new=-0.467215348689996;
alphay_new=2.756895075987282;

T_m=[50e-15 55e-15];
T_m=[-10 10];


%dist0.G=dist0.G-mean(dist0.G)+10170;
distplot(1631,dist0);

%%
dist=dist0;

T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));

dist.X=dist.X-mean(dist.X(T_m_i));
dist.Y=dist.Y-mean(dist.Y(T_m_i));
dist.PX=dist.PX-mean(dist.PX(T_m_i));
dist.PY=dist.PY-mean(dist.PY(T_m_i));

dist=distcut(dist,[-150e-6 150e-6],[-2e-5 2e-5],[-100e-6 100e-6],[-2e-5 2e-5],[],[],'leave');

T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));

dist.X=dist.X-mean(dist.X(T_m_i));
dist.Y=dist.Y-mean(dist.Y(T_m_i));
dist.PX=dist.PX-mean(dist.PX(T_m_i));
dist.PY=dist.PY-mean(dist.PY(T_m_i));

dist.X=dist.X.*emittance_scale.^2;
dist.PX=dist.PX.*emittance_scale.^2;
dist.Y=dist.Y.*emittance_scale.^2;
dist.PY=dist.PY.*emittance_scale.^2;
distplot(1632,dist);
distplot(16321,dist,T_m_i);

%%
%T_m=[-10 10];
T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));
emitx_m=sqrt(mean(dist.X(T_m_i).^2).*mean(dist.PX(T_m_i).^2)-mean(dist.X(T_m_i).*dist.PX(T_m_i)).^2).*mean(dist.G(T_m_i));
emity_m=sqrt(mean(dist.Y(T_m_i).^2).*mean(dist.PY(T_m_i).^2)-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).^2).*mean(dist.G(T_m_i));
betax_m=mean(dist.X(T_m_i).*dist.X(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
betay_m=mean(dist.Y(T_m_i).*dist.Y(T_m_i)).*mean(dist.G(T_m_i))./emity_m;
alphax_m=-mean(dist.X(T_m_i).*dist.PX(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
alphay_m=-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).*mean(dist.G(T_m_i))./emity_m;


%remove old correlation
dist.PX=dist.PX+alphax_m.*dist.X./betax_m;
dist.PY=dist.PY+alphay_m.*dist.Y./betay_m;

%mean(dist.X.*dist.X).*mean(dist.G)./emitx_m

distplot(1633,dist);
%scale the beam
dist.X=dist.X.*sqrt(betax_new./betax_m);
dist.Y=dist.Y.*sqrt(betay_new./betay_m);
dist.PX=dist.PX.*sqrt(betax_m./betax_new);
dist.PY=dist.PY.*sqrt(betay_m./betay_new);

%mean(dist.X.*dist.X).*mean(dist.G)./emitx

dist=distcut(dist,[-150e-6 150e-6],[-1e-5 1e-5],[-150e-6 150e-6],[-1.5e-5 1.5e-5],[],[],'leave');

distplot(1634,dist);

%add new correlation
dist.PX=dist.PX-alphax_new.*dist.X./betax_new;
dist.PY=dist.PY-alphay_new.*dist.Y./betay_new;

%mean(dist.X.*dist.X).*mean(dist.G)./emitx



emitx=sqrt(mean(dist.X.^2).*mean(dist.PX.^2)-mean(dist.X.*dist.PX).^2).*mean(dist.G);
emity=sqrt(mean(dist.Y.^2).*mean(dist.PY.^2)-mean(dist.Y.*dist.PY).^2).*mean(dist.G);
betax=mean(dist.X.*dist.X).*mean(dist.G)./emitx;
betay=mean(dist.Y.*dist.Y).*mean(dist.G)./emity;
alphax=-mean(dist.X.*dist.PX).*mean(dist.G)./emitx;
alphay=-mean(dist.Y.*dist.PY).*mean(dist.G)./emity;

T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));

emitx_m=sqrt(mean(dist.X(T_m_i).^2).*mean(dist.PX(T_m_i).^2)-mean(dist.X(T_m_i).*dist.PX(T_m_i)).^2).*mean(dist.G(T_m_i));
emity_m=sqrt(mean(dist.Y(T_m_i).^2).*mean(dist.PY(T_m_i).^2)-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).^2).*mean(dist.G(T_m_i));
betax_m=mean(dist.X(T_m_i).*dist.X(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
betay_m=mean(dist.Y(T_m_i).*dist.Y(T_m_i)).*mean(dist.G(T_m_i))./emity_m;
alphax_m=-mean(dist.X(T_m_i).*dist.PX(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
alphay_m=-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).*mean(dist.G(T_m_i))./emity_m;

%dist=distcut(dist,[-100e-6 100e-6],[-2e-5 2e-5],[-100e-6 100e-6],[-2e-5 2e-5],[],[]);
 dist=distcut(dist,[],[],[],[],[1e-15 40e-15],[],'delete');
 dist=distcut(dist,[],[],[],[],[70e-15 95e-15],[],'delete');

distplot(1635,dist);
%%
distwrite(dist,'C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\slotted.dist_cut_rematched_cut');
%distwrite(dist,'C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\collim_m_central.dist');

