function distcollimate(nm_0,nm_1,Nsigm)

dist=distread(nm_0);

%T_m=[-10 10];

distplot(1631,dist);
%Nsigm=3;
%%


% T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));
% 
% dist.X=dist.X-mean(dist.X(T_m_i));
% dist.Y=dist.Y-mean(dist.Y(T_m_i));
% dist.PX=dist.PX-mean(dist.PX(T_m_i));
% dist.PY=dist.PY-mean(dist.PY(T_m_i));
% 
% dist=distcut(dist,[-150e-6 150e-6],[-2e-5 2e-5],[-100e-6 100e-6],[-2e-5 2e-5],[],[],'leave');
% 
% T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));
% 
% dist.X=dist.X-mean(dist.X(T_m_i));
% dist.Y=dist.Y-mean(dist.Y(T_m_i));
% dist.PX=dist.PX-mean(dist.PX(T_m_i));
% dist.PY=dist.PY-mean(dist.PY(T_m_i));
% 
% dist.X=dist.X.*emittance_scale.^2;
% dist.PX=dist.PX.*emittance_scale.^2;
% dist.Y=dist.Y.*emittance_scale.^2;
% dist.PY=dist.PY.*emittance_scale.^2;
% distplot(1632,dist);
% distplot(16321,dist,T_m_i);

%%
%T_m=[-10 10];
% T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));
% emitx_m=sqrt(mean(dist.X(T_m_i).^2).*mean(dist.PX(T_m_i).^2)-mean(dist.X(T_m_i).*dist.PX(T_m_i)).^2).*mean(dist.G(T_m_i));
% emity_m=sqrt(mean(dist.Y(T_m_i).^2).*mean(dist.PY(T_m_i).^2)-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).^2).*mean(dist.G(T_m_i));
% betax_m=mean(dist.X(T_m_i).*dist.X(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
% betay_m=mean(dist.Y(T_m_i).*dist.Y(T_m_i)).*mean(dist.G(T_m_i))./emity_m;
% alphax_m=-mean(dist.X(T_m_i).*dist.PX(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
% alphay_m=-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).*mean(dist.G(T_m_i))./emity_m;


%remove correlation
dist.PX=dist.PX+dist.alphax.*dist.X./dist.betax;
dist.PY=dist.PY+dist.alphay.*dist.Y./dist.betay;

%mean(dist.X.*dist.X).*mean(dist.G)./emitx_m

distplot(1633,dist);

% %scale the beam
% dist.X=dist.X.*sqrt(betax_new./betax_m);
% dist.Y=dist.Y.*sqrt(betay_new./betay_m);
% dist.PX=dist.PX.*sqrt(betax_m./betax_new);
% dist.PY=dist.PY.*sqrt(betay_m./betay_new);

%mean(dist.X.*dist.X).*mean(dist.G)./emitx

dist=distcut(dist,[mean(dist.X)-Nsigm*std(dist.X) mean(dist.X)+Nsigm*std(dist.X)],[mean(dist.PX)-Nsigm*std(dist.PX) mean(dist.PX)+Nsigm*std(dist.PX)],[mean(dist.Y)-Nsigm*std(dist.Y) mean(dist.Y)+Nsigm*std(dist.Y)],[mean(dist.PY)-Nsigm*std(dist.PY) mean(dist.PY)+Nsigm*std(dist.PY)],[],[],'leave');

distplot(1634,dist);

%add correlation
dist.PX=dist.PX-dist.alphax.*dist.X./dist.betax;
dist.PY=dist.PY-dist.alphay.*dist.Y./dist.betay;

%mean(dist.X.*dist.X).*mean(dist.G)./emitx



% emitx=sqrt(mean(dist.X.^2).*mean(dist.PX.^2)-mean(dist.X.*dist.PX).^2).*mean(dist.G);
% emity=sqrt(mean(dist.Y.^2).*mean(dist.PY.^2)-mean(dist.Y.*dist.PY).^2).*mean(dist.G);
% betax=mean(dist.X.*dist.X).*mean(dist.G)./emitx;
% betay=mean(dist.Y.*dist.Y).*mean(dist.G)./emity;
% alphax=-mean(dist.X.*dist.PX).*mean(dist.G)./emitx;
% alphay=-mean(dist.Y.*dist.PY).*mean(dist.G)./emity;
% 
% T_m_i=round(find(dist.T>=T_m(1) & dist.T<=T_m(2)));
% 
% emitx_m=sqrt(mean(dist.X(T_m_i).^2).*mean(dist.PX(T_m_i).^2)-mean(dist.X(T_m_i).*dist.PX(T_m_i)).^2).*mean(dist.G(T_m_i));
% emity_m=sqrt(mean(dist.Y(T_m_i).^2).*mean(dist.PY(T_m_i).^2)-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).^2).*mean(dist.G(T_m_i));
% betax_m=mean(dist.X(T_m_i).*dist.X(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
% betay_m=mean(dist.Y(T_m_i).*dist.Y(T_m_i)).*mean(dist.G(T_m_i))./emity_m;
% alphax_m=-mean(dist.X(T_m_i).*dist.PX(T_m_i)).*mean(dist.G(T_m_i))./emitx_m;
% alphay_m=-mean(dist.Y(T_m_i).*dist.PY(T_m_i)).*mean(dist.G(T_m_i))./emity_m;
% 
% dist=distcut(dist,[-100e-6 100e-6],[-2e-5 2e-5],[-100e-6 100e-6],[-2e-5 2e-5],[],[]);
%  dist=distcut(dist,[],[],[],[],[1e-15 40e-15],[],'delete');
%  dist=distcut(dist,[],[],[],[],[70e-15 95e-15],[],'delete');

distplot(1635,dist);

%%
distwrite(dist,nm_1);
%distwrite(dist,'C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\collim_m_central.dist');

