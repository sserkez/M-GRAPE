%%
% nm_0='c:\-D-\Work\LCLS\Thesis_source\geninp_1200_u1';
% nm='c:\-D-\Work\LCLS\Thesis_source\geninp';
% copyfile(nm_0,nm)
% % %copyfile('300.in','geninp')
% % xlamds_arr=1.77e-9:0.002e-9:1.78e-9;
% % 
% % P(1:numel(xlamds_arr))=0;
% % for k=1:numel(xlamds_arr)
% % fclose all;
% % delete([nm_out,'.dfl']);
% % fid = fopen(['c:\-D-\Work\LCLS\Thesis_source\geninp'],'w');
% % fd = fopen(nm,'r');
% % 
% % 	%fd=fopen(nm0,'r');
% % 	while 1
% %     		tline = fgetl(fd);
% %     		if ~ischar(tline), break, end
% %     		if findstr('xlamds',tline)
% % 			fprintf(fid,' xlamds = %d\n',xlamds_arr(k));
% %             else
% %                 fprintf(fid,[tline,'\n']);
% % 			continue
% %             end
% %     end
% % fclose all;
% 
% dos('genesis301');
%% reads .out file
%fclose all;
clear variables
clear all
N=1;


%      nm_p{1}='C:\-D-\Work\LCLS\tmp\1\1000_u1_tdp.out';
%       nm_p{1}='C:\-D-\Work\LCLS\tmp\1\500_u1_tdp.out';
%       nm_p{1}='C:\-D-\Work\LCLS\tmp\1\500_u2_tdp_s.out';
%       nm_p{2}='C:\-D-\Work\LCLS\tmp\1\500_u2_tdp_n.out';
%       nm_p{3}='C:\-D-\Work\LCLS\tmp\1\500_u2_tdp_f.out';
%       nm_p{1}='C:\-D-\Work\LCLS\tmp\1\1000_u2_tdp_s.out';
%       nm_p{2}='C:\-D-\Work\LCLS\tmp\1\1000_u2_tdp_n.out';
%       nm_p{3}='C:\-D-\Work\LCLS\tmp\1\1000_u2_tdp_f.out';
      
%         nm_p{1}='C:\-D-\Work\LCLS\tmp\1000_u1_tdp_1.1.out';
%        nm_p{2}='C:\-D-\Work\LCLS\tmp\1000_u2_tdp_1.1.out';
%        
%         nm_p{1}='C:\-D-\Work\LCLS\tmp\1\1000_u1_tdp.out';
%        nm_p{1}='C:\-D-\Work\LCLS\tmp\2\1000_u2_tdp_f.out';
        %nm_p{1}='C:\-D-\Work\LCLS\tmp\2\1000_u2_tdp_f_t.out';
%        nm_p{2}='C:\-D-\Work\LCLS\tmp\2\1000_u2_tdp_f_t1.out';
%        nm_p{2}='C:\-D-\Work\LCLS\tmp\2\1000_u2_tdp_fhorn.out';
 
%nm_p{1}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u1_tdp_1.out';      
%nm_p{1}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u2_tdp_s1.out';
%nm_p{2}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u2_tdp_1.out';
%nm_p{3}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u2_tdp_f1.out';

% nm_p{1}='C:\-D-\Work\LCLS\tmp\3\530_tdp\530_u1_tdp_lcut_1.out';
 %nm_p{2}='C:\-D-\Work\LCLS\tmp\3\530_u1_tdp_lcut_propwake.out';

% nm_p{2}='C:\-D-\Work\LCLS\tmp\3\stage1.out';
%  nm_p{7}='C:\-D-\Work\LCLS\tmp\3\530_tdp\530_u2_tdp_1.out';
%   nm_p{1}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u2_tdp_1.out';
 % nm_p{2}='C:\-D-\Work\LCLS\tmp\3\s-n study\530_u2_tdp_s.out';

% nm_p{1}='C:\-D-\Work\LCLS\tmp\3\950_tdp\2\950_u1_slotted.out';
% nm_p{2}='C:\-D-\Work\LCLS\tmp\3\950_tdp\2\950_u1_collim.out';

% nm_p{1}='C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u2_1f1.out';
%nm_p{1}='C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u2_1.out';

% nm_p{1}='C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u1_1_1.0.out';
% nm_p{2}='C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u1_1_1.1.out';
% nm_p{3}='C:\-D-\Work\LCLS\tmp\3\930_tdp\930_u1_1_1.2.out';

% nm_p{1}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1_ga1.out';
% nm_p{2}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1_ga1.1.out';
% nm_p{3}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1_ga1.2.out';
% nm_p{4}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1_ga1.4.out';
% nm_p{5}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.05_ga1.out';
% nm_p{6}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.05_ga1.1.out';
% nm_p{7}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.05_ga1.2.out';
% nm_p{8}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.05_ga1.4.out';
% nm_p{9}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.1_ga1.out';
% nm_p{10}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.1_ga1.1.out';
% nm_p{11}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.1_ga1.2.out';
% nm_p{12}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.1_ga1.4.out';
% nm_p{13}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.2_ga1.out';
% nm_p{14}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.2_ga1.1.out';
% nm_p{15}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.2_ga1.2.out';
% nm_p{16}='C:\-D-\Work\LCLS\tmp\3\beam_study\2\sl_em1.2_ga1.4.out';


 %nm_p{1}='C:\-D-\Work\LCLS\tmp\3\beam_study\cut_u-detune\cut_notdet.out';
 %nm_p{1}='C:\-D-\Work\LCLS\tmp\3\beam_study\cut_u-detune\cut_det.out';
%   nm_p{2}='C:\-D-\Work\LCLS\tmp\3\beam_study\cut_u-detune\cut_det_smgrid.out';
%   nm_p{1}='C:\-D-\Work\LCLS\tmp\3\beam_study\cut_u-detune\cut_notdet_smgrid.out';
 %nm_p{2}='C:\-D-\Work\LCLS\tmp\3\beam_study\cut_u-detune\notcut_det.out';
 
%   nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u2_2.out'];
%   nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u2_28_2.out'];
  %nm_p{3}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u2_28up_2.out'];
  
 %nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_1.out'];
 %nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted\930_sl_u2_2.out'];

%  nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_1.out'];
%  nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u2_1.out'];
%  nm_p{3}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_1.out'];
%  nm_p{4}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Collimated_matched\930_cl_u1_1.out'];

 %nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_1.out'];
% nm_p{3}=['C:\-D-\Work\LCLS\tmp\3\beam_study\central\930_sl_u1_central.out'];

%nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\damage\300_u1_1.out'];

%nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_1_sc10.out'];
nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_fc_nn_ps_1_sc10.out'];
%nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_mb_nn_ps_1_sc10.out'];
%nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_mm\930_sl_u1_1_x.out'];
%nm_p{3}=['C:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u2_1_sc10.out'];

%nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\300_u2_f.out'];
%nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\300_u2_n.out'];

% nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\300_u1.out'];
% nm_p{2}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\300_u2_n.out'];

% nm_p{1}='c:\-D-\Work\SASE3_SXRSS\700_u1_tdp.out';
% nm_p{1}='c:\-D-\Work\SASE3_SXRSS\300_u1_tdp.out';
% nm_p{2}='c:\-D-\Work\SASE3_SXRSS\500_u1.out';
% nm_p{3}='c:\-D-\Work\SASE3_SXRSS\700_u1.out';
% nm_p{4}='c:\-D-\Work\SASE3_SXRSS\1000_u1.out';
%nm_p{2}='c:\-D-\Work\SASE3_SXRSS\u1_700.out';

% nm_p{1}='c:\-D-\Work\LCLS\Thesis_source\300_u1_1.out';
% nm_p{2}='c:\-D-\Work\LCLS\Thesis_source\500_u1_1.out';
% nm_p{3}='c:\-D-\Work\LCLS\Thesis_source\700_u1_1.out';
% nm_p{4}='c:\-D-\Work\LCLS\Thesis_source\1200_u1_1.out';

% nm_p{1}='c:\-D-\Work\SASE3_SXRSS\300_u1_tdp.out';
% nm_p{2}='c:\-D-\Work\SASE3_SXRSS\500_u1_tdp.out';
% nm_p{3}='c:\-D-\Work\SASE3_SXRSS\700_u1_tdp.out';
% nm_p{4}='c:\-D-\Work\SASE3_SXRSS\1000_u1_tdp.out';

nm_p{1}='c:\-D-\Work\SASE3_SXRSS\1000_u1.out';
%nm_p{1}='c:\-D-\Work\SASE3_SXRSS\500_u1.out';
 %nm_p{1}='c:\-D-\Work\SASE3_SXRSS\700_u1.out';
% nm_p{1}='c:\-D-\Work\SASE3_SXRSS\1000_u1.out';
%nm_p{1}='c:\-D-\Work\SASE3_SXRSS\300_u1_tdp.out';

% nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_3\700_u1_tdp_SASE.out';
% % nm_p{2}='c:\-D-\Work\SASE3_SXRSS\tdp_3\700_u1_tdp_2.out';
%  nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_3\U2_t_n.1000.out';
% % nm_p{2}='c:\-D-\Work\SASE3_SXRSS\tdp_3\U2_t_n.2.out';
% % nm_p{3}='c:\-D-\Work\SASE3_SXRSS\tdp_3\U2_t_n1.2.out';
%  nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_4\U2_t.2.out';
% % nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_3\U1.3.out';
%   nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_3_tap\U2.2_new.out';
% %nm_p{1}='c:\-D-\Work\LCLS\Last transfered\OLD\New_matlab_2013.11.05\500.out';
% %nm_p{1}='c:\-D-\Work\LCLS\tmp\3\damage\1200_u1_2.out';
%    %break    

% nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\700_u2_',num2str(1),'_f.out'];
% nm_p{2}=['C:\-D-\Work\SASE3_SXRSS\tdp_4\U2_t.',num2str(1),'.out'];

nm_p{1}=['C:\-D-\Work\LCLS\tmp\3\s-n study\1\700_u1_',num2str(1),'.out'];
nm_p{1}='c:\-D-\Work\SASE3_SXRSS\tdp_3\U1.1.out';
 nm_p{1}='C:\-D-\Work\SASE3_chicane\630_250\statrun1\U1.81.out';
 %nm_p{2}='C:\-D-\Work\SASE3_chicane\630_250\statrun1\U2.81.out';
 
 nm_p{1}='C:\-D-\Work\SASE3_chicane\run7\U1.2.out';
 nm_p{1}='D:\Work\!PROJECTS\Phase_controlled_harmonics\SASE3_v3\stage_1.out';
 %nm_p{2}='D:\Work\!PROJECTS\Phase_controlled_harmonics\SASE3_v3\stage_2.out';
 
 nm_p{2}='D:\Work\!PROJECTS\ocelot_test\test_11\run_0\run.0.s1.gout';
 %nm_p{2}=[];
 %nm_p{2}='C:\-D-\Work\SASE3_chicane\run7\U2.2.out';
% nm_p{1}='C:\-D-\Work\SASE3_chicane\630_250\Gianluca\run.0.s1.gout';
%nm_p{1}='C:\-D-\Work\SASE3_chicane\630_250\Gianluca\run.0.s3.gout';
% nm_p{1}='c:\-D-\Work\LCLS\Thesis_source\1200_u1_1.out';
% nm_p{1}='c:\-D-\Work\LCLS\tmp\3\930_tdp\Slotted_matched\930_sl_u1_1.out';
% %nm_p{1}='c:\-D-\Work\SASE3_SXRSS\optimal_waist_reimage\prop.out';
% nm_p{1}=['C:\-D-\Work\SASE3_SXRSS\tdp_4\U2_t.',num2str(1),'.out'];
 DiN=size(nm_p,2);
% DiN=1;

for Di=1:DiN %Data index
    if ~isempty(nm_p{Di})
        d(Di)=outread(nm_p{Di},1,0,2); 
    else
        d(Di)=struct;
    end
end
disp ('------------end-of-import-------------');
%% General amplification process plots
Di=1;
H{1}=outpot_e(1,d(Di));
H{2}=outpot_ph(2,d(Di));

if DiN==2
    Di=2;
    H{11}=outpot_e(11,d(Di));
    H{12}=outpot_ph(21,d(Di));
end

% span=1;
% figure(2525)
% Di=1;
% P=smooth(d(Di).outp.power.max_S',span)';
% gain_pp=log(P./circshift(P,[0,1]));
% %gain_e=log(E./circshift(E,[0,1]));
% % gain_e(isnan(gain_e))=0;
%  scatter(d(Di).outp.Zscale,(d(Di).outp.Zscale(2)-d(Di).outp.Zscale(1))./gain_pp);
%  hold all
%  Di=2;
%  P=smooth(d(Di).outp.power.max_S',span)';
% gain_pp=log(P./circshift(P,[0,1]));
% %gain_e=log(E./circshift(E,[0,1]));
% % gain_e(isnan(gain_e))=0;
%  scatter(d(Di).outp.Zscale,(d(Di).outp.Zscale(2)-d(Di).outp.Zscale(1))./gain_pp);
%  hold off
% % axis tight
% % grid on
%  ylim([0 5])
% % 
% % figure(5061)
% % plot(d(Di).outp.Zscale,2.*d(Di).outp.spectrum_mid.std_lamd./d(Di).inp.xlamds);
% % hold all
% % plot(d(Di).outp.Zscale,2.*d(Di).outp.spectrum_mid.std_lamd1./d(Di).inp.xlamds);
% % plot(d(Di).outp.Zscale,2.*d(Di).outp.spectrum_mid.std_lamd2./d(Di).inp.xlamds);
% % hold off
% % 
% % figure(5062)
% % plot(d(Di).outp.Zscale,d(Di).outp.spectrum_mid.lamdpos);
% % hold all
% % plot(d(Di).outp.Zscale,d(Di).outp.spectrum_mid.lamdpos2);
% % hold off
% %% radiation longitudinal size and position
% 
% % figure
% % plot(d(1).outp.Zscale,d(1).outp.power.std.*2)
% % hold all
% % plot(d(1).outp.Zscale,d(1).outp.power.peakpos)
% % hold off
% 
%% Bunch profile parameters at given position

Di=1;
Z=100; %[m]
H{4}=outplot_z(4,d(Di),Z);
try
    Z=100;
    H{14}=outplot_z(5,d(Di+1),Z);
catch
end
%% Wigner distribution at given position
Di=1;
Z=10; %[m]
Zi=find(d(Di).outp.Zscale<=Z,1,'last');

Ex=sqrt(d(Di).outp.p_mid.v(:,Zi)).*exp(1i*d(Di).outp.phi_mid.v(:,Zi));
% figure(6995)
% plot(d(Di).outp.Sscale,abs(Ex).^2);

% addpath('c:\-D-\Work\LCLS\New_matlab\tftb-0.2\mfiles');

% [W,~,Freqscale] = tfrwv(Ex,d(Di).outp.Zscale*1e-6/3e8);
% break
% figure(6996)
% imagesc(d(Di).outp.Zscale,Freqscale,W);
% 
W=mywigner(Ex);
figure(6996)


Sn=d(Di).outp.Sn;
sc=linspace(-Sn/2,Sn/2,Sn);%original
%sc=-(outp.Sn-1)/2:1:(outp.Sn-1)/2;

k0=2*pi/d(Di).inp.xlamds;
dk=2*pi/(d(Di).outp.Sn*d(Di).inp.xlamds*d(Di).inp.zsep);
Kscale=k0+dk*sc';
Lamdscale=2*pi./Kscale;
clear Kscale
imagesc(d(Di).outp.Sscale,Lamdscale,W'.^2);

break
%% Make GIF
%break
Di=1;
Z=0.1; %[m]d
Nframes=100;
outplot_z(4,d(Di),Z);
% Z = peaks;
% surf(Z)
% axis tight
%set(gca,'nextplot','replacechildren','visible','off')
f = getframe(4);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;
for k = 1:Nframes-1
    Z=max(d(Di).outp.Zscale)/(Nframes-1)*k;
  outplot_z(4,d(Di),Z);
  f = getframe(4);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end
%nm_p{Di}
%imwrite(im,map,[nm_p{Di},'.gif'],'DelayTime',0.2,'LoopCount',0);

%% dfl file
Di=2;
Nn=floor(d(Di).inp.nslice+d(Di).inp.zstop/d(Di).inp.xlamd/d(Di).inp.zsep);
[XX,N]=fieldimport_all([nm_p{Di},'.dfl'],d(Di).inp.ncar,1);


%XX=prop_TF(XX,d(Di).inp.dgrid*2,d(Di).inp.xlamds,0);
dflSscale=linspace(0,d(Di).inp.xlamds*d(Di).inp.zsep*N,N);

P0=sum(sum(sum(abs(XX).^2)))/N*abs(dflSscale(end)-dflSscale(1))/3e8;

fieldplot3d_a(1113,XX,d(Di).inp.leng,d(Di).inp.leng,1,dflSscale,[nm_p{Di},'.dfl'],1);
%fieldplot3d(1113,XX,d(Di).inp.leng,1,dflSscale,[nm_p{Di},'.dfl'],1);
[X2,dflLamdscale]=dfl_time2freq(XX,dflSscale,d(Di).inp.xlamds);
fieldplot3d_a(1112,X2,d(Di).inp.leng,d(Di).inp.leng,0,dflLamdscale,[nm_p{Di},'_freq.dfl'],1);
%%
U1_to_G=1.265;      %First undulator to grating [m]
deltaU=3.87;
%deltaU=0;
X=prop_TF(XX,d(Di).inp.dgrid*2,d(Di).inp.xlamds,U1_to_G+deltaU-0.6);
fieldplot3d(111,X,d(Di).inp.dgrid*2,1,dflSscale,[nm_p{Di},'.dfl'],1);
footprintplot(112,X,d(Di).inp.dgrid*2,1,'x',1,dflSscale,[nm_p{Di},'.dfl'],1);
break
%%
figure(67)
plot(reshape(sum(sum(abs(X2).^2,1),2),1,[]));
%% Source charactrization

Di=1;
disp(d(Di).nm_p);
[XX,N]=fieldimport_all([nm_p{Di},'.dfl'],d(Di).inp.ncar,1);
% XX(:,:,floor(N*2/3):end)=[];
% N=size(XX,3);
XX=prop_TF(XX,d(Di).inp.dgrid*2,d(Di).inp.xlamds,-1.09);
figure(5326)
plot(shiftdim(sum(sum(abs(XX).^2))));
fieldplot(611,XX,d(Di).inp.leng,'orig_dfl',1);

%zz=[0];
n1=300;
n2=700;
if d(Di).inp.itdp==1
    X=XX(:,:,n1:n2);
elseif d(Di).inp.itdp==0
    X=XX;
else
    error('dfl dimentions contradict itdf flag in .out');
end
interpM=1;
interpL=1;
for i=1:size(X,3)
X2(:,:,i)=fieldinterpolate_a(X(:,:,i),d(Di).inp.leng,d(Di).inp.leng,1,interpM,interpM,interpL,interpL,'spline');
end
X=single(X2); clear X2 XX
leng2=d(Di).inp.leng*interpL;

fieldplot(65,X,leng2,'undulator end',1);
    set(figure(65),'Position',[1151,222,568,563]);

%  X=prop_TF(X,leng2,d(Di).inp.xlamds,-(12-11.28));
%  waistscan_1(68,X,leng2,d(Di).inp.xlamds,zz);

Xff=fftshift(fft2(ifftshift(X)));

fieldplot_ff(66,Xff,leng2,d(Di).inp.xlamds,'far field',1)
set(figure(66),'Position',[1151+10,222+10,568+10,563+10]);
%%
zz=[-3:0.2:1];
zz=[-0.3];
%X=prop_TF(X,leng2,d(Di).inp.xlamds,-(12-11.28));
waistscan_1(68,X,leng2,d(Di).inp.xlamds,zz);


h1=figure(68);
for i=1:3
subplot(3,1,i)
xlim([zz(1) zz(end)]);
end
%%

h2=figure(69);
clf(69)
objects=allchild(h1);
copyobj(get(h1,'children'),h2);
%copyobj(objects,h2);

sigma=16e-6; %1000
z=-1; 
% sigma=17.4e-6;
% z=-0.8; %500

%650
% sigma=16e-6;
% z=-2.2; %650
%  sigma=1.3e-5;
%  z=-2;

Xg=fieldgaussian(d(Di).inp.ncar,leng2,sigma,sigma,z,z,d(Di).inp.xlamds,1);
waistscan_g(69,Xg,leng2,d(Di).inp.xlamds,zz);

set(figure(69), 'Position', [1372, 58, 540, 926]);
%[-8:0.5:3]
break
%%
sigma=1.5e-5;
z=-2.0; %650
% sigma=1.4e-5;
% z=-3.5;
X=fieldgaussian(d(Di).inp.ncar,d(Di).inp.leng,sigma,sigma,z,z,d(Di).inp.xlamds,1);
waistscan_g(X,d(Di).inp.leng,d(Di).inp.xlamds,[-6:0.5:2]);
%% Imaging characterization

%% Comparison plots

Z1=24;
Z2=24;
a1=1;
a2=2;

figure(1432)
plot(d(a1).outp.Sscale,d(a1).outp.e_spread.v(:,find(d(a1).outp.Zscale>=Z1,1,'first')),'linewidth',0.5,'linestyle','--','color','b');
hold on
plot(d(a2).outp.Sscale,d(a2).outp.e_spread.v(:,find(d(a2).outp.Zscale>=Z2,1,'first')),'linewidth',0.5,'linestyle','--','color','r');
hold off
%%
    Z=100;

% if DiN==2
% 
% 
%     for Di=1:DiN
%         Zi(Di)=find(d(Di).outp.Zscale<=Z,1,'last');
%         lbl{Di}=nm_p{Di};
%     end
% 
%     figure(56);
%     plot(d(1).outp.Lamdscale,d(1).outp.spectrum_mid.v(:,Zi(1)),'linewidth',2);%/max(d(1).outp.spectrum_mid.v(:,Zi(1)))
%     hold all
%     plot(d(2).outp.Lamdscale,d(2).outp.spectrum_mid.v(:,Zi(2)),'linewidth',2,'color','r','linestyle','--');%/max(d(2).outp.spectrum_mid.v(:,Zi(2)))
%     hold off
%     l=legend(lbl{1},lbl{2},'Location','northwest');
%     set(l, 'Interpreter', 'none','FontSize',7);
%     %legend('direct\newlineapproach','phenomenological\newlineapproach','Location','northwest')
%     %xlim([1.2385 1.2395]);
%     %xlim([1.764 1.766]);
%     ylabel('P(\lambda) [a.u.]');
%     xlabel('\lambda [nm]');
% 
%     figure(59);
%     semilogy(d(1).outp.Zscale,d(1).outp.power.max_S,'LineWidth',2);
%     hold all
%     semilogy(d(2).outp.Zscale,d(2).outp.power.max_S,'LineWidth',2,'color','r','linestyle','--');
%     hold off
%     %     legend('direct approach','phenomenological approach','Location','northwest')
%     l=legend(lbl{1},lbl{2},'Location','northwest');
%     set(l, 'Interpreter', 'none','FontSize',7);
%     ylabel('P [W]');
%     xlabel('z [m]');
%     %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
%     axis tight
% 
%     figure(591);
%     semilogy(d(1).outp.Zscale,d(1).outp.power.E,'LineWidth',2);
%     hold all
%     semilogy(d(2).outp.Zscale,d(2).outp.power.E,'LineWidth',2,'color','r','linestyle','--');
%     hold off
%     %     legend('direct approach','phenomenological approach','Location','northwest')
%     l=legend(lbl{1},lbl{2},'Location','northwest');
%     set(l, 'Interpreter', 'none','FontSize',7);
%     ylabel('E [J]');
%     xlabel('z [m]');
%     %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
%     axis tight
%     clear l
% 
% elseif DiN==3
% 
% 
%     for Di=1:DiN
%         Zi(Di)=find(d(Di).outp.Zscale<=Z,1,'last');
%         lbl{Di}=nm_p{Di};
%     end
% %lbl{1}='signal'; lbl{2}='noise'; lbl{3}='total'; lbl{4}='signal+noise';
% 
%     figure(56);  
%     plot(d(1).outp.Lamdscale,d(1).outp.spectrum_mid.v(:,Zi(1)),'linewidth',2);%/max(d(1).outp.spectrum_mid.v(:,Zi(1)))
%     hold all
%     plot(d(2).outp.Lamdscale,d(2).outp.spectrum_mid.v(:,Zi(2)),'linewidth',2,'color','r','linestyle','--');%/max(d(2).outp.spectrum_mid.v(:,Zi(2)))
%     plot(d(3).outp.Lamdscale,d(3).outp.spectrum_mid.v(:,Zi(3)),'linewidth',2,'color',[0 0.5 0],'linestyle','-.');%/max(d(2).outp.spectrum_mid.v(:,Zi(2)))
%     %plot(d(1).outp.Lamdscale,d(1).outp.spectrum_mid.v(:,Zi(1))+d(2).outp.spectrum_mid.v(:,Zi(1)),'linewidth',1,'color','k');%/max(d(1).outp.spectrum_mid.v(:,Zi(1)))
%     hold off
%     l=legend(lbl,'Location','northwest');
%     set(l, 'Interpreter', 'none');
%     %legend('direct\newlineapproach','phenomenological\newlineapproach','Location','northwest')
%     %xlim([1.2385 1.2395]);
%     %xlim([1.764 1.766]);
%     ylabel('P(\lambda) [a.u.]');
%     xlabel('\lambda [nm]');
% 
%     figure(59);
%     semilogy(d(1).outp.Zscale,d(1).outp.power.max_S,'LineWidth',2);
%     hold all
%     semilogy(d(2).outp.Zscale,d(2).outp.power.max_S,'LineWidth',2,'color','r','linestyle','--');
%     semilogy(d(3).outp.Zscale,d(3).outp.power.max_S,'LineWidth',2,'color',[0 0.5 0],'linestyle','-.');
%     %semilogy(d(1).outp.Zscale,d(1).outp.power.max_S+d(2).outp.power.max_S,'LineWidth',1,'color','k');
%     hold off
%     %     legend('direct approach','phenomenological approach','Location','northwest')
%     l=legend(lbl,'Location','southeast');
%     set(l, 'Interpreter', 'none');
%     ylabel('P [W]');
%     xlabel('z [m]');
%     %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
%     axis tight
% 
%     figure(591);
%     semilogy(d(1).outp.Zscale,d(1).outp.power.E,'LineWidth',2);
%     hold all
%     semilogy(d(2).outp.Zscale,d(2).outp.power.E,'LineWidth',2,'color','r','linestyle','--');
%     semilogy(d(3).outp.Zscale,d(3).outp.power.E,'LineWidth',2,'color',[0 0.5 0],'linestyle','-.');
%     %semilogy(d(1).outp.Zscale,d(1).outp.power.E+d(2).outp.power.E,'LineWidth',1,'color','k');
%     hold off
%     %     legend('direct approach','phenomenological approach','Location','northwest')
%     l=legend(lbl,'Location','southeast');
%     set(l, 'Interpreter', 'none');
%     ylabel('E [J]');
%     xlabel('z [m]');
%     %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
%     axis tight
%     clear l
% 
%     elseif DiN>=3

        
    for Di=1:DiN
        Zi(Di)=find(d(Di).outp.Zscale<=Z,1,'last');
        lbl{Di}=nm_p{Di};
    end
%lbl{1}='signal+noise'; lbl{2}='noise';% lbl{3}='total'; lbl{4}='signal+noise';
    
    figure(56);  
    plot(d(1).outp.Lamdscale,d(1).outp.spectrum_mid.v(:,Zi(1)),'linewidth',2);%/max(d(1).outp.spectrum_mid.v(:,Zi(1)))
    hold all
    for i=2:DiN
    plot(d(i).outp.Lamdscale,d(i).outp.spectrum_mid.v(:,Zi(i)),'linewidth',2);%/max(d(2).outp.spectrum_mid.v(:,Zi(2)))
    end
    hold off

    l=legend(lbl,'Location','northwest');
    set(l, 'Interpreter', 'none','FontSize',7);
    %legend('direct\newlineapproach','phenomenological\newlineapproach','Location','northwest')
    %xlim([1.2385 1.2395]);
    %xlim([1.764 1.766]);
    ylabel('P(\lambda) [a.u.]');
    xlabel('\lambda [nm]');

    figure(59);
    semilogy(d(1).outp.Zscale,d(1).outp.power.max_S,'LineWidth',2);
    hold all
       for i=2:DiN
    semilogy(d(i).outp.Zscale,d(i).outp.power.max_S,'LineWidth',2);
       end
       hold off
    %     legend('direct approach','phenomenological approach','Location','northwest')
    l=legend(lbl,'Location','southeast');
    set(l, 'Interpreter', 'none','FontSize',7);
    ylabel('P [W]');
    xlabel('z [m]');
    %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
    axis tight

    figure(591);
    semilogy(d(1).outp.Zscale,d(1).outp.power.E,'LineWidth',2);
    hold all
       for i=2:DiN
    semilogy(d(i).outp.Zscale,d(i).outp.power.E,'LineWidth',2);
       end
    hold off
    %     legend('direct approach','phenomenological approach','Location','northwest')
    l=legend(lbl,'Location','southeast');
    set(l, 'Interpreter', 'none','FontSize',7);
    ylabel('E [J]');
    xlabel('z [m]');
    %xlim([d(1).outp.Zscale(1) d(1).outp.Zscale(end)]);
    axis tight
    clear l
    
    
% end
%% spectrum and power evolution plots
Di=1;

if d(Di).inp.itdp~=1
    error('steady state: cannot plot spectrum & power')
end
    
    %plot pulse evolution
    Zgrid_p=repmat(d(Di).outp.Zscale,d(Di).outp.Sn,1);
    Sgrid=repmat(d(Di).outp.Sscale,1,d(Di).outp.Zn);
    Lamdgrid=repmat(d(Di).outp.Lamdscale,1,d(Di).outp.Zn);
    P=double(d(Di).outp.power.v);
    S=double(d(Di).outp.spectrum_mid.v);
%     figure(611);
%     h=surf(Zgrid,Sgrid,P,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
    Zgrid_s=repmat(d(Di).outp.Zscale,numel(d(Di).outp.Lamdscale),1);
     figure(612);
%     h=surf(Zgrid,Lamdgrid,S,'linestyle','none');
%     set(get(h,'Parent'),'ZScale','log');
%     axis tight
     Zm=10; Sm=5;
     [Zgrid_p_n, Sgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).outp.Zn/Zm), linspace(min(d(Di).outp.Sscale),max(d(Di).outp.Sscale),d(Di).outp.Sn/Sm));
     P_n=interp2(Zgrid_p,Sgrid,P,Zgrid_p_n,Sgrid_n);
     %P_n(P_n<1e5)=NaN;
     H.h5=surf(Zgrid_p_n,Sgrid_n,P_n,log(P_n),'linestyle','none');
     %set(get(H.h3,'Parent'),'ZScale','lin');
     %set(get(H.h3,'Parent'),'ZTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
     title('power');
     axis tight
     xlabel('z [m]');
     ylabel('s [m]');
     zlabel('P [W]');
     
     Zm=10; Lm=2;
    [Zgrid_s_n, Lamdgrid_n]=meshgrid(linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),d(Di).outp.Zn/Zm), linspace(min(d(Di).outp.Lamdscale),max(d(Di).outp.Lamdscale),numel(d(Di).outp.Lamdscale)/Lm));
    S_n=interp2(Zgrid_s,Lamdgrid,S,Zgrid_s_n,Lamdgrid_n);
    S_n(:,1)=NaN;
    figure(622);
    H.h6=surf(Zgrid_s_n,Lamdgrid_n.*1e9,S_n,log(S_n),'linestyle','none');
    %set(get(H.h4,'Parent'),'ZScale','lin');
    %set(get(H.h4,'Parent'),'ZTick',[1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12 1e13 1e14 1e15]);
    title('spectrum');
    axis tight
     xlabel('z [m]');
     ylabel('\lambda [nm]');
     zlabel('P(\lambda) [a.u.]');
    %xlim([d(Di).outp.Zscale(1) d(Di).outp.Zscale(end)]);
    
    clear Zm Lm Sm Zgrid_p Zgrid_s Sgrid Lamdgrid Zgrid_p_n Zgrid_s_n Sgrid_n Lamdgrid_n P P_n S S_n

%% Ebeam focusing plot
break
%     Xmin=min(d(Di).outp.x.mean_S-2*d(Di).outp.xrms.mean_S);
%     Xmax=max(d(Di).outp.x.mean_S+2*d(Di).outp.xrms.mean_S);
%     Ymin=min(d(Di).outp.y.mean_S-2*d(Di).outp.yrms.mean_S);
%     Ymax=max(d(Di).outp.y.mean_S+2*d(Di).outp.yrms.mean_S);
%    [X,Y]=meshgrid(linspace(Xmin,Xmax,50),linspace(Ymin,Ymax,50));
Zn=linspace(min(d(Di).outp.Zscale),max(d(Di).outp.Zscale),30);

x_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.x.mean_S,Zn);
y_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.y.mean_S,Zn);
xrms_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.xrms.mean_S,Zn);
yrms_mean_S=interp1(d(Di).outp.Zscale,d(Di).outp.yrms.mean_S,Zn);

    figure(710);
    for i=1:100
        t=2*pi/100*i;
    x=x_mean_S+1*xrms_mean_S*cos(t);
    y=y_mean_S+1*yrms_mean_S*sin(t);
    plot3(x,y,Zn,'linewidth',1.5);
    hold all
    end
    hold off
    grid on
    axis tight
    
    clear x y x_mean_S y_mean_S xrms_mean_S yrms_mean_S Zn i t
% [haxes,hline(1),hline(2)] = plotyy(outp.Zscale,outp.power.mean_S,outp.Zscale,outp.power.mean_S,'semilogy','plot');
% %set(haxes,{'ycolor'},{'b';'r'})
% ylabel(haxes(1),'P [W] (log)');
% ylabel(haxes(2),'P [W] (lin)');
% xlabel(haxes(2),'z [m]');
% set(hline(1),'LineWidth',2);
% set(hline(2),'LineWidth',1,'color',[0 0.5 0]);
% for i=1:2
%     ax(i) = get(hline(i),'Parent');
%     set(ax(i),'XLim',[outp.Zscale(1) outp.Zscale(end)]);
% end
% clear ax

%% tmp

X=d(Di).outp.Sscale*1e6;
E=d(Di).outp.energy.v(:,Zi)+d(Di).inp.gamma0;
I=d(Di).outp.current;
dE=d(Di).outp.e_spread.v(:,Zi);
dE(E<9240)=NaN;
E(E<9240)=NaN;

figure(3456)
plot(X,E,'linewidth',1.5);
xlabel('s [um]');
ylabel('\gamma');
axis tight
xlim([min(X) max(X)]);
set(figure(3456), 'Position', [100, 100, 550, 450]);

figure(3457)
plot(X,dE,'linewidth',1.5);
xlabel('s [um]');
ylabel('\sigma_\gamma');
axis tight
xlim([min(X) max(X)]);
set(figure(3457), 'Position', [100, 100, 550, 450]);

figure(3458)
plot(X,I,'linewidth',1.5);
xlabel('s [um]');
ylabel('I [A]');
axis tight
xlim([min(X) max(X)]);
set(figure(3458), 'Position', [100, 100, 550, 450]);