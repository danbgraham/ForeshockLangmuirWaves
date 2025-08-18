%% Initialize
Tint1 = irf.tint('2017-12-18T12:02:37.545Z/2017-12-18T12:02:39.545Z'); 
ic1 = 1;

Tint2 = irf.tint('2018-12-10T05:27:50.733Z/2018-12-10T05:27:52.733Z'); 
ic2 = 1;

Tint3 = irf.tint('2022-03-10T09:20:50.306Z/2022-03-10T09:20:52.306Z');
ic3 = 1;

%% Load data

Tint1l = Tint1+[-0.03 0.03]; 
Tint2l = Tint2+[-0.03 0.03]; 
Tint3l = Tint3+[-0.03 0.03]; 

Bxyz1=mms.get_data('B_dmpa_brst_l2',Tint1l,ic1);
ne1 = mms.get_data('Ne_fpi_brst_l2',Tint1l,ic1);
c_eval('Exyz1 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint1l);',ic1);

Bxyz2=mms.get_data('B_dmpa_brst_l2',Tint2l,ic2);
ne2 = mms.get_data('Ne_fpi_brst_l2',Tint2l,ic2);
c_eval('Exyz2 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint2l);',ic2);

Bxyz3=mms.get_data('B_dmpa_brst_l2',Tint3l,ic3);
ne3 = mms.get_data('Ne_fpi_brst_l2',Tint3l,ic3);
c_eval('Exyz3 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint3l);',ic3);

%% Data analysis
fmin = 1e2;
dfE = 1/median(diff(Exyz1.time.epochUnix));

Exyzfac1 = irf_convert_fac(Exyz1,Bxyz1,[1 0 0]);
Exyzfac1 = Exyzfac1.tlim(Tint1);
Exyzfachf1 = Exyzfac1.filt(fmin,0,dfE,5);

Exyzfac2 = irf_convert_fac(Exyz2,Bxyz2,[1 0 0]);
Exyzfac2 = Exyzfac2.tlim(Tint2);
Exyzfachf2 = Exyzfac2.filt(fmin,0,dfE,5);

Exyzfac3 = irf_convert_fac(Exyz3,Bxyz3,[1 0 0]);
Exyzfac3 = Exyzfac3.tlim(Tint3);
Exyzfachf3 = Exyzfac3.filt(fmin,0,dfE,5);

Units = irf_units;
qe = Units.e;
me = Units.me;
epso=Units.eps0;

fpe1 = sqrt(qe^2*ne1.data*1e6/(me*epso))/2/pi;
fpe2 = sqrt(qe^2*ne2.data*1e6/(me*epso))/2/pi;
fpe3 = sqrt(qe^2*ne3.data*1e6/(me*epso))/2/pi;

fpe1 = irf.ts_scalar(ne1.time,fpe1/1e3);
fpe2 = irf.ts_scalar(ne2.time,fpe2/1e3);
fpe3 = irf.ts_scalar(ne3.time,fpe3/1e3);

Ewavelet1 = irf_wavelet(Exyzfac1,'linear',100,'wavelet_width',5.36*30);
Ewavelet2 = irf_wavelet(Exyzfac2,'linear',100,'wavelet_width',5.36*30);
Ewavelet3 = irf_wavelet(Exyzfac3,'linear',100,'wavelet_width',5.36*30);

% Int 1
nc = 20;
nf = length(Ewavelet1.f);
idx = [nc/2:nc:length(Ewavelet1.t)-nc/2];
Ewavelettimes = Ewavelet1.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet1.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet1.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet1.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE1=struct('t',Ewavelettimes);
specE1.f=Ewavelet1.f/1000;
specE1.p=Ewaveletx+Ewavelety+Ewaveletz;
specE1.f_label='';
specE1.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

maxpower1 = max(max(Ewaveletx+Ewavelety+Ewaveletz));

specFE1=struct('t',Ewavelettimes);
specFE1.f=Ewavelet1.f/1000;
specFE1.p=(Ewaveletx+Ewavelety)./(Ewaveletx+Ewavelety+Ewaveletz);
specFE1.f_label='';
specFE1.p_label={'F_E'};
specFE1.p(specE1.p < maxpower1/1e5) = NaN;

% Int 2
idx = [nc/2:nc:length(Ewavelet2.t)-nc/2];
Ewavelettimes = Ewavelet2.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet2.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet2.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet2.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE2=struct('t',Ewavelettimes);
specE2.f=Ewavelet2.f/1000;
specE2.p=Ewaveletx+Ewavelety+Ewaveletz;
specE2.f_label='';
specE2.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

maxpower2 = max(max(Ewaveletx+Ewavelety+Ewaveletz));

specFE2=struct('t',Ewavelettimes);
specFE2.f=Ewavelet2.f/1000;
specFE2.p=(Ewaveletx+Ewavelety)./(Ewaveletx+Ewavelety+Ewaveletz);
specFE2.f_label='';
specFE2.p_label={'F_E'};
specFE2.p(specE2.p < maxpower2/1e5) = NaN;

% Int 3
idx = [nc/2:nc:length(Ewavelet3.t)-nc/2];
Ewavelettimes = Ewavelet3.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet3.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet3.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet3.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE3=struct('t',Ewavelettimes);
specE3.f=Ewavelet3.f/1000;
specE3.p=Ewaveletx+Ewavelety+Ewaveletz;
specE3.f_label='';
specE3.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

maxpower3 = max(max(Ewaveletx+Ewavelety+Ewaveletz));

specFE3=struct('t',Ewavelettimes);
specFE3.f=Ewavelet3.f/1000;
specFE3.p=(Ewaveletx+Ewavelety)./(Ewaveletx+Ewavelety+Ewaveletz);
specFE3.f_label='';
specFE3.p_label={'F_E'};
specFE3.p(specE3.p < maxpower3/1e5) = NaN;

%% Plot Figure

h=irf_plot(9,'newfigure');
xSize=1400; ySize=500;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.27;
ywidth = 0.29;
set(h(1),'position',[0.04 0.96-ywidth xwidth ywidth]);
set(h(2),'position',[0.04 0.96-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.04 0.96-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.37 0.96-ywidth xwidth ywidth]);
set(h(5),'position',[0.37 0.96-2*ywidth xwidth ywidth]);
set(h(6),'position',[0.37 0.96-3*ywidth xwidth ywidth]);
set(h(7),'position',[0.70 0.96-ywidth xwidth ywidth]);
set(h(8),'position',[0.70 0.96-2*ywidth xwidth ywidth]);
set(h(9),'position',[0.70 0.96-3*ywidth xwidth ywidth]);

h(1)=irf_panel('Efac1');
irf_plot(h(1),Exyzfachf1);
hold(h(1),'on')
irf_plot(h(1),Exyzfachf1.x,'k');
irf_plot(h(1),Exyzfachf1.y,'b');
hold(h(1),'off')
irf_legend(h(1),{'E_{\perp1}','E_{\perp2}','E_{||}'},[0.01 0.94],'fontsize',14)
ylabel(h(1),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',14)
c_eval('title(h(1),''MMS?'')',ic1);

h(2)=irf_panel('Espec1');
irf_spectrogram(h(2),specE1,'log');
hold(h(2),'on')
irf_plot(h(2),fpe1,'k')
hold(h(2),'off')
irf_legend(h(2),'(b)',[0.99 0.94],'color','w','fontsize',14)
irf_legend(h(2),'f_{pe}',[0.2 0.6],'color','k','fontsize',14)
clim(h(2),[log10(maxpower1)-9 log10(maxpower1)]);
ylabel(h(2),{'f (kHz)'},'fontsize',14,'Interpreter','tex');

h(3)=irf_panel('FEspec1');
irf_spectrogram(h(3),specFE1,'lin');
hold(h(3),'on')
irf_plot(h(3),fpe1,'k')
hold(h(3),'off')
irf_legend(h(3),'(c)',[0.99 0.94],'color','k','fontsize',14)
clim(h(3),[0 1]);
ylabel(h(3),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
xtickangle(h(3),0)

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint1);


h(4)=irf_panel('Efac2');
irf_plot(h(4),Exyzfachf2);
hold(h(4),'on')
irf_plot(h(4),Exyzfachf2.x,'k');
irf_plot(h(4),Exyzfachf2.y,'b');
hold(h(4),'off')
ylabel(h(4),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',14)
c_eval('title(h(4),''MMS?'')',ic2);

h(5)=irf_panel('Espec2');
irf_spectrogram(h(5),specE2,'log');
hold(h(5),'on')
irf_plot(h(5),fpe2,'k')
hold(h(5),'off')
irf_legend(h(5),'(e)',[0.99 0.94],'color','w','fontsize',14)
clim(h(5),[log10(maxpower2)-9 log10(maxpower2)]);
ylabel(h(5),{'f (kHz)'},'fontsize',14,'Interpreter','tex');

h(6)=irf_panel('FEspec2');
irf_spectrogram(h(6),specFE2,'lin');
hold(h(6),'on')
irf_plot(h(6),fpe2,'k')
hold(h(6),'off')
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',14)
clim(h(6),[0 1]);
ylabel(h(6),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
xtickangle(h(6),0)

irf_plot_axis_align(h(4:6));
irf_zoom(h(4:6),'x',Tint2);

h(7)=irf_panel('Efac3');
irf_plot(h(7),Exyzfachf3);
hold(h(7),'on')
irf_plot(h(7),Exyzfachf3.x,'k');
irf_plot(h(7),Exyzfachf3.y,'b');
hold(h(7),'off')
ylabel(h(7),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(7),'(g)',[0.99 0.94],'color','k','fontsize',14)
c_eval('title(h(7),''MMS?'')',ic3);

h(8)=irf_panel('Espec3');
irf_spectrogram(h(8),specE3,'log');
hold(h(8),'on')
irf_plot(h(8),fpe3,'k')
hold(h(8),'off')
irf_legend(h(8),'(h)',[0.99 0.94],'color','w','fontsize',14)
clim(h(8),[log10(maxpower3)-9 log10(maxpower3)]);
ylabel(h(8),{'f (kHz)'},'fontsize',14,'Interpreter','tex');

h(9)=irf_panel('FEspec3');
irf_spectrogram(h(9),specFE3,'lin');
hold(h(9),'on')
irf_plot(h(9),fpe3,'k')
hold(h(9),'off')
irf_legend(h(9),'(i)',[0.99 0.94],'color','k','fontsize',14)
clim(h(9),[0 1]);
ylabel(h(9),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
xtickangle(h(9),0)

irf_plot_axis_align(h(7:9));
irf_zoom(h(7:9),'x',Tint3);

colormap(h(2),'jet');
colormap(h(3),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(8),'jet');
colormap(h(9),'jet');

%irf_plot_axis_align(h([1:3 7:9]));
%irf_plot_axis_align(h([4:6 10:12]));

set(h(1:9),'fontsize',14);



