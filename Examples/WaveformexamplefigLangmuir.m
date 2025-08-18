%% Initialize
Tint1 = irf.tint('2017-12-07T06:23:47.694Z/2017-12-07T06:23:49.694Z'); 
ic1 = 2;

Tint2 = irf.tint('2019-01-23T14:48:42.989Z/2019-01-23T14:48:44.989Z'); 
ic2 = 1;

Tint3 = irf.tint('2017-12-18T14:53:12.148Z/2017-12-18T14:53:14.148Z');
ic3 = 3;

Tint4 = irf.tint('2017-10-26T18:18:58.523Z/2017-10-26T18:19:00.523Z');
ic4 = 1;

%% Load data

Tint1l = Tint1+[-0.03 0.03]; 
Tint2l = Tint2+[-0.03 0.03]; 
Tint3l = Tint3+[-0.03 0.03]; 
Tint4l = Tint4+[-0.03 0.03]; 

Bxyz1=mms.get_data('B_dmpa_brst_l2',Tint1l,ic1);
c_eval('Exyz1 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint1l);',ic1);

Bxyz2=mms.get_data('B_dmpa_brst_l2',Tint2l,ic2);
c_eval('Exyz2 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint2l);',ic2);

Bxyz3=mms.get_data('B_dmpa_brst_l2',Tint3l,ic3);
c_eval('Exyz3 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint3l);',ic3);

Bxyz4=mms.get_data('B_dmpa_brst_l2',Tint4l,ic4);
c_eval('Exyz4 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint4l);',ic4);

%% Data analysis
fmin = 1e4;
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

Exyzfac4 = irf_convert_fac(Exyz4,Bxyz4,[1 0 0]);
Exyzfac4 = Exyzfac4.tlim(Tint4);
Exyzfachf4 = Exyzfac4.filt(fmin,0,dfE,5);


nf = 200;
nc = 20;


Ewavelet1 = irf_wavelet(Exyzfac1,'nf',nf,'f',[1.66e4 1.85e4],'wavelet_width',5.36*60);
Ewavelet2 = irf_wavelet(Exyzfac2,'nf',nf,'f',[2.6e4 2.85e4],'wavelet_width',5.36*60);
Ewavelet3 = irf_wavelet(Exyzfac3,'nf',nf,'f',[1.25e4 1.5e4],'wavelet_width',5.36*60);
Ewavelet4 = irf_wavelet(Exyzfac4,'nf',nf,'f',[1.6e4 1.85e4],'wavelet_width',5.36*60);

% Int 1
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
specFE1.p(specE1.p < maxpower1/1e6) = NaN;

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
specFE2.p(specE2.p < maxpower2/1e6) = NaN;

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
specFE3.p(specE3.p < maxpower3/1e6) = NaN;

% Int 4
idx = [nc/2:nc:length(Ewavelet4.t)-nc/2];
Ewavelettimes = Ewavelet4.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet4.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet4.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet4.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE4=struct('t',Ewavelettimes);
specE4.f=Ewavelet4.f/1000;
specE4.p=Ewaveletx+Ewavelety+Ewaveletz;
specE4.f_label='';
specE4.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

maxpower4 = max(max(Ewaveletx+Ewavelety+Ewaveletz));

specFE4=struct('t',Ewavelettimes);
specFE4.f=Ewavelet4.f/1000;
specFE4.p=(Ewaveletx+Ewavelety)./(Ewaveletx+Ewavelety+Ewaveletz);
specFE4.f_label='';
specFE4.p_label={'F_E'};
specFE4.p(specE4.p < maxpower4/1e6) = NaN;

%% Plot Figure

h=irf_plot(12,'newfigure');
xSize=1000; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.45;
ywidth = 0.14;
set(h(1),'position',[0.05 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.05 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.05 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.55 0.97-ywidth xwidth ywidth]);
set(h(5),'position',[0.55 0.97-2*ywidth xwidth ywidth]);
set(h(6),'position',[0.55 0.97-3*ywidth xwidth ywidth]);
set(h(7),'position',[0.05 0.89-4*ywidth xwidth ywidth]);
set(h(8),'position',[0.05 0.89-5*ywidth xwidth ywidth]);
set(h(9),'position',[0.05 0.89-6*ywidth xwidth ywidth]);
set(h(10),'position',[0.55 0.89-4*ywidth xwidth ywidth]);
set(h(11),'position',[0.55 0.89-5*ywidth xwidth ywidth]);
set(h(12),'position',[0.55 0.89-6*ywidth xwidth ywidth]);

h(1)=irf_panel('Efac1');
irf_plot(h(1),Exyzfachf1);
hold(h(1),'on')
irf_plot(h(1),Exyzfachf1.x,'k');
irf_plot(h(1),Exyzfachf1.y,'b');
hold(h(1),'off')
irf_legend(h(1),{'E_{\perp1}','E_{\perp2}','E_{||}'},[0.01 0.94],'fontsize',12)
ylabel(h(1),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',12)
c_eval('title(h(1),''MMS?'')',ic1);

h(2)=irf_panel('Espec1');
irf_spectrogram(h(2),specE1,'log');
irf_legend(h(2),'(b)',[0.99 0.94],'color','w','fontsize',12)
clim(h(2),[log10(maxpower1)-9 log10(maxpower1)]);
ylabel(h(2),{'f (kHz)'},'fontsize',12,'Interpreter','tex');

h(3)=irf_panel('FEspec1');
irf_spectrogram(h(3),specFE1,'lin');
irf_legend(h(3),'(c)',[0.99 0.94],'color','k','fontsize',12)
clim(h(3),[0 1]);
ylabel(h(3),{'f (kHz)'},'fontsize',12,'Interpreter','tex');
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
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',12)
c_eval('title(h(4),''MMS?'')',ic2);

h(5)=irf_panel('Espec2');
irf_spectrogram(h(5),specE2,'log');
irf_legend(h(5),'(e)',[0.99 0.94],'color','w','fontsize',12)
clim(h(5),[log10(maxpower2)-9 log10(maxpower2)]);
ylabel(h(5),{'f (kHz)'},'fontsize',12,'Interpreter','tex');

h(6)=irf_panel('FEspec2');
irf_spectrogram(h(6),specFE2,'lin');
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',12)
clim(h(6),[0 1]);
ylabel(h(6),{'f (kHz)'},'fontsize',12,'Interpreter','tex');
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
irf_legend(h(7),'(g)',[0.99 0.94],'color','k','fontsize',12)
c_eval('title(h(7),''MMS?'')',ic3);

h(8)=irf_panel('Espec3');
irf_spectrogram(h(8),specE3,'log');
irf_legend(h(8),'(h)',[0.99 0.94],'color','w','fontsize',12)
clim(h(8),[log10(maxpower3)-9 log10(maxpower3)]);
ylabel(h(8),{'f (kHz)'},'fontsize',12,'Interpreter','tex');

h(9)=irf_panel('FEspec3');
irf_spectrogram(h(9),specFE3,'lin');
irf_legend(h(9),'(i)',[0.99 0.94],'color','k','fontsize',12)
clim(h(9),[0 1]);
ylabel(h(9),{'f (kHz)'},'fontsize',12,'Interpreter','tex');
xtickangle(h(9),0)

irf_plot_axis_align(h(7:9));
irf_zoom(h(7:9),'x',Tint3);


h(10)=irf_panel('Efac4');
irf_plot(h(10),Exyzfachf4);
hold(h(10),'on')
irf_plot(h(10),Exyzfachf4.x,'k');
irf_plot(h(10),Exyzfachf4.y,'b');
hold(h(10),'off')
ylabel(h(10),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(10),'(j)',[0.99 0.94],'color','k','fontsize',12)
c_eval('title(h(10),''MMS?'')',ic4);

h(11)=irf_panel('Espec4');
irf_spectrogram(h(11),specE4,'log');
irf_legend(h(11),'(k)',[0.99 0.94],'color','w','fontsize',12)
clim(h(11),[log10(maxpower4)-9 log10(maxpower4)]);
ylabel(h(11),{'f (kHz)'},'fontsize',12,'Interpreter','tex');

h(12)=irf_panel('FEspec4');
irf_spectrogram(h(12),specFE4,'lin');
irf_legend(h(12),'(l)',[0.99 0.94],'color','k','fontsize',12)
clim(h(12),[0 1]);
ylabel(h(12),{'f (kHz)'},'fontsize',12,'Interpreter','tex');
xtickangle(h(12),0)

irf_plot_axis_align(h(10:12));
irf_zoom(h(10:12),'x',Tint4);

colormap(h(2),'jet');
colormap(h(3),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(8),'jet');
colormap(h(9),'jet');
colormap(h(11),'jet');
colormap(h(12),'jet');

irf_plot_axis_align(h([1:3 7:9]));
irf_plot_axis_align(h([4:6 10:12]));

set(h(1:12),'fontsize',12);

