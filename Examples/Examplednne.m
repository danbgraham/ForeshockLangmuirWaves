%% Initialize
Tint1 = irf.tint('2017-12-02T15:52:51.779Z/2017-12-02T15:52:53.779Z'); 
ic1 = 1;

Tint2 = irf.tint('2019-12-22T19:34:55.075Z/2019-12-22T19:34:57.075Z'); 
ic2 = 1;

%% Load data

Tint1l = Tint1+[-0.03 0.03]; 
Tint2l = Tint2+[-0.03 0.03];  

Bxyz1=mms.get_data('B_dmpa_brst_l2',Tint1l,ic1);
c_eval('Exyz1 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint1l);',ic1);

Bxyz2=mms.get_data('B_dmpa_brst_l2',Tint2l,ic2);
c_eval('Exyz2 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint2l);',ic2);


%% Parameters
Units = irf_units;
qe = Units.e;
me = Units.me;
eps0 = Units.eps0;

Te1 = 8;
Te2 = 9;

fpe1 = 2.0868e+04;
fpe2 = 2.4765e+04;

ne1 = 4*pi^2*fpe1^2*me*eps0/qe^2*1e-6;
ne2 = 4*pi^2*fpe2^2*me*eps0/qe^2*1e-6;

FE1 = 0.4911;
FE2 = 0.8126;

%% Data analysis
dfE = 1/median(diff(Exyz1.time.epochUnix));

Exyzfac1 = irf_convert_fac(Exyz1,Bxyz1,[1 0 0]);
Exyzfachf1 = Exyzfac1.filt(fpe1-2e3,fpe1+2e3,dfE,5);

Exyzfac2 = irf_convert_fac(Exyz2,Bxyz2,[1 0 0]);
Exyzfachf2 = Exyzfac2.filt(fpe2-2e3,fpe2+2e3,dfE,5);

%%
nf = 200;
nc = 20;

Ewavelet1 = irf_wavelet(Exyzfac1,'nf',nf,'f',[fpe1-1e3 fpe1+1e3],'wavelet_width',5.36*60);
Ewavelet2 = irf_wavelet(Exyzfac2,'nf',nf,'f',[fpe2-1e3 fpe2+1e3],'wavelet_width',5.36*60);

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

%% Calculate envelopes
[Eenvpar1,Eenvperp11,Eenvperp21] = calcenvelopeppp(Exyzfac1,fpe1/1.5);
Eenv1 = calcenvelope(Exyzfac1,fpe1/1.5);

[Eenvpar2,Eenvperp12,Eenvperp22] = calcenvelopeppp(Exyzfac2,fpe2/1.5);
Eenv2 = calcenvelope(Exyzfac2,fpe2/1.5);

Eenvperp1s1 = smoothdata(Eenvperp11.data,'movmean',2000);
Eenvperp2s1 = smoothdata(Eenvperp21.data,'movmean',2000);
Eenvpars1 = smoothdata(Eenvpar1.data,'movmean',2000);

Eenvperp1s2 = smoothdata(Eenvperp12.data,'movmean',2000);
Eenvperp2s2 = smoothdata(Eenvperp22.data,'movmean',2000);
Eenvpars2 = smoothdata(Eenvpar2.data,'movmean',2000);


%% Calculate instantaneous frequencies

instfrequency1 = irf_instantaneousfrequency(Exyzfachf1,'numptsav',2000,'avmethod','movmean');
instfrequency2 = irf_instantaneousfrequency(Exyzfachf2,'numptsav',2000,'avmethod','movmean');

instfreq1 = (Eenvperp1s1.^2.*instfrequency1.data(:,1)+Eenvperp2s1.^2.*instfrequency1.data(:,2)+Eenvpars1.^2.*instfrequency1.data(:,3))./(Eenvperp1s1.^2+Eenvperp2s1.^2+Eenvpars1.^2);
instfreq1 = irf.ts_scalar(instfrequency1.time,instfreq1);

instfreq2 = (Eenvperp1s2.^2.*instfrequency2.data(:,1)+Eenvperp2s2.^2.*instfrequency2.data(:,2)+Eenvpars2.^2.*instfrequency2.data(:,3))./(Eenvperp1s2.^2+Eenvperp2s2.^2+Eenvpars2.^2);
instfreq2 = irf.ts_scalar(instfrequency2.time,instfreq2);

%% Density fluctuations
dnne1 = 2*(instfreq1.data-mean(instfreq1.data))/mean(instfreq1.data);
dnne1 = irf.ts_scalar(instfreq1.time,dnne1);

dnne2 = 2*(instfreq2.data-mean(instfreq2.data))/mean(instfreq2.data);
dnne2 = irf.ts_scalar(instfreq2.time,dnne2);

%% Plot Figure

h=irf_plot(6,'newfigure');
xSize=1000; ySize=500;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.40;
ywidth = 0.29;
set(h(1),'position',[0.06 0.96-ywidth xwidth ywidth]);
set(h(2),'position',[0.06 0.96-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.06 0.96-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.57 0.96-ywidth xwidth ywidth]);
set(h(5),'position',[0.57 0.96-2*ywidth xwidth ywidth]);
set(h(6),'position',[0.57 0.96-3*ywidth xwidth ywidth]);

h(1)=irf_panel('Efac1');
irf_plot(h(1),Exyzfachf1);
hold(h(1),'on')
irf_plot(h(1),Exyzfachf1.x,'k');
irf_plot(h(1),Exyzfachf1.y,'b');
hold(h(1),'off')
irf_legend(h(1),{'E_{\perp1}',' E_{\perp2}',' E_{||}'},[0.01 0.94],'fontsize',14)
ylabel(h(1),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',14)
irf_legend(h(1),['F_E = ' num2str(FE1,2)],[0.01 0.1],'fontsize',12)
c_eval('title(h(1),''MMS?'')',ic1);

h(2)=irf_panel('Espec1');
irf_spectrogram(h(2),specE1,'log');
hold(h(2),'on')
irf_plot(h(2),instfreq1/1e3,'k','linewidth',1.5)
hold(h(2),'off')
irf_legend(h(2),'(b)',[0.99 0.94],'color','w','fontsize',14)
irf_legend(h(2),'f_{L}',[0.21 0.6],'color','k','fontsize',14)
clim(h(2),[-10 log10(maxpower1)]);
ylabel(h(2),{'f (kHz)'},'fontsize',14,'Interpreter','tex');

h(3)=irf_panel('dnne1');
irf_plot(h(3),dnne1)
irf_legend(h(3),'(c)',[0.99 0.94],'color','k','fontsize',14)
ylabel(h(3),{'\delta n_e/n_e'},'fontsize',14,'Interpreter','tex');
irf_zoom(h(3),'y',[-0.025 0.025])
xtickangle(h(3),0)

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint1);

h(4)=irf_panel('Efac2');
irf_plot(h(4),Exyzfachf2);
%hold(h(4),'on')
%irf_plot(h(4),Exyzfachf2.x,'k');
%irf_plot(h(4),Exyzfachf2.y,'b');
%irf_plot(h(4),Exyzfachf2);
%hold(h(4),'off')
ylabel(h(4),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',14)
c_eval('title(h(4),''MMS?'')',ic2);
irf_legend(h(4),['F_E = ' num2str(FE2,2)],[0.01 0.1],'fontsize',12)
irf_zoom(h(4),'y',[-11 11])

h(5)=irf_panel('Espec2');
irf_spectrogram(h(5),specE2,'log');
hold(h(5),'on')
irf_plot(h(5),instfreq2/1e3,'k','linewidth',1.5)
hold(h(5),'off')
irf_legend(h(5),'(e)',[0.99 0.94],'color','w','fontsize',14)
clim(h(5),[-10 log10(maxpower2)]);
ylabel(h(5),{'f (kHz)'},'fontsize',14,'Interpreter','tex');

h(6)=irf_panel('dnne2');
irf_plot(h(6),dnne2)
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',14)
ylabel(h(6),{'\delta n_e/n_e'},'fontsize',14,'Interpreter','tex');
xtickangle(h(6),0)
irf_zoom(h(6),'y',[-1.01e-2 1.01e-2])

irf_plot_axis_align(h(4:6));
irf_zoom(h(4:6),'x',Tint2);

colormap(h(2),'jet');
colormap(h(5),'jet');

%irf_plot_axis_align(h([1:3 7:9]));
%irf_plot_axis_align(h([4:6 10:12]));

set(h(1:6),'fontsize',14);



