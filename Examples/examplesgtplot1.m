%%
Tint1 = irf.tint('2022-02-16T10:35:46.182Z/2022-02-16T10:35:48.182Z');
Tint2 = irf.tint('2019-02-21T02:53:47.300Z/2019-02-21T02:53:49.300Z');

%% Load data
ic1 = 1;
ic2 = 1;

Tint1l = Tint1+[-0.03 0.03]; 
Tint2l = Tint2+[-0.03 0.03]; 

Bxyz1=mms.get_data('B_dmpa_brst_l2',Tint1l,ic1);
ne1 = mms.get_data('Ne_fpi_brst_l2',Tint1l,ic1);
c_eval('Exyz1 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint1l);',ic1);

Bxyz2=mms.get_data('B_dmpa_brst_l2',Tint2l,ic2);
ne2 = mms.get_data('Ne_fpi_brst_l2',Tint2l,ic2);
c_eval('Exyz2 = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tint2l);',ic2);

%%

Exyzfac1 = irf_convert_fac(Exyz1,Bxyz1,[1 0 0]);
Exyzfac2 = irf_convert_fac(Exyz2,Bxyz2,[1 0 0]);

%%

nf = 200;
nc = 20;

Ewavelet = irf_wavelet(Exyzfac1,'nf',nf,'f',[21.5 23.8]*1e3,'wavelet_width',5.36*50);
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE1=struct('t',Ewavelettimes);
specE1.f=Ewavelet.f/1000;
specE1.p=Ewaveletx+Ewavelety+Ewaveletz;
specE1.f_label='';
specE1.p_label={'log_{10} E^2 (mV^2 m^{-2} Hz^{-1})'};

Ewavelet = irf_wavelet(Exyzfac2,'nf',nf,'f',[23.5 26]*1e3,'wavelet_width',5.36*50);
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE2=struct('t',Ewavelettimes);
specE2.f=Ewavelet.f/1000;
specE2.p=Ewaveletx+Ewavelety+Ewaveletz;
specE2.f_label='';
specE2.p_label={'log_{10} E^2 (mV^2 m^{-2} Hz^{-1})'};

maxpower1 = max(max(specE1.p));
maxpower2 = max(max(specE2.p));

%%

Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
epso=Units.eps0;

nemedian1 = median(ne1.data);
nemedian2 = median(ne2.data);

fpemedian1 = sqrt(nemedian1*1e6*e^2/Me/epso)/2/pi;
fpemedian2 = sqrt(nemedian2*1e6*e^2/Me/epso)/2/pi;

dfE1 = 1/median(diff(Exyzfac1.time.epochUnix));
Exyzfachf1 = Exyzfac1.filt(fpemedian1/1.5,0,dfE1,5);
Eenv1 = calcenvelope(Exyzfac1,fpemedian1/1.5);

dfE2 = 1/median(diff(Exyzfac2.time.epochUnix));
Exyzfachf2 = Exyzfac2.filt(fpemedian2/1.5,0,dfE1,5);
Eenv2 = calcenvelope(Exyzfac2,fpemedian1/1.5);

sgtstruct1 = sgtstatdata(Exyzfac1,fpemedian1);
sgtstruct2 = sgtstatdata(Exyzfac2,fpemedian2);

%%

h=irf_plot(8,'newfigure');
xSize=1000; ySize=600;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.43;
ywidth = 0.17;
set(h(1),'position',[0.055 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.055 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.055 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.555 0.97-ywidth xwidth ywidth]);
set(h(5),'position',[0.555 0.97-2*ywidth xwidth ywidth]);
set(h(6),'position',[0.555 0.97-3*ywidth xwidth ywidth]);

set(h(7),'position',[0.05 0.07 0.43 0.31]);
set(h(8),'position',[0.555 0.07 0.43 0.31]);


h(1)=irf_panel('Efac1');
irf_plot(h(1),Exyzfachf1.z,'r');
hold(h(1),'on')
irf_plot(h(1),Exyzfachf1.x,'k');
irf_plot(h(1),Exyzfachf1.y,'b');
hold(h(1),'off')
irf_legend(h(1),{'E_{\perp1}','E_{\perp2}','E_{||}'},[0.01 0.94],'fontsize',14)
ylabel(h(1),{'E (mV m^{-1})'},'fontsize',14,'Interpreter','tex');
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',14)
title(h(1),'MMS1')

h(2)=irf_panel('Eenv1');
irf_plot(h(2),Eenv1);
ylabel(h(2),{'E_{env} (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',14)

h(3)=irf_panel('Espec1');
irf_spectrogram(h(3),specE1,'log');
irf_legend(h(3),'(c)',[0.99 0.94],'color','w','fontsize',14)
clim(h(3),[log10(maxpower1)-9 log10(maxpower1)]);
ylabel(h(3),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
colormap(h(3),'jet');
xtickangle(h(3),0)

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint1);
set(h(1:3),'fontsize',14);


h(4)=irf_panel('Efac2');
irf_plot(h(4),Exyzfachf2.z,'r');
hold(h(4),'on')
irf_plot(h(4),Exyzfachf2.x,'k');
irf_plot(h(4),Exyzfachf2.y,'b');
hold(h(4),'off')
ylabel(h(4),{'E (mV m^{-1})'},'fontsize',14,'Interpreter','tex');
irf_legend(h(4),'(e)',[0.99 0.94],'color','k','fontsize',14)
title(h(4),'MMS1')

h(5)=irf_panel('Eenv2');
irf_plot(h(5),Eenv2);
ylabel(h(5),{'E_{env} (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(5),'(f)',[0.99 0.94],'color','k','fontsize',14)

h(6)=irf_panel('Espec2');
irf_spectrogram(h(6),specE2,'log');
irf_legend(h(6),'(g)',[0.99 0.94],'color','w','fontsize',14)
clim(h(6),[log10(maxpower1)-9 log10(maxpower1)]);
ylabel(h(6),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
colormap(h(6),'jet');
xtickangle(h(6),0)

irf_plot_axis_align(h(4:6));
irf_zoom(h(4:6),'x',Tint2);
set(h(4:6),'fontsize',14);


logEvec1 = linspace(min(sgtstruct1.logE)-0.5,max(sgtstruct1.logE)+0.5,200);

sigmaSGT1 = sgtstruct1.sigmafitESGT;
muSGT1 = sgtstruct1.mufitESGT;
AmpE1 = sgtstruct1.AmpfitESGT;

sigmaSGTNL1 = sgtstruct1.sigmafitESGTNL;
muSGTNL1 = sgtstruct1.mufitESGTNL;
log10Ec1 = sgtstruct1.logEfitESGTNL;
AmpENL1 = sgtstruct1.AmpfitESGTNL;

PlogESGT1 = AmpE1*exp(-(logEvec1 - muSGT1).^2/(2*sigmaSGT1.^2));
PlogESGTNL1 = AmpENL1*(exp(-(logEvec1 - muSGTNL1).^2/(2*sigmaSGTNL1.^2)) - exp(-(2*log10Ec1 -logEvec1 - muSGTNL1).^2/(2*sigmaSGTNL1.^2)));
PlogESGTNL1(PlogESGTNL1 < 0) = NaN;

beta11 = sgtstruct1.mu3^2/sgtstruct1.mu2^3;
beta21 = sgtstruct1.mu4/sgtstruct1.mu2^2;

h(7)=irf_panel('PlogElog1');
plot(h(7),sgtstruct1.logE,sgtstruct1.PlogE,'kx');
hold(h(7),'on')
plot(h(7),logEvec1,PlogESGTNL1,'b','linewidth',1.5);
plot(h(7),logEvec1,PlogESGT1,'r','linewidth',1.5);
plot(h(7),logEvec1,PlogESGTNL1,'b--','linewidth',1.5);
hold(h(7),'off')
ylabel(h(7),'P(log E)','fontsize',14);
xlabel(h(7),'log E','fontsize',14);
set(h(7),'yscale','log')
irf_legend(h(7),'(d)',[0.99 0.98],'color','k','fontsize',14)
axis(h(7),[min(sgtstruct1.logE)-0.2 max(sgtstruct1.logE)+0.2 min(sgtstruct1.PlogE)*0.5 max(sgtstruct1.PlogE)*1.2])
irf_legend(h(7),['\beta_1 = ' num2str(round(beta11,3))],[0.5 0.05],'color','k','fontsize',14)
irf_legend(h(7),['\beta_2 = ' num2str(round(beta21,2))],[0.5 0.16],'color','k','fontsize',14)
irf_legend(h(7),['\sigma = ' num2str(round(sqrt(sgtstruct1.mu2),2))],[0.5 0.27],'color','k','fontsize',14)
irf_legend(h(7),['\mu = ' num2str(round(sgtstruct1.mu1,2))],[0.5 0.38],'color','k','fontsize',14)
irf_legend(h(7),['E_c = ' num2str(round(10^log10Ec1,1)) ' mV m^{-1}'],[0.6 0.6],'color','k','fontsize',14)


logEvec2 = linspace(min(sgtstruct2.logE)-0.5,max(sgtstruct2.logE)+0.5,200);

sigmaSGT2 = sgtstruct2.sigmafitESGT;
muSGT2 = sgtstruct2.mufitESGT;
AmpE2 = sgtstruct2.AmpfitESGT;

sigmaSGTNL2 = sgtstruct2.sigmafitESGTNL;
muSGTNL2 = sgtstruct2.mufitESGTNL;
log10Ec2 = sgtstruct2.logEfitESGTNL;
AmpENL2 = sgtstruct2.AmpfitESGTNL;

PlogESGT2 = AmpE2*exp(-(logEvec2 - muSGT2).^2/(2*sigmaSGT2.^2));
PlogESGTNL2 = AmpENL2*(exp(-(logEvec2 - muSGTNL2).^2/(2*sigmaSGTNL2.^2)) - exp(-(2*log10Ec2 -logEvec2 - muSGTNL2).^2/(2*sigmaSGTNL2.^2)));
PlogESGTNL2(PlogESGTNL2 < 0) = NaN;

beta12 = sgtstruct2.mu3^2/sgtstruct2.mu2^3;
beta22 = sgtstruct2.mu4/sgtstruct2.mu2^2;

h(8)=irf_panel('PlogElog2');
plot(h(8),sgtstruct2.logE,sgtstruct2.PlogE,'kx');
hold(h(8),'on')
plot(h(8),logEvec2,PlogESGTNL2,'b','linewidth',1.5);
plot(h(8),logEvec2,PlogESGT2,'r','linewidth',1.5);
plot(h(8),logEvec2,PlogESGTNL2,'b--','linewidth',1.5);
hold(h(8),'off')
ylabel(h(8),'P(log E)','fontsize',14);
xlabel(h(8),'log E','fontsize',14);
set(h(8),'yscale','log')
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',14)
axis(h(8),[-0.3 2.3 min(sgtstruct2.PlogE)*0.5 max(sgtstruct2.PlogE)*1.2])
irf_legend(h(8),['\beta_1 = ' num2str(round(beta12,3))],[0.5 0.05],'color','k','fontsize',14)
irf_legend(h(8),['\beta_2 = ' num2str(round(beta22,2))],[0.5 0.16],'color','k','fontsize',14)
irf_legend(h(8),['\sigma = ' num2str(round(sqrt(sgtstruct2.mu2),2))],[0.5 0.27],'color','k','fontsize',14)
irf_legend(h(8),['\mu = ' num2str(round(sgtstruct2.mu1,2))],[0.5 0.38],'color','k','fontsize',14)
irf_legend(h(8),['E_c = ' num2str(round(10^log10Ec2,1)) ' mV m^{-1}'],[0.6 0.6],'color','k','fontsize',14)



