ic = 1;
Tint = irf.tint('2022-02-16T10:35:17.901Z/2022-02-16T10:35:19.901Z');

%%

Tintl = Tint+[-0.03 0.03]; 

Bxyz=mms.get_data('B_dmpa_brst_l2',Tintl,ic);
ne = mms.get_data('Ne_fpi_brst_l2',Tintl,ic);
c_eval('Exyz = mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',Tintl);',ic);

Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);

%%

nf = 200;
nc = 20;

Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[21.5 24.5]*1e3,'wavelet_width',5.36*50);
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
specE=struct('t',Ewavelettimes);
specE.f=Ewavelet.f/1000;
specE.p=Ewaveletx+Ewavelety+Ewaveletz;
specE.f_label='';
specE.p_label={'log_{10} E^2 (mV^2 m^{-2} Hz^{-1})'};

specFE=struct('t',Ewavelettimes);
specFE.f=Ewavelet.f/1000;
specFE.p=(Ewaveletx+Ewavelety)./(Ewaveletx+Ewavelety+Ewaveletz);
specFE.f_label='';
specFE.p_label={'F_E'};

specFE.p(specE.p < 1e-7) = NaN;

maxpower = max(max(specE.p));

%% 
Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
epso=Units.eps0;

nemedian = median(ne.data);

fpemedian = sqrt(nemedian*1e6*e^2/Me/epso)/2/pi;

dfE = 1/median(diff(Exyzfac.time.epochUnix));
Exyzfachf = Exyzfac.filt(fpemedian/1.5,0,dfE,5);
Eenv = calcenvelope(Exyzfac,fpemedian/1.5);
[Eenvpar, Eenvperp] = calcenvelopepp(Exyzfac, fpemedian/1.5);

sgtstruct = sgtstatdata(Exyzfac,fpemedian);

FE = sum(Exyzfachf.x.data.^2+Exyzfachf.y.data.^2)./sum(Exyzfachf.x.data.^2+Exyzfachf.y.data.^2+Exyzfachf.z.data.^2);

%%
h=irf_plot(5,'newfigure');
xSize=600; ySize=700;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.90;
ywidth = 0.13;
set(h(1),'position',[0.085 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.085 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.085 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.085 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.085 0.06 0.90 0.32]);

h(1)=irf_panel('Efac');
irf_plot(h(1),Exyzfachf.z,'r');
hold(h(1),'on')
irf_plot(h(1),Exyzfachf.x,'k');
irf_plot(h(1),Exyzfachf.y,'b');
hold(h(1),'off')
irf_legend(h(1),{'E_{\perp1}', 'E_{\perp2}',' E_{||}'},[0.01 0.94],'fontsize',14)
ylabel(h(1),{'E (mV m^{-1})'},'fontsize',14,'Interpreter','tex');
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',14)
title(h(1),'MMS1')

h(2)=irf_panel('Eenv');
irf_plot(h(2),Eenv,'k');
hold(h(2),'on')
irf_plot(h(2),Eenvpar,'r');
irf_plot(h(2),Eenvperp,'b');
hold(h(1),'off')
ylabel(h(2),{'E_{env} (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(2),{'E',' E_{\perp}','E_{||}'},[0.01 0.94],'fontsize',14)
irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',14)

h(3)=irf_panel('Espec');
irf_spectrogram(h(3),specE,'log');
irf_legend(h(3),'(c)',[0.99 0.94],'color','w','fontsize',14)
clim(h(3),[log10(maxpower)-9 log10(maxpower)]);
ylabel(h(3),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
colormap(h(3),'jet');

h(4)=irf_panel('FEspec');
irf_spectrogram(h(4),specFE,'lin');
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',14)
clim(h(4),[0 1]);
ylabel(h(4),{'f (kHz)'},'fontsize',14,'Interpreter','tex');
colormap(h(4),'jet');
xtickangle(h(4),0)

irf_plot_axis_align(h(1:4));
irf_zoom(h(1:4),'x',Tint);

logEvec = linspace(min(sgtstruct.logE)-0.8,max(sgtstruct.logE)+0.5,200);

sigmaSGT = sgtstruct.sigmafitESGT;
muSGT = sgtstruct.mufitESGT;
AmpE = sgtstruct.AmpfitESGT;

sigmaSGTpar = sgtstruct.sigmafitEparSGT;
muSGTpar = sgtstruct.mufitEparSGT;
AmpEpar = sgtstruct.AmpfitEparSGT;

sigmaSGTperp = sgtstruct.sigmafitEperpSGT;
muSGTperp = sgtstruct.mufitEperpSGT;
AmpEperp = sgtstruct.AmpfitEperpSGT;

PlogESGT = AmpE*exp(-(logEvec - muSGT).^2/(2*sigmaSGT.^2));
PlogEparSGT = AmpEpar*exp(-(logEvec - muSGTpar).^2/(2*sigmaSGTpar.^2));
PlogEperpSGT = AmpEperp*exp(-(logEvec - muSGTperp).^2/(2*sigmaSGTperp.^2));

sigmaSGTNL = sgtstruct.sigmafitESGTNL;
muSGTNL = sgtstruct.mufitESGTNL;
log10Ec = sgtstruct.logEfitESGTNL;
AmpENL = sgtstruct.AmpfitESGTNL;

sigmaSGTNLpar = sgtstruct.sigmafitEparSGTNL;
muSGTNLpar = sgtstruct.mufitEparSGTNL;
log10Ecpar = sgtstruct.logEfitEparSGTNL;
AmpENLpar = sgtstruct.AmpfitEparSGTNL;

sigmaSGTNLperp = sgtstruct.sigmafitEperpSGTNL;
muSGTNLperp = sgtstruct.mufitEperpSGTNL;
log10Ecperp = sgtstruct.logEfitEperpSGTNL;
AmpENLperp = sgtstruct.AmpfitEperpSGTNL;

PlogESGTNL = AmpENL*(exp(-(logEvec - muSGTNL).^2/(2*sigmaSGTNL.^2)) - exp(-(2*log10Ec -logEvec - muSGTNL).^2/(2*sigmaSGTNL.^2)));
PlogESGTNL(PlogESGTNL < 0) = NaN;

PlogEparSGTNL = AmpENLpar*(exp(-(logEvec - muSGTNLpar).^2/(2*sigmaSGTNLpar.^2)) - exp(-(2*log10Ecpar -logEvec - muSGTNLpar).^2/(2*sigmaSGTNLpar.^2)));
PlogEparSGTNL(PlogEparSGTNL < 0) = NaN;

PlogEperpSGTNL = AmpENLperp*(exp(-(logEvec - muSGTNLperp).^2/(2*sigmaSGTNLperp.^2)) - exp(-(2*log10Ecperp -logEvec - muSGTNLperp).^2/(2*sigmaSGTNLperp.^2)));
PlogEperpSGTNL(PlogEperpSGTNL < 0) = NaN;

beta1 = sgtstruct.mu3^2/sgtstruct.mu2^3;
beta2 = sgtstruct.mu4/sgtstruct.mu2^2;
beta1par = sgtstruct.mu3par^2/sgtstruct.mu2par^3;
beta2par = sgtstruct.mu4par/sgtstruct.mu2par^2;
beta1perp = sgtstruct.mu3perp^2/sgtstruct.mu2perp^3;
beta2perp = sgtstruct.mu4perp/sgtstruct.mu2perp^2;


h(5)=irf_panel('PlogElog1');
plot(h(5),sgtstruct.logE,sgtstruct.PlogE,'kx');
hold(h(5),'on')
plot(h(5),logEvec,PlogESGT,'k','linewidth',1.5);
plot(h(5),logEvec,PlogESGTNL,'k--','linewidth',1.5);
plot(h(5),sgtstruct.logEperp,sgtstruct.PlogEperp,'bx');
plot(h(5),logEvec,PlogEperpSGT,'b','linewidth',1.5);
plot(h(5),logEvec,PlogEperpSGTNL,'b--','linewidth',1.5);
plot(h(5),sgtstruct.logEpar,sgtstruct.PlogEpar,'rx');
plot(h(5),logEvec,PlogEparSGT,'r','linewidth',1.5);
plot(h(5),logEvec,PlogEparSGTNL,'r--','linewidth',1.5);
hold(h(5),'off')
ylabel(h(5),'P(log E)','fontsize',14);
xlabel(h(5),'log E','fontsize',14);
set(h(5),'yscale','log')
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',14)
axis(h(5),[-1 2 1e-3 2])
irf_legend(h(5),['\beta_1 = ' num2str(round(beta1,3))],[0.7 0.05],'color','k','fontsize',14)
irf_legend(h(5),['\beta_2 = ' num2str(round(beta2,2))],[0.7 0.16],'color','k','fontsize',14)
irf_legend(h(5),['\sigma = ' num2str(round(sqrt(sgtstruct.mu2),2))],[0.7 0.27],'color','k','fontsize',14)
irf_legend(h(5),['\mu = ' num2str(round(sgtstruct.mu1,2))],[0.7 0.38],'color','k','fontsize',14)

irf_legend(h(5),['\beta_1 = ' num2str(round(beta1par,3))],[0.40 0.05],'color','r','fontsize',14)
irf_legend(h(5),['\beta_2 = ' num2str(round(beta2par,2))],[0.40 0.16],'color','r','fontsize',14)
irf_legend(h(5),['\sigma = ' num2str(round(sqrt(sgtstruct.mu2par),2))],[0.40 0.27],'color','r','fontsize',14)
irf_legend(h(5),['\mu = ' num2str(round(sgtstruct.mu1par,2))],[0.40 0.38],'color','r','fontsize',14)

irf_legend(h(5),['\beta_1 = ' num2str(round(beta1perp,3))],[0.02 0.65],'color','b','fontsize',14)
irf_legend(h(5),['\beta_2 = ' num2str(round(beta2perp,2))],[0.02 0.76],'color','b','fontsize',14)
irf_legend(h(5),['\sigma = ' num2str(round(sqrt(sgtstruct.mu2perp),2))],[0.02 0.87],'color','b','fontsize',14)
irf_legend(h(5),['\mu = ' num2str(round(sgtstruct.mu1perp,2))],[0.02 0.98],'color','b','fontsize',14)

