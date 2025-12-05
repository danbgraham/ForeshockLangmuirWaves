%% Load data
MMS1data = loaddata(1);
MMS2data = loaddata(2);
MMS3data = loaddata(3);
MMS4data = loaddata(4);

Wmax = [MMS1data.Wmax; MMS2data.Wmax; MMS3data.Wmax; MMS4data.Wmax];
Emax = [MMS1data.Emax; MMS2data.Emax; MMS3data.Emax; MMS4data.Emax];
Region = [MMS1data.regiontype; MMS2data.regiontype; MMS3data.regiontype; MMS4data.regiontype];
Temedian = [MMS1data.Temedian; MMS2data.Temedian; MMS3data.Temedian; MMS4data.Temedian];
nemedian = [MMS1data.nemedian; MMS2data.nemedian; MMS3data.nemedian; MMS4data.nemedian];
starttimesall = [MMS1data.starttimes; MMS2data.starttimes; MMS3data.starttimes; MMS4data.starttimes];
icall = [ones(size(MMS1data.Wmax))*1; ones(size(MMS2data.Wmax))*2; ones(size(MMS3data.Wmax))*1; ones(size(MMS4data.Wmax))*4];

%clear MMS1data MMS2data MMS3data MMS4data

idxSW = Region == 1;

WmaxSW = Wmax(idxSW);
EmaxSW = Emax(idxSW);

load('sgtdata.mat');

%% 
%Te = sgtdata.Temedian;
%ne = sgtdata.nemedian;

Chi2nl = sgtdata.Chi2ESGTNLr;
Chi2 = sgtdata.Chi2ESGTr;
logEc = sgtdata.logEfitESGTNL;


idx = Chi2nl < Chi2 & Chi2nl < 10 & logEc > 0;

logEcES = logEc(idx);
WmaxES = WmaxSW(idx);


%% Make histograms
[PWmax,xWmax] = hist_dg(log10(WmaxSW),'range',[-6.5 -1.5],'nbins',100);
[PEmax,xEmax] = hist_dg(log10(EmaxSW),'range',[log10(5) 2.5],'nbins',100);
[PEc,xEc] = hist_dg(logEcES,'range',[log10(5) 2.5],'nbins',100);
[PWc,xWc] = hist_dg(log10(WmaxES),'range',[-6.5 -1.5],'nbins',100);

%% plot Figure
h=irf_plot(2,'newfigure');
xSize=800; ySize=300;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.42;
ywidth = 0.80;
set(h(1),'position',[0.07 0.17 xwidth ywidth]);
set(h(2),'position',[0.57 0.17 xwidth ywidth]);

plot(h(1),10.^xWmax,PWmax,'LineWidth',2)
yscale(h(1),'log')
xscale(h(1),'log')
xticks(h(1),[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])
hold(h(1),'on')
patch(h(1),[2e-2 1 1 2e-2], [1 1 2e3 2e3],[1 0.5 0.5],'edgecolor','none','FaceAlpha',.3)
patch(h(1),[1e-5 1e-4 1e-4 1e-5], [1 1 2e3 2e3],[0.5 0.5 1],'edgecolor','none','FaceAlpha',.3)
plot(h(1),10.^xWmax,PWmax,'LineWidth',2,'Color','k')
plot(h(1),10.^xWc,PWc,'LineWidth',2,'Color','r')
hold(h(1),'off')
xlabel(h(1),'W_{max}')
ylabel(h(1),'P(log W_{max})')
axis(h(1),[5e-7 5e-2 1 2e3])
irf_legend(h(1),'W_{max}',[0.8 0.99],'color','k','fontsize',14)
irf_legend(h(1),'W_{c}',[0.76 0.90],'color','r','fontsize',14)
irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',14)

plot(h(2),10.^xEmax,PEmax,'LineWidth',2)
yscale(h(2),'log')
xscale(h(2),'log')
hold(h(2),'on')
patch(h(2),[150 400 400 150], [1 1 2e3 2e3],[1 0.5 0.5],'edgecolor','none','FaceAlpha',.3)
patch(h(2),[5 30 30 5], [1 1 2e3 2e3],[0.5 0.5 1],'edgecolor','none','FaceAlpha',.3)
plot(h(2),10.^xEmax,PEmax,'LineWidth',2,'Color','k')
plot(h(2),10.^xEc,PEc,'r','LineWidth',2)
hold(h(2),'off')
irf_legend(h(2),'E_{max}',[0.8 0.99],'color','k','fontsize',14)
irf_legend(h(2),'E_{c}',[0.76 0.90],'color','r','fontsize',14)
axis(h(2),[5 300 1 2e3])
xlabel(h(2),'E_{max} (mV m^{-1})')
ylabel(h(2),'P(log E_{max})')
irf_legend(h(2),'(b)',[0.99 0.99],'color','k','fontsize',14)