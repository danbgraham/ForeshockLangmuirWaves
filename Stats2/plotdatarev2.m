load('data/MMS1data.mat')
load('data/MMS2data.mat')
load('data/MMS3data.mat')
load('data/MMS4data.mat')

starttimes = [MMS1data.starttimes; MMS2data.starttimes; MMS3data.starttimes; MMS4data.starttimes];
endtimes = [MMS1data.endtimes; MMS2data.endtimes; MMS3data.endtimes; MMS4data.endtimes];
F = [MMS1data.F; MMS2data.F; MMS3data.F; MMS4data.F];
R = [MMS1data.R; MMS2data.R; MMS3data.R; MMS4data.R];
FScoords = [MMS1data.FScoords; MMS2data.FScoords; MMS3data.FScoords; MMS4data.FScoords];
FStanpoint = [MMS1data.FStanpoint; MMS2data.FStanpoint; MMS3data.FStanpoint; MMS4data.FStanpoint];
Emax = [MMS1data.Emax; MMS2data.Emax; MMS3data.Emax; MMS4data.Emax];
Emaxpar = [MMS1data.Emaxpar; MMS2data.Emaxpar; MMS3data.Emaxpar; MMS4data.Emaxpar];
Emaxperp = [MMS1data.Emaxperp; MMS2data.Emaxperp; MMS3data.Emaxperp; MMS4data.Emaxperp];
Ermspar = [MMS1data.Ermspar; MMS2data.Ermspar; MMS3data.Ermspar; MMS4data.Ermspar];
Ermsperp1 = [MMS1data.Ermsperp1; MMS2data.Ermsperp1; MMS3data.Ermsperp1; MMS4data.Ermsperp1];
Ermsperp2 = [MMS1data.Ermsperp2; MMS2data.Ermsperp2; MMS3data.Ermsperp2; MMS4data.Ermsperp2];
fpe = [MMS1data.fpe; MMS2data.fpe; MMS3data.fpe; MMS4data.fpe];
nemedian = [MMS1data.nemedian; MMS2data.nemedian; MMS3data.nemedian; MMS4data.nemedian];
Temedian = [MMS1data.Temedian; MMS2data.Temedian; MMS3data.Temedian; MMS4data.Temedian];
nimedian = [MMS1data.nimedian; MMS2data.nimedian; MMS3data.nimedian; MMS4data.nimedian];
Timedian = [MMS1data.Timedian; MMS2data.Timedian; MMS3data.Timedian; MMS4data.Timedian];
Vimedian = [MMS1data.Vimedian; MMS2data.Vimedian; MMS3data.Vimedian; MMS4data.Vimedian];
BGSEmedian = [MMS1data.BGSEmedian; MMS2data.BGSEmedian; MMS3data.BGSEmedian; MMS4data.BGSEmedian];
fpeak = [MMS1data.fpeak; MMS2data.fpeak; MMS3data.fpeak; MMS4data.fpeak];
Pfpeak = [MMS1data.Pfpeak; MMS2data.Pfpeak; MMS3data.Pfpeak; MMS4data.Pfpeak];
Pfmedian = [MMS1data.Pfmedian; MMS2data.Pfmedian; MMS3data.Pfmedian; MMS4data.Pfmedian];
flagASPOC = [MMS1data.flagASPOC; MMS2data.flagASPOC; MMS3data.flagASPOC; MMS4data.flagASPOC];
flagEDI = [MMS1data.flagEDI; MMS2data.flagEDI; MMS3data.flagEDI; MMS4data.flagEDI];

%Tintcutoff = irf_time([2019 03 19 00 00 00.0]);
timelengths = endtimes-starttimes;
idx = timelengths > 1.1;
starttimes = starttimes(idx);
endtimes = endtimes(idx);
F = F(idx);
R = R(idx,:);
FScoords = FScoords(idx,:);
FStanpoint = FStanpoint(idx,:);
Emax = Emax(idx);
Emaxpar = Emaxpar(idx);
Emaxperp = Emaxperp(idx);
Ermspar = Ermspar(idx);
Ermsperp1 = Ermsperp1(idx);
Ermsperp2 = Ermsperp2(idx);
fpe = fpe(idx);
nemedian = nemedian(idx);
Temedian = Temedian(idx);
nimedian = nimedian(idx);
Timedian = Timedian(idx);
Vimedian = Vimedian(idx,:);
BGSEmedian = BGSEmedian(idx,:);
fpeak = fpeak(idx);
Pfpeak = Pfpeak(idx);
Pfmedian = Pfmedian(idx);
flagASPOC = flagASPOC(idx);
flagEDI = flagEDI(idx);

Ermstot = sqrt(Ermspar.^2 + Ermsperp1.^2+ Ermsperp2.^2);
Ermsperp = sqrt(Ermsperp1.^2+ Ermsperp2.^2);

%% Plot probably of large electric field in the electron foreshock

idx2 = abs(FScoords(:,1)) < 50 & abs(FScoords(:,2)) < 50;
idx2E = idx2 & Emax > 0.1;
idx2E2 = idx2 & Emax > 1;

Xedges = [-50:1:50];
Yedges = [-50:1:50];
[N,Xedges,Yedges] = histcounts2(FScoords(idx2,1),FScoords(idx2,2),Xedges,Yedges);
[NE,Xedges,Yedges] = histcounts2(FScoords(idx2E,1),FScoords(idx2E,2),Xedges,Yedges);
[NE2,Xedges,Yedges] = histcounts2(FScoords(idx2E2,1),FScoords(idx2E2,2),Xedges,Yedges);
Xpos = [-49.5:1:49.5];
Ypos = [-49.5:1:49.5];

percentW = NE'./N'*1e2;
percentW2 = NE2'./N'*1e2;
percentW(N' < 5) = NaN;
percentW2(N' < 5) = NaN;

N(N < 1) = NaN;


Dfvec = [-14.75:0.5:49.75];
Emaxmean = zeros(size(Dfvec));
Emaxmedian = zeros(size(Dfvec));
Emaxparmean = zeros(size(Dfvec));
Emaxparmedian = zeros(size(Dfvec));
Emaxperpmean = zeros(size(Dfvec));
Emaxperpmedian = zeros(size(Dfvec));

idx2 = abs(FScoords(:,1)) < 50 & abs(FScoords(:,2)) < 50;

for ii = 1:length(Dfvec)
  idxtemp = FScoords(:,1) > Dfvec(ii)-0.25 & FScoords(:,1) < Dfvec(ii)+0.25;
  Emaxmean(ii) = mean(Emax(idxtemp));
  Emaxmedian(ii) = median(Emax(idxtemp));
  Emaxparmean(ii) = mean(Emaxpar(idxtemp));
  Emaxparmedian(ii) = median(Emaxpar(idxtemp));
  Emaxperpmean(ii) = mean(Emaxperp(idxtemp));
  Emaxperpmedian(ii) = median(Emaxperp(idxtemp));
end

%% Figure

fn=figure;
set(fn,'Position',[10 10 700 900])
h(1)=axes('position',[0.07 0.55 0.25 0.41]); % [x y dx dy]
h(2)=axes('position',[0.40 0.55 0.25 0.41]); % [x y dx dy]
h(3)=axes('position',[0.73 0.55 0.25 0.41]); % [x y dx dy]
h(4)=axes('position',[0.11 0.35 0.87 0.14]); % [x y dx dy]
h(5)=axes('position',[0.11 0.205 0.87 0.14]); % [x y dx dy]
h(6)=axes('position',[0.11 0.06 0.87 0.14]); % [x y dx dy]
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 
set(fn,'defaultAxesFontSize',14)

pcolor(h(1),Xpos,Ypos,log10(N'))
shading(h(1),'flat');
c=colorbar('peer',h(1),'ver');
colormap(h(1),'jet');
c.Label.String = 'log_{10}(counts)';
hold(h(1),'on');
plot(h(1),[0 0],[-50 50],'m--','linewidth',2)
hold(h(1),'off');
xlabel(h(1),'D_f (R_E)','fontsize',14)
ylabel(h(1),'R (R_E)','fontsize',14)
title(h(1),'Counts','fontsize',14);
set(gcf,'color','w')
irf_legend(h(1),'(a)',[0.04 0.98],'fontsize',14)
axis(h(1),[-20 50 -50 50])

pcolor(h(2),Xpos,Ypos,percentW)
shading(h(2),'flat');
c=colorbar('peer',h(2),'ver');
colormap(h(2),'jet');
c.Label.String = '%';
hold(h(2),'on');
plot(h(2),[0 0],[-50 50],'m--','linewidth',2)
hold(h(2),'off');
xlabel(h(2),'D_f (R_E)','fontsize',14)
ylabel(h(2),'R (R_E)','fontsize',14)
title(h(2),'E_{max} > 0.1 mV m^{-1} Percentage','fontsize',14);
set(gcf,'color','w')
irf_legend(h(2),'(b)',[0.04 0.98],'fontsize',14)
axis(h(2),[-20 50 -50 50])

pcolor(h(3),Xpos,Ypos,percentW2)
shading(h(3),'flat');
c=colorbar('peer',h(3),'ver');
colormap(h(3),'jet');
c.Label.String = '%';
hold(h(3),'on');
plot(h(3),[0 0],[-50 50],'m--','linewidth',2)
hold(h(3),'off');
xlabel(h(3),'D_f (R_E)','fontsize',14)
ylabel(h(3),'R (R_E)','fontsize',14)
title(h(3),'E_{max} > 1 mV m^{-1} Percentage','fontsize',14);
set(gcf,'color','w')
irf_legend(h(3),'(c)',[0.04 0.98],'fontsize',14)
axis(h(3),[-20 50 -50 50])


plot(h(4),FScoords(idx2,1),Emax(idx2),'k.')
hold(h(4),'on')
plot(h(4),Dfvec,Emaxmean,'r','linewidth',2)
plot(h(4),Dfvec,Emaxmedian,'g','linewidth',2)
plot(h(4),[0 0],[1e-2 1e3],'m--','linewidth',2)
hold(h(4),'off')
%xlabel(h(4),'D_f (R_E)','interpreter','tex','fontsize',14)
ylabel(h(4),{'E_{max}','(mV m^{-1})'},'interpreter','tex','fontsize',14)
set(h(4),'fontsize',14)
irf_legend(h(4),'mean',[0.90, 0.98],'fontsize',14,'color','r');
irf_legend(h(4),'median',[0.90, 0.85],'fontsize',14,'color','g');
irf_legend(h(4),'(d)',[0.98, 0.98],'fontsize',14,'color','k');
set(h(4),'yscale','log')
axis(h(4),[-15 50 5e-3 3e2])
set(h(4), 'XTickLabel', [])
yticks(h(4),[1e-2 1e0 1e2])

plot(h(5),FScoords(idx2,1),Emaxpar(idx2),'k.')
hold(h(5),'on')
plot(h(5),Dfvec,Emaxparmean,'r','linewidth',2)
plot(h(5),Dfvec,Emaxparmedian,'g','linewidth',2)
plot(h(5),[0 0],[1e-2 1e3],'m--','linewidth',2)
hold(h(5),'off')
%xlabel(h(5),'D_f (R_E)','interpreter','tex','fontsize',14)
ylabel(h(5),{'E_{||,max}','(mV m^{-1})'},'interpreter','tex','fontsize',14)
set(h(5),'fontsize',14)
irf_legend(h(5),'(e)',[0.98, 0.98],'fontsize',14,'color','k');
set(h(5),'yscale','log')
axis(h(5),[-15 50 5e-3 3e2])
set(h(5), 'XTickLabel', [])
yticks(h(5),[1e-2 1e0 1e2])

plot(h(6),FScoords(idx2,1),Emaxperp(idx2),'k.')
hold(h(6),'on')
plot(h(6),Dfvec,Emaxperpmean,'r','linewidth',2)
plot(h(6),Dfvec,Emaxperpmedian,'g','linewidth',2)
plot(h(6),[0 0],[1e-2 1e3],'m--','linewidth',2)
hold(h(6),'off')
xlabel(h(6),'D_f (R_E)','interpreter','tex','fontsize',14)
ylabel(h(6),{'E_{\perp,max}','(mV m^{-1})'},'interpreter','tex','fontsize',14)
set(h(6),'fontsize',14)
irf_legend(h(6),'(f)',[0.98, 0.98],'fontsize',14,'color','k');
set(h(6),'yscale','log')
axis(h(6),[-15 50 5e-3 3e2])
yticks(h(6),[1e-2 1e0 1e2])

set(gcf,'color','w')

%% Histograms of all electric fields

Dfall = FScoords(:,1);
Emaxf = Emax(Emax > 1);
Emaxperpf = Emaxperp(Emaxperp > 1);
Emaxparf = Emaxpar(Emaxpar > 1);

Ermstotf = Ermstot(Emax > 1);
Ermsperpf = Ermsperp(Emaxperp > 1);
Ermsparf = Ermspar(Emaxpar > 1);

[Pl,El] = hist_dg(log10(Emax),'range',[-4 2.5],'nbins',200);
[Plperp,Elperp] = hist_dg(log10(Emaxperp),'range',[-4 2.5],'nbins',200);
[Plpar,Elpar] = hist_dg(log10(Emaxpar),'range',[-4 2.5],'nbins',200);

[Plr,Elr] = hist_dg(log10(Ermstot),'range',[-4 2.5],'nbins',200);
[Plperpr,Elperpr] = hist_dg(log10(Ermsperp),'range',[-4 2.5],'nbins',200);
[Plparr,Elparr] = hist_dg(log10(Ermspar),'range',[-4 2.5],'nbins',200);

idxtemp = El < 0.3 & El > 0.3;
grad = -0.7;
ycross = 3.65;
xx = 0:0.01:2.5;

grad2 = -1.0;
ycross2 = 3.6;

[Prat,Erat] = hist_dg(log10(Emaxf./Ermstotf),'range',[0 3],'nbins',100);
[Pperprat,Eperprat] = hist_dg(log10(Emaxperpf./Ermsperpf),'range',[0 3],'nbins',100);
[Pparrat,Eparrat] = hist_dg(log10(Emaxparf./Ermsparf),'range',[0 3],'nbins',100);

%% Distribution of electric fields for ranges of Df
%Total electric field
idxDf1 = FScoords(:,1) < 2 & FScoords(:,1) > 0 & Emax > 0.1;
idxDf2 = FScoords(:,1) < 4 & FScoords(:,1) > 2 & Emax > 0.1;
idxDf3 = FScoords(:,1) < 6 & FScoords(:,1) > 4 & Emax > 0.1;
idxDf4 = FScoords(:,1) < 8 & FScoords(:,1) > 6 & Emax > 0.1;
idxDf5 = FScoords(:,1) < 10 & FScoords(:,1) > 8 & Emax > 0.1;
idxDf6 = FScoords(:,1) < 12 & FScoords(:,1) > 10 & Emax > 0.1;
idxDf7 = FScoords(:,1) < 14 & FScoords(:,1) > 12 & Emax > 0.1;
idxDf8 = FScoords(:,1) < 16 & FScoords(:,1) > 14 & Emax > 0.1;
idxDf9 = FScoords(:,1) < 18 & FScoords(:,1) > 16 & Emax > 0.1;
idxDf10 = FScoords(:,1) < 20 & FScoords(:,1) > 18 & Emax > 0.1;

idxDf1x = FScoords(:,1) < 2 & FScoords(:,1) > 0;
idxDf2x = FScoords(:,1) < 4 & FScoords(:,1) > 2;
idxDf3x = FScoords(:,1) < 6 & FScoords(:,1) > 4;
idxDf4x = FScoords(:,1) < 8 & FScoords(:,1) > 6;
idxDf5x = FScoords(:,1) < 10 & FScoords(:,1) > 8;
idxDf6x = FScoords(:,1) < 12 & FScoords(:,1) > 10;
idxDf7x = FScoords(:,1) < 14 & FScoords(:,1) > 12;
idxDf8x = FScoords(:,1) < 16 & FScoords(:,1) > 14;
idxDf9x = FScoords(:,1) < 18 & FScoords(:,1) > 16;
idxDf10x = FScoords(:,1) < 20 & FScoords(:,1) > 18;

[PlD1,ElD1] = hist_dg(log10(Emax(idxDf1)),'range',[-3 3],'nbins',50);
[PlD2,ElD2] = hist_dg(log10(Emax(idxDf2)),'range',[-3 3],'nbins',50);
[PlD3,ElD3] = hist_dg(log10(Emax(idxDf3)),'range',[-3 3],'nbins',50);
[PlD4,ElD4] = hist_dg(log10(Emax(idxDf4)),'range',[-3 3],'nbins',50);
[PlD5,ElD5] = hist_dg(log10(Emax(idxDf5)),'range',[-3 3],'nbins',50);

logEdom = [-2.5 2.5];
logEvec = logEdom(1):0.01:logEdom(2);

meanlogE1 = mean(log10(Emax(idxDf1)));
stdlogE1 = std(log10(Emax(idxDf1)));
PSGTDf1 = 1/(sqrt(2*pi)*stdlogE1)*exp(-(logEvec - meanlogE1).^2/(2*stdlogE1^2));
PSGTDf1 = max(PlD1)/max(PSGTDf1)*PSGTDf1;

meanlogE2 = mean(log10(Emaxpar(idxDf2)));
stdlogE2 = std(log10(Emaxpar(idxDf2)));
PSGTDf2 = 1/(sqrt(2*pi)*stdlogE2)*exp(-(logEvec - meanlogE2).^2/(2*stdlogE2^2));
PSGTDf2 = max(PlD2)/max(PSGTDf2)*PSGTDf2;

meanlogE3 = mean(log10(Emax(idxDf3)));
stdlogE3 = std(log10(Emax(idxDf3)));
PSGTDf3 = 1/(sqrt(2*pi)*stdlogE3)*exp(-(logEvec - meanlogE3).^2/(2*stdlogE3^2));
PSGTDf3 = max(PlD3)/max(PSGTDf3)*PSGTDf3;

meanlogE4 = mean(log10(Emax(idxDf4)));
stdlogE4 = std(log10(Emax(idxDf4)));
PSGTDf4 = 1/(sqrt(2*pi)*stdlogE4)*exp(-(logEvec - meanlogE4).^2/(2*stdlogE4^2));
PSGTDf4 = max(PlD4)/max(PSGTDf4)*PSGTDf4;

meanlogE5 = mean(log10(Emax(idxDf5)));
stdlogE5 = std(log10(Emax(idxDf5)));
PSGTDf5 = 1/(sqrt(2*pi)*stdlogE5)*exp(-(logEvec - meanlogE5).^2/(2*stdlogE5^2));
PSGTDf5 = max(PlD5)/max(PSGTDf5)*PSGTDf5;

%parallel
[PlD1par,ElD1par] = hist_dg(log10(Emaxpar(idxDf1)),'range',[-3 3],'nbins',50);
[PlD2par,ElD2par] = hist_dg(log10(Emaxpar(idxDf2)),'range',[-3 3],'nbins',50);
[PlD3par,ElD3par] = hist_dg(log10(Emaxpar(idxDf3)),'range',[-3 3],'nbins',50);
[PlD4par,ElD4par] = hist_dg(log10(Emaxpar(idxDf4)),'range',[-3 3],'nbins',50);
[PlD5par,ElD5par] = hist_dg(log10(Emaxpar(idxDf5)),'range',[-3 3],'nbins',50);

logEdom = [-2.5 2.5];
logEvec = logEdom(1):0.01:logEdom(2);

meanlogE1par = mean(log10(Emaxpar(idxDf1)));
stdlogE1par = std(log10(Emaxpar(idxDf1)));
PSGTDf1par = 1/(sqrt(2*pi)*stdlogE1par)*exp(-(logEvec - meanlogE1par).^2/(2*stdlogE1par^2));
PSGTDf1par = max(PlD1par)/max(PSGTDf1par)*PSGTDf1par;

meanlogE2par = mean(log10(Emaxpar(idxDf2)));
stdlogE2par = std(log10(Emaxpar(idxDf2)));
PSGTDf2par = 1/(sqrt(2*pi)*stdlogE2par)*exp(-(logEvec - meanlogE2par).^2/(2*stdlogE2par^2));
PSGTDf2par = max(PlD2par)/max(PSGTDf2par)*PSGTDf2par;

meanlogE3par = mean(log10(Emaxpar(idxDf3)));
stdlogE3par = std(log10(Emaxpar(idxDf3)));
PSGTDf3par = 1/(sqrt(2*pi)*stdlogE3par)*exp(-(logEvec - meanlogE3par).^2/(2*stdlogE3par^2));
PSGTDf3par = max(PlD3par)/max(PSGTDf3par)*PSGTDf3par;

meanlogE4par = mean(log10(Emaxpar(idxDf4)));
stdlogE4par = std(log10(Emaxpar(idxDf4)));
PSGTDf4par = 1/(sqrt(2*pi)*stdlogE4par)*exp(-(logEvec - meanlogE4par).^2/(2*stdlogE4par^2));
PSGTDf4par = max(PlD4par)/max(PSGTDf4par)*PSGTDf4par;

meanlogE5par = mean(log10(Emaxpar(idxDf5)));
stdlogE5par = std(log10(Emaxpar(idxDf5)));
PSGTDf5par = 1/(sqrt(2*pi)*stdlogE5par)*exp(-(logEvec - meanlogE5par).^2/(2*stdlogE5par^2));
PSGTDf5par = max(PlD5par)/max(PSGTDf5par)*PSGTDf5par;

% Perpendicular
[PlD1perp,ElD1perp] = hist_dg(log10(Emaxperp(idxDf1)),'range',[-3 3],'nbins',50);
[PlD2perp,ElD2perp] = hist_dg(log10(Emaxperp(idxDf2)),'range',[-3 3],'nbins',50);
[PlD3perp,ElD3perp] = hist_dg(log10(Emaxperp(idxDf3)),'range',[-3 3],'nbins',50);
[PlD4perp,ElD4perp] = hist_dg(log10(Emaxperp(idxDf4)),'range',[-3 3],'nbins',50);
[PlD5perp,ElD5perp] = hist_dg(log10(Emaxperp(idxDf5)),'range',[-3 3],'nbins',50);

logEdom = [-2.5 2.5];
logEvec = logEdom(1):0.01:logEdom(2);

meanlogE1perp = mean(log10(Emaxperp(idxDf1)));
stdlogE1perp = std(log10(Emaxperp(idxDf1)));
PSGTDf1perp = 1/(sqrt(2*pi)*stdlogE1perp)*exp(-(logEvec - meanlogE1perp).^2/(2*stdlogE1perp^2));
PSGTDf1perp = max(PlD1perp)/max(PSGTDf1perp)*PSGTDf1perp;

meanlogE2perp = mean(log10(Emaxperp(idxDf2)));
stdlogE2perp = std(log10(Emaxperp(idxDf2)));
PSGTDf2perp = 1/(sqrt(2*pi)*stdlogE2perp)*exp(-(logEvec - meanlogE2perp).^2/(2*stdlogE2perp^2));
PSGTDf2perp = max(PlD2perp)/max(PSGTDf2perp)*PSGTDf2perp;

meanlogE3perp = mean(log10(Emaxperp(idxDf3)));
stdlogE3perp = std(log10(Emaxperp(idxDf3)));
PSGTDf3perp = 1/(sqrt(2*pi)*stdlogE3perp)*exp(-(logEvec - meanlogE3perp).^2/(2*stdlogE3perp^2));
PSGTDf3perp = max(PlD3perp)/max(PSGTDf3perp)*PSGTDf3perp;

meanlogE4perp = mean(log10(Emaxperp(idxDf4)));
stdlogE4perp = std(log10(Emaxperp(idxDf4)));
PSGTDf4perp = 1/(sqrt(2*pi)*stdlogE4perp)*exp(-(logEvec - meanlogE4perp).^2/(2*stdlogE4perp^2));
PSGTDf4perp = max(PlD4perp)/max(PSGTDf4perp)*PSGTDf4perp;

meanlogE5perp = mean(log10(Emaxperp(idxDf5)));
stdlogE5perp = std(log10(Emaxperp(idxDf5)));
PSGTDf5perp = 1/(sqrt(2*pi)*stdlogE5perp)*exp(-(logEvec - meanlogE5perp).^2/(2*stdlogE5perp^2));
PSGTDf5perp = max(PlD5perp)/max(PSGTDf5perp)*PSGTDf5perp;

color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.9290, 0.6940, 0.1250];
color4 = [0.4940, 0.1840, 0.5560];
color5 = [0.4660, 0.6740, 0.1880];
color6 = [0.3010, 0.7450, 0.9330];


%% Figure
fn=figure;
set(fn,'Position',[10 10 1000 600])
h(1)=axes('position',[0.05 0.59 0.275 0.4]); % [x y dx dy]
h(2)=axes('position',[0.385 0.59 0.275 0.4]); % [x y dx dy]
h(3)=axes('position',[0.72 0.59 0.275 0.4]); % [x y dx dy]
h(4)=axes('position',[0.05 0.090 0.275 0.4]); % [x y dx dy]
h(5)=axes('position',[0.385 0.090 0.275 0.4]); % [x y dx dy]
h(6)=axes('position',[0.72 0.090 0.275 0.4]); % [x y dx dy]
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 
set(fn,'defaultAxesFontSize',13)

loglog(h(1),10.^El,Pl,'linewidth',2,'color','k')
hold(h(1),'on')
loglog(h(1),10.^Elpar,Plpar,'linewidth',2,'color','r')
loglog(h(1),10.^Elperp,Plperp,'linewidth',2,'color','b')
loglog(h(1),10.^xx,10.^(ycross+grad*xx),'m--','linewidth',2)
loglog(h(1),10.^xx,10.^(ycross2+grad2*xx),'g--','linewidth',2)
hold(h(1),'off')
xlabel(h(1),'E_{max} (mV m^{-1})','interpreter','tex','fontsize',14)
ylabel(h(1),'P(log E_{max})','interpreter','tex','fontsize',14)
irf_legend(h(1),'\alpha_{||} = -0.7',[0.88, 0.98],'fontsize',14,'color','m');
irf_legend(h(1),'\alpha_{\perp} = -1.0',[0.88, 0.88],'fontsize',14,'color','g');
irf_legend(h(1),'(a)',[0.98, 0.98],'fontsize',14,'color','k');
axis(h(1),[1e-3 4e2 1 4e4])
legend(h(1),{'E','E_{||}','E_{\perp}'},'location','southwest')

loglog(h(2),10.^Elr,Plr,'linewidth',2,'color','k')
hold(h(2),'on')
loglog(h(2),10.^Elparr,Plparr,'linewidth',2,'color','r')
loglog(h(2),10.^Elperpr,Plperpr,'linewidth',2,'color','b')
%plot(xx,3+grad*xx,'r--','linewidth',2)
hold(h(2),'off')
xlabel(h(2),'E_{rms} (mV m^{-1})','interpreter','tex','fontsize',14)
ylabel(h(2),'P(log E_{rms})','interpreter','tex','fontsize',14)
irf_legend(h(2),'(b)',[0.98, 0.98],'fontsize',14,'color','k');
axis(h(2),[1e-3 4e2 1 4e4])

loglog(h(3),10.^Erat,Prat,'linewidth',2,'color','k')
hold(h(3),'on')
loglog(h(3),10.^Eparrat,Pparrat,'linewidth',2,'color','r')
loglog(h(3),10.^Eperprat,Pperprat,'linewidth',2,'color','b')
hold(h(3),'off')
set(h(3),'yscale','log')
xlabel(h(3),'E_{max}/E_{rms}','interpreter','tex','fontsize',14)
ylabel(h(3),'P(log E_{max}/E_{rms})','interpreter','tex','fontsize',14)
irf_legend(h(3),'(c)',[0.98, 0.98],'fontsize',14,'color','k');
irf_legend(h(3),'E_{max} > 1 mV m^{-1}',[0.02, 0.98],'fontsize',14,'color','k');
axis(h(3),[1 4e2 1 5e3])


loglog(h(4),10.^ElD1,PlD1,'linewidth',2)
hold(h(4),'on')
loglog(h(4),10.^ElD2,PlD2,'linewidth',2)
loglog(h(4),10.^ElD3,PlD3,'linewidth',2)
loglog(h(4),10.^ElD4,PlD4,'linewidth',2)
loglog(h(4),10.^ElD5,PlD5,'linewidth',2)
%loglog(h(4),10.^ElD6,PlD6,'linewidth',2)
loglog(h(4),10.^logEvec,PSGTDf1,'--','linewidth',2,'color',color1)
loglog(h(4),10.^logEvec,PSGTDf2,'--','linewidth',2,'color',color2)
loglog(h(4),10.^logEvec,PSGTDf3,'--','linewidth',2,'color',color3)
loglog(h(4),10.^logEvec,PSGTDf4,'--','linewidth',2,'color',color4)
loglog(h(4),10.^logEvec,PSGTDf5,'--','linewidth',2,'color',color5)
%loglog(h(4),10.^logEvec,PSGTDf6,'--','linewidth',2,'color',color6)
hold(h(4),'off')
xlabel(h(4),'E_{max} (mV m^{-1})','interpreter','tex','fontsize',14)
ylabel(h(4),'P(log E_{max})','interpreter','tex','fontsize',14)
axis(h(4),[1e-2 3e2 1 2e3])
legend(h(4),{'0 < D_f/R_E < 2','2 < D_f/R_E < 4','4 < D_f/R_E < 6','6 < D_f/R_E < 8','8 < D_f/R_E < 10'},'location','southwest')
irf_legend(h(4),'(d)',[0.99 0.99],'color','k','fontsize',14)

loglog(h(5),10.^ElD1par,PlD1par,'linewidth',2)
hold(h(5),'on')
loglog(h(5),10.^ElD2par,PlD2par,'linewidth',2)
loglog(h(5),10.^ElD3par,PlD3par,'linewidth',2)
loglog(h(5),10.^ElD4par,PlD4par,'linewidth',2)
loglog(h(5),10.^ElD5par,PlD5par,'linewidth',2)
%loglog(h(5),10.^ElD6,PlD6,'linewidth',2)
loglog(h(5),10.^logEvec,PSGTDf1par,'--','linewidth',2,'color',color1)
loglog(h(5),10.^logEvec,PSGTDf2par,'--','linewidth',2,'color',color2)
loglog(h(5),10.^logEvec,PSGTDf3par,'--','linewidth',2,'color',color3)
loglog(h(5),10.^logEvec,PSGTDf4par,'--','linewidth',2,'color',color4)
loglog(h(5),10.^logEvec,PSGTDf5par,'--','linewidth',2,'color',color5)
%loglog(h(5),10.^logEvec,PSGTDf6,'--','linewidth',2,'color',color6)
hold(h(5),'off')
xlabel(h(5),'E_{||,max} (mV m^{-1})','interpreter','tex','fontsize',14)
ylabel(h(5),'P(log E_{||,max})','interpreter','tex','fontsize',14)
axis(h(5),[1e-2 3e2 1 2e3])
irf_legend(h(5),'(e)',[0.99 0.99],'color','k','fontsize',14)

loglog(h(6),10.^ElD1perp,PlD1perp,'linewidth',2)
hold(h(6),'on')
loglog(h(6),10.^ElD2perp,PlD2perp,'linewidth',2)
loglog(h(6),10.^ElD3perp,PlD3perp,'linewidth',2)
loglog(h(6),10.^ElD4perp,PlD4perp,'linewidth',2)
loglog(h(6),10.^ElD5perp,PlD5perp,'linewidth',2)
%loglog(h(6),10.^ElD6,PlD6,'linewidth',2)
loglog(h(6),10.^logEvec,PSGTDf1perp,'--','linewidth',2,'color',color1)
loglog(h(6),10.^logEvec,PSGTDf2perp,'--','linewidth',2,'color',color2)
loglog(h(6),10.^logEvec,PSGTDf3perp,'--','linewidth',2,'color',color3)
loglog(h(6),10.^logEvec,PSGTDf4perp,'--','linewidth',2,'color',color4)
loglog(h(6),10.^logEvec,PSGTDf5perp,'--','linewidth',2,'color',color5)
%loglog(h(6),10.^logEvec,PSGTDf6,'--','linewidth',2,'color',color6)
hold(h(6),'off')
xlabel(h(6),'E_{\perp,max} (mV m^{-1})','interpreter','tex','fontsize',14)
ylabel(h(6),'P(log E_{\perp,max})','interpreter','tex','fontsize',14)
axis(h(6),[1e-2 3e2 1 2e3])
irf_legend(h(6),'(f)',[0.99 0.99],'color','k','fontsize',14)
%legend(h(6),{'-4 < D_f < -2','-2 < D_f < 0','0 < D_f < 2','2 < D_f < 4','4 < D_f < 6'},'location','northwest')
set(gcf,'color','w')