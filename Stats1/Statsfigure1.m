%% Load data
MMS1data = loaddata(1);
MMS2data = loaddata(2);
MMS3data = loaddata(3);
MMS4data = loaddata(4);

%% 
starttimes = [MMS1data.starttimes; MMS2data.starttimes; MMS3data.starttimes; MMS4data.starttimes];
endtimes = [MMS1data.endtimes; MMS2data.endtimes; MMS3data.endtimes; MMS4data.endtimes];
F = [MMS1data.F; MMS2data.F; MMS3data.F; MMS4data.F];
lmax = [MMS1data.lmax; MMS2data.lmax; MMS3data.lmax; MMS4data.lmax];
lint = [MMS1data.lint; MMS2data.lint; MMS3data.lint; MMS4data.lint];
lmin = [MMS1data.lmin; MMS2data.lmin; MMS3data.lmin; MMS4data.lmin];
thmax = [MMS1data.thmax; MMS2data.thmax; MMS3data.thmax; MMS4data.thmax];
thint = [MMS1data.thint; MMS2data.thint; MMS3data.thint; MMS4data.thint];
thmin = [MMS1data.thmin; MMS2data.thmin; MMS3data.thmin; MMS4data.thmin];
fpemedian = [MMS1data.fpemedian; MMS2data.fpemedian; MMS3data.fpemedian; MMS4data.fpemedian];
fpemax = [MMS1data.fpemax; MMS2data.fpemax; MMS3data.fpemax; MMS4data.fpemax];
fcemedian = [MMS1data.fcemedian; MMS2data.fcemedian; MMS3data.fcemedian; MMS4data.fcemedian];
fcemax = [MMS1data.fcemax; MMS2data.fcemax; MMS3data.fcemax; MMS4data.fcemax];
lambdaDmedian = [MMS1data.lambdaDmedian; MMS2data.lambdaDmedian; MMS3data.lambdaDmedian; MMS4data.lambdaDmedian];
lambdaDmax = [MMS1data.lambdaDmax; MMS2data.lambdaDmax; MMS3data.lambdaDmax; MMS4data.lambdaDmax];
rhoemedian = [MMS1data.rhoemedian; MMS2data.rhoemedian; MMS3data.rhoemedian; MMS4data.rhoemedian];
rhoemax = [MMS1data.rhoemax; MMS2data.rhoemax; MMS3data.rhoemax; MMS4data.rhoemax];
Wmax = [MMS1data.Wmax; MMS2data.Wmax; MMS3data.Wmax; MMS4data.Wmax];
Emax = [MMS1data.Emax; MMS2data.Emax; MMS3data.Emax; MMS4data.Emax];
fpeak = [MMS1data.fpeak; MMS2data.fpeak; MMS3data.fpeak; MMS4data.fpeak];
nemedian = [MMS1data.nemedian; MMS2data.nemedian; MMS3data.nemedian; MMS4data.nemedian];
Temedian = [MMS1data.Temedian; MMS2data.Temedian; MMS3data.Temedian; MMS4data.Temedian];
nemax = [MMS1data.nemax; MMS2data.nemax; MMS3data.nemax; MMS4data.nemax];
Temax = [MMS1data.Temax; MMS2data.Temax; MMS3data.Temax; MMS4data.Temax];
nimedian = [MMS1data.nimedian; MMS2data.nimedian; MMS3data.nimedian; MMS4data.nimedian];
Vimedian = [MMS1data.Vimedian; MMS2data.Vimedian; MMS3data.Vimedian; MMS4data.Vimedian];
Timedian = [MMS1data.Timedian; MMS2data.Timedian; MMS3data.Timedian; MMS4data.Timedian];
Posgse = [MMS1data.Posgse; MMS2data.Posgse; MMS3data.Posgse; MMS4data.Posgse];
Posgsm = [MMS1data.Posgsm; MMS2data.Posgsm; MMS3data.Posgsm; MMS4data.Posgse];
Bgsemedian = [MMS1data.Bgsemedian; MMS2data.Bgsemedian; MMS3data.Bgsemedian; MMS4data.Bgsemedian];
Bgsmmedian = [MMS1data.Bgsmmedian; MMS2data.Bgsmmedian; MMS3data.Bgsmmedian; MMS4data.Bgsmmedian];
Bgsemax = [MMS1data.Bgsemax; MMS2data.Bgsemax; MMS3data.Bgsemax; MMS4data.Bgsemax];
Bgsmmax = [MMS1data.Bgsmmax; MMS2data.Bgsmmax; MMS3data.Bgsmmax; MMS4data.Bgsmmax];
Bswgse = [MMS1data.Bswgse; MMS2data.Bswgse; MMS3data.Bswgse; MMS4data.Bswgse];
Bswgsm = [MMS1data.Bswgsm; MMS2data.Bswgsm; MMS3data.Bswgsm; MMS4data.Bswgsm];
Vsw = [MMS1data.Vsw; MMS2data.Vsw; MMS3data.Vsw; MMS4data.Vsw];
nsw = [MMS1data.nsw; MMS2data.nsw; MMS3data.nsw; MMS4data.nsw];
Psw = [MMS1data.Psw; MMS2data.Psw; MMS3data.Psw; MMS4data.Psw];
ASPOCflag = [MMS1data.ASPOCflag; MMS2data.ASPOCflag; MMS3data.ASPOCflag; MMS4data.ASPOCflag];
EDIflag = [MMS1data.EDIflag; MMS2data.EDIflag; MMS3data.EDIflag; MMS4data.EDIflag];
regiontype = [MMS1data.regiontype; MMS2data.regiontype; MMS3data.regiontype; MMS4data.regiontype];
timefromSW = [MMS1data.timefromSW; MMS2data.timefromSW; MMS3data.timefromSW; MMS4data.timefromSW];

Bmedian = sqrt(Bgsemax(:,1).^2+Bgsemax(:,2).^2+Bgsemax(:,3).^2);

%% Solar wind plots
Units = irf_units;
RE = Units.R_Earth/1e3;
PosgsmRE = Posgsm/RE;
Ryz = sqrt(PosgsmRE(:,2).^2+PosgsmRE(:,3).^2);

idxUH = F > 0.5;
idxL = F < 0.5;

idxtailL = PosgsmRE(:,1) < 0 & Ryz < 10 & F < 0.5;
idxtailUH = PosgsmRE(:,1) < 0 & Ryz < 10 & F > 0.5;

idxB = (abs(Bgsemax(:,1)-Bgsemedian(:,1))+abs(Bgsemax(:,2)-Bgsemedian(:,2))+abs(Bgsemax(:,3)-Bgsemedian(:,3)))./Bmedian < 0.1;
idxSW = regiontype == 1; % & idxB;
idxtail = PosgsmRE(:,1) < -5 & Ryz < 10;
idxMS = ~idxtail & ~idxSW;

[PFSW,FSW] = hist_dg(F(idxSW),'range',[0 1],'nbins',100);
[PFMS,FMS] = hist_dg(F(idxMS),'range',[0 1],'nbins',100);
[PFtail,Ftail] = hist_dg(F(idxtail),'range',[0 1],'nbins',100);

[FS,vbcS,thetaS] = getpropertiesSTEREO;
[PFS,xPS] = hist_dg(FS,'range',[0 1],'nbins',50);

%% Make figure
fn=figure;
set(fn,'Position',[10 10 350 500])
h(1)=axes('position',[0.15 0.56 0.83 0.40]); % [x y dx dy]
h(2)=axes('position',[0.12 0.08 0.85 0.40]); % [x y dx dy]
%h(3)=axes('position',[0.07 0.08 0.40 0.4]); % [x y dx dy]
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 
set(fn,'defaultAxesFontSize',13)

load('swsurf.mat');

zlimvals = [0.99 1.05];
surf(h(1),swsurf.kperp,swsurf.kar,real(swsurf.freq),swsurf.FE)
set(h(1),'ylim',[0 1e0]);
set(h(1),'xlim',[0 1e0]);
set(h(1),'yscale','log');
set(h(1),'xscale','log');
set(h(1),'xtick',[1e-4 1e-3 1e-2 1e-1 1]);
set(h(1),'ytick',[1e-4 1e-3 1e-2 1e-1 1]);
hold(h(1),'on')
plot3(h(1),swsurf.kperp(1)*ones(length(swsurf.kar)),swsurf.kar,real(swsurf.freq(:,1)),'color','k','linewidth',1.5)
plot3(h(1),swsurf.kperp,swsurf.kar(1)*ones(length(swsurf.kperp)),real(swsurf.freq(1,:)),'color','k','linewidth',1.5)
hold(h(1),'off')
zlim(h(1),zlimvals)
clim(h(1),[0 1])
ylabel(h(1),'k_{||}\lambda_D','fontsize',12)
xlabel(h(1),'k_{\perp}\lambda_D','fontsize',12)
zlabel(h(1),'f/f_{pe}','fontsize',12)
title(h(1),'F_E','fontsize',12)
set(h(1),'fontsize',12);
c=colorbar('peer',h(1),'ver');
rr = interp1([1 88 128 168 256],[0.2  0.0 0.25 1.0 1.0],1:256);
gg = interp1([1 88 128 168 256],[0.2  1.0 1.0 1.0 0.2],1:256);
bb = interp1([1 88 128 168 256],[1.0   1.0 0.25 0.0 0.2],1:256);
jet2 = [rr' gg' bb'];  
colormap(h(1),jet2);
irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',14)
xtickangle(h(1),0)
ytickangle(h(1),0)

set(gcf,'color','w')

plot(h(2),FMS,PFMS/max(PFMS),'linewidth',2,'color','r')
hold(h(2),'on')
%plot(h(2),xPS,PFS/max(PFS),'linewidth',2,'color',[0.0 0.7 0.0])
plot(h(2),FSW,PFSW/max(PFSW),'linewidth',2,'color','k')
plot(h(2),Ftail,PFtail/max(PFtail),'linewidth',2,'color',[0 0.7 0])
plot(h(2),xPS,PFS/max(PFS),'linewidth',2,'color',[0.3 0.3 1])
plot(h(2),FSW,PFSW/max(PFSW),'linewidth',2,'color','k')
hold(h(2),'off')
xlabel(h(2),'F_E','interpreter','tex','fontsize',12)
ylabel(h(2),'Counts (normalized)','interpreter','tex','fontsize',12)
irf_legend(h(2),'Electron foreshock',[0.93 0.98],'fontsize',12,'color','k')
irf_legend(h(2),'Magnetopause',[0.93 0.90],'fontsize',12,'color','r')
irf_legend(h(2),'Magnetotail',[0.93 0.82],'fontsize',12,'color',[0 0.7 0])
irf_legend(h(2),'STEREO: Type III',[0.93 0.74],'fontsize',12,'color',[0.3 0.3 1])
%irf_legend(h(2),'STEREO: Type III',[0.3 0.74],'fontsize',12,'color',[0 0.7 0])
set(h(2),'fontsize',12)
irf_legend(h(2),'(b)',[0.2 0.99],'color','k','fontsize',14)

set(gcf,'color','w')