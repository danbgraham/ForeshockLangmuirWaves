%% Load data
load('sgtdata.mat');

mu1 = sgtdata.mu1;
mu2 = sgtdata.mu2;
mu3 = sgtdata.mu3;
mu4 = sgtdata.mu4;

F = sgtdata.F;



%% Data analysis
beta1 = mu3.^2./mu2.^3;
beta2 = mu4./mu2.^2;

kappa = beta1.*(beta2+3).^2./(4*(2*beta2 - 3*beta1-6).*(4*beta2 - 3*beta1));

beta2t1 = 1*beta1;
beta2t3 = (3*beta1+6)/2;
at5 = beta1 - 32; 
bt5 = 78*beta1+96;
ct5 = -36*beta1.^2 - 63*beta1;
beta2t5 = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);

idxistypeI = beta2 < beta2t3;
idxistypeVI = beta2 > beta2t3 & beta2 < beta2t5;
idxistypeIV = beta2 > beta2t5;

histstructbeta12 = mean2dhist(beta1,beta2,beta2,[0 2],[1 5],50);

[Pbeta1,xbeta1] = hist_dg(beta1,'range',[0 2],'nbins',100);
[Pbeta2,xbeta2] = hist_dg(beta2,'range',[1 5],'nbins',100);

[Pbeta1log,xbeta1log] = hist_dg(log10(beta1),'range',[-5 1],'nbins',100);

%idxlowF = F < 0.2;
%[Pbeta1f1,xbeta1f1] = hist_dg(beta1(idxlowF),'range',[0 2],'nbins',100);
%[Pbeta2f1,xbeta2f1] = hist_dg(beta2(idxlowF),'range',[1 5],'nbins',100);
%[Pbeta1f2,xbeta1f2] = hist_dg(beta1(~idxlowF),'range',[0 2],'nbins',100);
%[Pbeta2f2,xbeta2f2] = hist_dg(beta2(~idxlowF),'range',[1 5],'nbins',100);
%[Pbeta1logf1,xbeta1logf1] = hist_dg(log10(beta1(idxlowF)),'range',[-5 1],'nbins',100);
%[Pbeta1logf2,xbeta1logf2] = hist_dg(log10(beta1(~idxlowF)),'range',[-5 1],'nbins',100);

%% Figure: Statistics of the total electric field.
fn=figure;
xwidth = 0.38;
ywidth = 0.38;

set(fn,'Position',[10 10 600 500])
    h(1)=axes('position',[0.06 0.58 xwidth+0.03 ywidth]); % [x y dx dy]
    h(2)=axes('position',[0.58 0.58 xwidth ywidth]);
    h(3)=axes('position',[0.085 0.095 xwidth ywidth]);
    h(4)=axes('position',[0.58 0.095 xwidth ywidth]);
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2);
    set(fn,'defaultAxesFontSize',16)

beta1b1 = [0 2];
beta2b1 = [1 3];

beta1t3 = [0:0.1:2];
beta2t3 = (3*beta1t3+6)/2;

beta1t5 = [0:0.01:2];
at5 = beta1t5 - 32; 
bt5 = 78*beta1t5+96;
ct5 = -36*beta1t5.^2 - 63*beta1t5;

beta2t5 = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);
    
countsb12 = histstructbeta12.countshist;
countsb12(countsb12 == 0) = NaN;
pcolor(h(1),histstructbeta12.xvec,histstructbeta12.yvec,countsb12);
shading(h(1),'flat');
c=colorbar('peer',h(1),'ver');
c.Label.String = 'Counts';
xlabel(h(1),'\beta_1','interpreter','tex')
ylabel(h(1),'\beta_2','interpreter','tex')
irf_legend(h(1),'(a)',[0.98 0.94],'fontsize',16)
irf_legend(h(1),'I',[0.7 0.5],'fontsize',16)
irf_legend(h(1),'VI',[0.56 1.00],'fontsize',16,'color','r')
irf_legend(h(1),'IV',[0.3 1.00],'fontsize',16,'color',[0 0.7 0])
hold(h(1),'on')
plot(h(1),beta1b1,beta2b1,'k','Linewidth',2)
plot(h(1),beta1t3,beta2t3,'r','Linewidth',2)
plot(h(1),beta1t5,beta2t5,'Linewidth',2,'color',[0 0.7 0])
fill(h(1), [0 2 2 0],[1 3 0 0],[0.7 0.7 0.7]);
hold(h(1),'off')

plot(h(2),xbeta2,Pbeta2,'k','linewidth',2)
xlabel(h(2),'\beta_2','interpreter','tex','fontsize',16)
ylabel(h(2),'Counts','fontsize',16)
irf_legend(h(2),'(b)',[0.98 0.98],'fontsize',16)
axis(h(2),[1 5 0 200])

plot(h(3),xbeta1,Pbeta1,'k','linewidth',2)
xlabel(h(3),'\beta_1','interpreter','tex','fontsize',16)
ylabel(h(3),'Counts','fontsize',16)
irf_legend(h(3),'(c)',[0.98 0.98],'fontsize',16)

plot(h(4),xbeta1log,Pbeta1log,'k','linewidth',2)
xlabel(h(4),'log_{10}(\beta_1)','interpreter','tex','fontsize',16)
ylabel(h(4),'Counts','fontsize',16)
irf_legend(h(4),'(d)',[0.98 0.98],'fontsize',16)
axis(h(4),[-5 1 0 200])

set(gcf,'color','w')

%% Parallel electric field analysis 
Eenvparmedian = sgtdata.Eenvparmedian;
Eenvparmax = sgtdata.Eenvparmax;
Eenvparmin = sgtdata.Eenvparmin;

Eenvperpmedian = sgtdata.Eenvperpmedian;
Eenvperpmax = sgtdata.Eenvperpmax;
Eenvperpmin = sgtdata.Eenvperpmin;

idxperp = Eenvparmin > 0.01 & Eenvparmedian > 0.1 & Eenvperpmin > 0.01 & Eenvperpmedian > 0.1;

mu1par = sgtdata.mu1par(idxperp);
mu2par = sgtdata.mu2par(idxperp);
mu3par = sgtdata.mu3par(idxperp);
mu4par = sgtdata.mu4par(idxperp);

beta1par = mu3par.^2./mu2par.^3;
beta2par = mu4par./mu2par.^2;

kappapar = beta1par.*(beta2par+3).^2./(4*(2*beta2par - 3*beta1par-6).*(4*beta2par - 3*beta1par));

beta2t1par = 1*beta1par;
beta2t3par = (3*beta1par+6)/2;
at5 = beta1par - 32; 
bt5 = 78*beta1par+96;
ct5 = -36*beta1par.^2 - 63*beta1par;
beta2t5par = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);

idxistypeIpar = beta2par < beta2t3par;
idxistypeVIpar = beta2par > beta2t3par & beta2par < beta2t5par;
idxistypeIVpar = beta2par > beta2t5par;

beta2t1 = 1*beta1;
beta2t3 = (3*beta1+6)/2;
at5 = beta1 - 32; 
bt5 = 78*beta1+96;
ct5 = -36*beta1.^2 - 63*beta1;
beta2t5 = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);

idxistypeI = beta2 < beta2t3;
idxistypeVI = beta2 > beta2t3 & beta2 < beta2t5;
idxistypeIV = beta2 > beta2t5;

histstructbeta12par = mean2dhist(beta1par,beta2par,beta2par,[0 2],[1 5],50);

[Pbeta1par,xbeta1par] = hist_dg(beta1par,'range',[0 2],'nbins',50);
[Pbeta2par,xbeta2par] = hist_dg(beta2par,'range',[1 5],'nbins',50);

[Pbeta2,xbeta2] = hist_dg(beta2,'range',[1 5],'nbins',50);

[Pbeta1logpar,xbeta1logpar] = hist_dg(log10(beta1par),'range',[-5 1],'nbins',50);

[Pbeta1log,xbeta1log] = hist_dg(log10(beta1),'range',[-5 1],'nbins',50);

%% Perpendicular electric field analysis
mu1perp = sgtdata.mu1perp(idxperp);
mu2perp = sgtdata.mu2perp(idxperp);
mu3perp = sgtdata.mu3perp(idxperp);
mu4perp = sgtdata.mu4perp(idxperp);

beta1perp = mu3perp.^2./mu2perp.^3;
beta2perp = mu4perp./mu2perp.^2;

kappaperp = beta1perp.*(beta2perp+3).^2./(4*(2*beta2perp - 3*beta1perp-6).*(4*beta2perp - 3*beta1perp));

beta2t1perp = 1*beta1perp;
beta2t3perp = (3*beta1perp+6)/2;
at5 = beta1perp - 32; 
bt5 = 78*beta1perp+96;
ct5 = -36*beta1perp.^2 - 63*beta1perp;
beta2t5perp = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);

idxistypeIperp = beta2perp < beta2t3perp;
idxistypeVIperp = beta2perp > beta2t3perp & beta2perp < beta2t5perp;
idxistypeIVperp = beta2perp > beta2t5perp;

histstructbeta12perp = mean2dhist(beta1perp,beta2perp,beta2perp,[0 2],[1 5],50);

[Pbeta1perp,xbeta1perp] = hist_dg(beta1perp,'range',[0 2],'nbins',50);
[Pbeta2perp,xbeta2perp] = hist_dg(beta2perp,'range',[1 5],'nbins',50);

[Pbeta1logperp,xbeta1logperp] = hist_dg(log10(beta1perp),'range',[-5 1],'nbins',50);

%% Figure: Statistics of the parallel and perpendicular components of the electric field.
fn=figure;
xwidth = 0.39;
ywidth = 0.38;

set(fn,'Position',[10 10 600 520])
    h(1)=axes('position',[0.06 0.57 xwidth+0.03 ywidth]); % [x y dx dy]
    h(2)=axes('position',[0.555 0.57 xwidth+0.03 ywidth]);
    h(3)=axes('position',[0.085 0.095 xwidth ywidth]);
    h(4)=axes('position',[0.58 0.095 xwidth ywidth]);
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2);
    set(fn,'defaultAxesFontSize',16)

beta1b1 = [0 2];
beta2b1 = [1 3];

beta1t3 = [0:0.1:2];
beta2t3 = (3*beta1t3+6)/2;

beta1t5 = [0:0.01:2];
at5 = beta1t5 - 32; 
bt5 = 78*beta1t5+96;
ct5 = -36*beta1t5.^2 - 63*beta1t5;

beta2t5 = (-bt5 - sqrt(bt5.^2 - 4*at5.*ct5))./(2*at5);
    
countsb12 = histstructbeta12par.countshist;
countsb12(countsb12 == 0) = NaN;
pcolor(h(1),histstructbeta12par.xvec,histstructbeta12par.yvec,countsb12);
shading(h(1),'flat');
c=colorbar('peer',h(1),'ver');
c.Label.String = 'Counts';
xlabel(h(1),'\beta_1','interpreter','tex')
ylabel(h(1),'\beta_2','interpreter','tex')
irf_legend(h(1),'(a)',[0.98 0.94],'fontsize',16)
irf_legend(h(1),'I',[0.7 0.5],'fontsize',16)
irf_legend(h(1),'VI',[0.56 1.00],'fontsize',16,'color','r')
irf_legend(h(1),'IV',[0.3 1.00],'fontsize',16,'color',[0 0.7 0])
hold(h(1),'on')
plot(h(1),beta1b1,beta2b1,'k','Linewidth',2)
plot(h(1),beta1t3,beta2t3,'r','Linewidth',2)
plot(h(1),beta1t5,beta2t5,'Linewidth',2,'color',[0 0.7 0])
fill(h(1), [0 2 2 0],[1 3 0 0],[0.7 0.7 0.7]);
hold(h(1),'off')
title(h(1),'E_{||}','interpreter','tex')

countsb12 = histstructbeta12perp.countshist;
countsb12(countsb12 == 0) = NaN;
pcolor(h(2),histstructbeta12perp.xvec,histstructbeta12perp.yvec,countsb12);
shading(h(2),'flat');
c=colorbar('peer',h(2),'ver');
c.Label.String = 'Counts';
xlabel(h(2),'\beta_1','interpreter','tex')
ylabel(h(2),'\beta_2','interpreter','tex')
irf_legend(h(2),'(b)',[0.98 0.94],'fontsize',16)
irf_legend(h(2),'I',[0.7 0.5],'fontsize',16)
irf_legend(h(2),'VI',[0.56 1.00],'fontsize',16,'color','r')
irf_legend(h(2),'IV',[0.3 1.00],'fontsize',16,'color',[0 0.7 0])
hold(h(2),'on')
plot(h(2),beta1b1,beta2b1,'k','Linewidth',2)
plot(h(2),beta1t3,beta2t3,'r','Linewidth',2)
plot(h(2),beta1t5,beta2t5,'Linewidth',2,'color',[0 0.7 0])
fill(h(2), [0 2 2 0],[1 3 0 0],[0.7 0.7 0.7]);
hold(h(2),'off')
title(h(2),'E_{\perp}','interpreter','tex')

plot(h(3),xbeta1logpar,Pbeta1logpar,'r','linewidth',2)
hold(h(3),'on')
plot(h(3),xbeta1logperp,Pbeta1logperp,'b','linewidth',2)
hold(h(3),'off')
xlabel(h(3),'log_{10}(\beta_1)','interpreter','tex','fontsize',16)
ylabel(h(3),'Counts','fontsize',16)
irf_legend(h(3),'(c)',[0.98 0.98],'fontsize',16)
axis(h(3),[-5 1 0 200])
legend(h(3),{'E_{||}','E_{\perp}'},'location','northwest','interpreter','tex','fontsize',16)

plot(h(4),xbeta2par,Pbeta2par,'r','linewidth',2)
hold(h(4),'on')
plot(h(4),xbeta2perp,Pbeta2perp,'b','linewidth',2)
hold(h(4),'off')
xlabel(h(4),'\beta_2','interpreter','tex','fontsize',16)
ylabel(h(4),'Counts','fontsize',16)
irf_legend(h(4),'(d)',[0.98 0.98],'fontsize',16)
axis(h(4),[1 5 0 200])

set(gcf,'color','w')

%%
