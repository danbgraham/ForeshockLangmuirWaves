load('sgtdata.mat');

NN = sgtdata.NN;
dlogE = sgtdata.difflogE;
dlogEpar = sgtdata.difflogEpar;
dlogEperp = sgtdata.difflogEperp;

idxperp = sgtdata.Eenvparmin > 0.01 & sgtdata.Eenvparmedian > 0.1 & sgtdata.Eenvperpmin > 0.01 & sgtdata.Eenvperpmedian > 0.1;

Chi2SGTr = sgtdata.Chi2ESGTr;
Chi2SGTNLr = sgtdata.Chi2ESGTNLr;

Chi2SGTr = Chi2SGTr.*NN.*dlogE; % Use Chi2 as defined for Chi2 test
Chi2SGTNLr = Chi2SGTNLr.*NN.*dlogE;

Chi2SGTparr = sgtdata.Chi2EparSGTr(idxperp);
Chi2SGTNLparr = sgtdata.Chi2EparSGTNLr(idxperp);

Chi2SGTparr = Chi2SGTparr.*NN(idxperp).*dlogEpar(idxperp); % Use Chi2 as defined for Chi2 test
Chi2SGTNLparr = Chi2SGTNLparr.*NN(idxperp).*dlogEpar(idxperp);

Chi2SGTperpr = sgtdata.Chi2EperpSGTr(idxperp);
Chi2SGTNLperpr = sgtdata.Chi2EperpSGTNLr(idxperp);

Chi2SGTperpr = Chi2SGTperpr.*NN(idxperp).*dlogEperp(idxperp); % Use Chi2 as defined for Chi2 test
Chi2SGTNLperpr = Chi2SGTNLperpr.*NN(idxperp).*dlogEperp(idxperp);

%% Calculate histograms
histChi2SGTtot = mean2dhist(log10(Chi2SGTr),log10(Chi2SGTNLr),log10(Chi2SGTNLr),[0.5 4.5],[0.5 4.5],80);
countsSGTtot = histChi2SGTtot.countshist;
countsSGTtot(countsSGTtot == 0) = NaN;
xvecSGTtot = 10.^histChi2SGTtot.xvec;
yvecSGTtot = 10.^histChi2SGTtot.yvec;

histChi2SGTpar = mean2dhist(log10(Chi2SGTparr),log10(Chi2SGTNLparr),log10(Chi2SGTNLparr),[0.5 4.5],[0.5 4.5],80);
countsSGTpar = histChi2SGTpar.countshist;
countsSGTpar(countsSGTpar == 0) = NaN;
xvecSGTpar = 10.^histChi2SGTpar.xvec;
yvecSGTpar = 10.^histChi2SGTpar.yvec;

histChi2SGTperp = mean2dhist(log10(Chi2SGTperpr),log10(Chi2SGTNLperpr),log10(Chi2SGTNLperpr),[0.5 4.5],[0.5 4.5],80);
countsSGTperp = histChi2SGTperp.countshist;
countsSGTperp(countsSGTperp == 0) = NaN;
xvecSGTperp = 10.^histChi2SGTperp.xvec;
yvecSGTperp = 10.^histChi2SGTperp.yvec;

xline = 10.^(0:0.1:5);
yline = xline;

%% 
fn=figure;
xwidth = 0.25;
ywidth = 0.73;

set(fn,'Position',[10 10 1000 300])
    h(1)=axes('position',[0.06 0.18 xwidth ywidth]); % [x y dx dy]
    h(2)=axes('position',[0.395 0.18 xwidth ywidth]);
    h(3)=axes('position',[0.73 0.18 xwidth ywidth]);
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2);
    set(fn,'defaultAxesFontSize',16);

pcolor(h(1),xvecSGTtot,yvecSGTtot,countsSGTtot);
shading(h(1),'flat');
c=colorbar('peer',h(1),'ver');
c.Label.String = 'Counts';
set(h(1),'xscale','log')
set(h(1),'yscale','log')
hold(h(1),'on')
plot(h(1),xline,yline,'r--')
hold(h(1),'off')
xlabel(h(1),'\chi^2_{r,SGT}','Interpreter','tex')
ylabel(h(1),'\chi^2_{r,SGTNL}','Interpreter','tex')
irf_legend(h(1),'(a)',[0.90, 0.99],'color','k','fontsize',16);
xticks(h(1),[1e0 1e1 1e2 1e3 1e4])
title(h(1),'E_{env}')

pcolor(h(2),xvecSGTpar,yvecSGTpar,countsSGTpar);
shading(h(2),'flat');
c=colorbar('peer',h(2),'ver');
c.Label.String = 'Counts';
set(h(2),'xscale','log')
set(h(2),'yscale','log')
hold(h(2),'on')
plot(h(2),xline,yline,'r--')
hold(h(2),'off')
xlabel(h(2),'\chi^2_{r,SGT}','Interpreter','tex')
ylabel(h(2),'\chi^2_{r,SGTNL}','Interpreter','tex')
irf_legend(h(2),'(b)',[0.90, 0.99],'color','k','fontsize',16);
xticks(h(2),[1e0 1e1 1e2 1e3 1e4])
title(h(2),'E_{||,env}')

pcolor(h(3),xvecSGTperp,yvecSGTperp,countsSGTperp);
shading(h(3),'flat');
c=colorbar('peer',h(3),'ver');
c.Label.String = 'Counts';
set(h(3),'xscale','log')
set(h(3),'yscale','log')
hold(h(3),'on')
plot(h(3),xline,yline,'r--')
hold(h(3),'off')
xlabel(h(3),'\chi^2_{r,SGT}','Interpreter','tex')
ylabel(h(3),'\chi^2_{r,SGTNL}','Interpreter','tex')
irf_legend(h(3),'(c)',[0.90, 0.99],'color','k','fontsize',16);
xticks(h(3),[1e0 1e1 1e2 1e3 1e4])
title(h(3),'E_{\perp,env}')

set(gcf,'color','w')