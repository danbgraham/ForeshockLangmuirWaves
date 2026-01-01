%% Load data

load('pearsondataEtot.mat'); 

NN = pearsondataEtot.NN;
dlogE = pearsondataEtot.difflogE;

Chi2SGTmomsr = pearsondataEtot.Chi2SGTmomsr;
Chi2SGTfitr = pearsondataEtot.Chi2SGTfitr;
Chi2pearsmomsr = pearsondataEtot.Chi2pearsmomsr;
Chi2pearsfitr = pearsondataEtot.Chi2pearsfitr;

Chi2SGTmomsr = Chi2SGTmomsr.*NN.*dlogE; % Use Chi2 as defined for Chi2 test
Chi2SGTfitr = Chi2SGTfitr.*NN.*dlogE;
Chi2pearsmomsr = Chi2pearsmomsr.*NN.*dlogE;
Chi2pearsfitr = Chi2pearsfitr.*NN.*dlogE;

Chi2SGTmoms = pearsondataEtot.Chi2SGTmoms;
Chi2SGTfit = pearsondataEtot.Chi2SGTfit;

%% Calculate 2D histograms
histChi2SGT = mean2dhist(log10(Chi2SGTmomsr),log10(Chi2SGTfitr),log10(Chi2SGTfitr),[0.5 4.5],[0.5 4.5],80);
countsSGT = histChi2SGT.countshist;
countsSGT(countsSGT == 0) = NaN;
xvecSGT = 10.^histChi2SGT.xvec;
yvecSGT = 10.^histChi2SGT.yvec;

histChi2pears = mean2dhist(log10(Chi2pearsmomsr),log10(Chi2pearsfitr),log10(Chi2pearsfitr),[0.5 4.5],[0.5 4.5],80);
countspears = histChi2pears.countshist;
countspears(countspears == 0) = NaN;
xvecpears = 10.^histChi2pears.xvec;
yvecpears = 10.^histChi2pears.yvec;

histChi2moms = mean2dhist(log10(Chi2SGTmomsr),log10(Chi2pearsmomsr),log10(Chi2pearsmomsr),[0.5 4.5],[0.5 4.5],80);
countsmoms = histChi2moms.countshist;
countsmoms(countsmoms == 0) = NaN;
xvecmoms = 10.^histChi2moms.xvec;
yvecmoms = 10.^histChi2moms.yvec;

histChi2fit = mean2dhist(log10(Chi2SGTfitr),log10(Chi2pearsfitr),log10(Chi2pearsfitr),[0.5 4.5],[0.5 4.5],80);
countsfit = histChi2fit.countshist;
countsfit(countsfit == 0) = NaN;
xvecfit = 10.^histChi2fit.xvec;
yvecfit = 10.^histChi2fit.yvec;


xline = 10.^(0:0.1:5);
yline = xline;

%% Plot Figure
fn=figure;
xwidth = 0.375;
ywidth = 0.38;

set(fn,'Position',[10 10 750 600])
    h(1)=axes('position',[0.085 0.615 xwidth ywidth]); % [x y dx dy]
    h(2)=axes('position',[0.59 0.615 xwidth ywidth]);
    h(3)=axes('position',[0.085 0.11 xwidth ywidth]);
    h(4)=axes('position',[0.59 0.11 xwidth ywidth]);
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2);
    set(fn,'defaultAxesFontSize',18);

pcolor(h(1),xvecSGT,yvecSGT,countsSGT);
shading(h(1),'flat');
c=colorbar('peer',h(1),'ver');
c.Label.String = 'Counts';
set(h(1),'xscale','log')
set(h(1),'yscale','log')
hold(h(1),'on')
plot(h(1),xline,yline,'r--')
hold(h(1),'off')
xlabel(h(1),'\chi^2_{r,SGT} moms','Interpreter','tex')
ylabel(h(1),'\chi^2_{r,SGT} fit','Interpreter','tex')
xticks(h(1),[1e0 1e1 1e2 1e3 1e4])
irf_legend(h(1),'(a)',[0.90, 0.99],'color','k','fontsize',18);

pcolor(h(2),xvecpears,yvecpears,countspears);
shading(h(2),'flat');
c=colorbar('peer',h(2),'ver');
c.Label.String = 'Counts';
set(h(2),'xscale','log')
set(h(2),'yscale','log')
hold(h(2),'on')
plot(h(2),xline,yline,'r--')
hold(h(2),'off')
xlabel(h(2),'\chi^2_{r,P} moms','Interpreter','tex')
ylabel(h(2),'\chi^2_{r,P} fit','Interpreter','tex')
xticks(h(2),[1e0 1e1 1e2 1e3 1e4])
irf_legend(h(2),'(b)',[0.90, 0.99],'color','k','fontsize',18);

pcolor(h(3),xvecmoms,yvecmoms,countsmoms);
shading(h(3),'flat');
c=colorbar('peer',h(3),'ver');
c.Label.String = 'Counts';
set(h(3),'xscale','log')
set(h(3),'yscale','log')
hold(h(3),'on')
plot(h(3),xline,yline,'r--')
hold(h(3),'off')
xlabel(h(3),'\chi^2_{r,P} moms','Interpreter','tex')
ylabel(h(3),'\chi^2_{r,SGT} mons','Interpreter','tex')
xticks(h(3),[1e0 1e1 1e2 1e3 1e4])
irf_legend(h(3),'(c)',[0.90, 0.99],'color','k','fontsize',18);

pcolor(h(4),xvecfit,yvecfit,countsfit);
shading(h(4),'flat');
c=colorbar('peer',h(4),'ver');
c.Label.String = 'Counts';
set(h(4),'xscale','log')
set(h(4),'yscale','log')
hold(h(4),'on')
plot(h(4),xline,yline,'r--')
hold(h(4),'off')
xlabel(h(4),'\chi^2_{r,P} fit','Interpreter','tex')
ylabel(h(4),'\chi^2_{r,SGT} fit','Interpreter','tex')
xticks(h(4),[1e0 1e1 1e2 1e3 1e4])
irf_legend(h(4),'(d)',[0.90, 0.99],'color','k','fontsize',18);


set(gcf,'color','w')