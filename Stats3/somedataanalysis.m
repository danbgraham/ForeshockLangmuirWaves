%%
load('sgtdata.mat');

mu1 = sgtdata.mu1;
mu2 = sgtdata.mu2;
mu3 = sgtdata.mu3;
mu4 = sgtdata.mu4;
mu5 = sgtdata.mu5;
mu6 = sgtdata.mu6;
mu7 = sgtdata.mu7;
mu8 = sgtdata.mu8;

F = sgtdata.F;

Chir2E = sgtdata.Chi2ESGTr.*sgtdata.NN.*sgtdata.difflogE;
Chir2ENL = sgtdata.Chi2ESGTNLr.*sgtdata.NN.*sgtdata.difflogE;

idx = Chir2ENL < Chir2E;
sum(idx);

Emax = sgtdata.Eenvmax;
EmaxES = Emax(idx);

logEc = sgtdata.logEfitESGTNL;

Ec = 10.^logEc(idx);

[PE,xE] = hist_dg(log10(Emax),'range',log10([5 300]),'nbins',20);
[PEES,xEES] = hist_dg(log10(EmaxES),'range',log10([5 300]),'nbins',20);


Chir2Epar = sgtdata.Chi2EparSGTr;
Chir2EparNL = sgtdata.Chi2EparSGTNLr;

idx = Chir2EparNL < Chir2Epar;

Emaxpar = sgtdata.Eenvparmax;
[PEpar,xEpar] = hist_dg(log10(Emaxpar),'range',log10([5 300]),'nbins',20);
[PEparES,xEparES] = hist_dg(log10(Emaxpar(idx)),'range',log10([5 300]),'nbins',20);


%% Calculate
load('sgtdata.mat');

mu1 = sgtdata.mu1;
mu2 = sgtdata.mu2;
mu3 = sgtdata.mu3;
mu4 = sgtdata.mu4;
mu5 = sgtdata.mu5;
mu6 = sgtdata.mu6;
mu7 = sgtdata.mu7;
mu8 = sgtdata.mu8;
NN = sgtdata.NN;

beta1 = mu3.^2./mu2.^3;
beta2 = mu4./mu2.^2;
beta3 = mu3.*mu5./mu2.^4;
beta4 = mu6./mu2.^3;
beta5 = mu7.*mu3./mu2.^5; 
beta6 = mu8./mu2.^4;


dbeta1 = beta1.*(4*beta4 - 24*beta2 + 36 +9*beta1.*beta2 -12*beta3 + 35*beta1);
dbeta2 = beta6-4*beta2.*beta4+4*beta2.^3-beta2.^2+16*beta2.*beta1 - 8*beta3 + 16*beta1;

dbeta1 = sqrt(dbeta1)./sqrt(NN);
dbeta2 = sqrt(dbeta2)./sqrt(NN);

idx = beta1 < 0.4 & beta2 < 3.25 & beta2 > 2.78;

idx2 = beta1-dbeta1 < 0.4 & beta2-dbeta2 < 3.25 & beta2+dbeta2 > 2.78;

%% 
Eenvparmedian = sgtdata.Eenvparmedian;
Eenvparmax = sgtdata.Eenvparmax;
Eenvparmin = sgtdata.Eenvparmin;

Eenvperpmedian = sgtdata.Eenvperpmedian;
Eenvperpmax = sgtdata.Eenvperpmax;
Eenvperpmin = sgtdata.Eenvperpmin;

idxperp = Eenvparmin > 0.01 & Eenvparmedian > 0.1 & Eenvperpmin > 0.01 & Eenvperpmedian > 0.1;

load('sgtdata.mat');

mu1 = sgtdata.mu1perp;
mu2 = sgtdata.mu2perp;
mu3 = sgtdata.mu3perp;
mu4 = sgtdata.mu4perp;
mu5 = sgtdata.mu5perp;
mu6 = sgtdata.mu6perp;
mu7 = sgtdata.mu7perp;
mu8 = sgtdata.mu8perp;
NN = sgtdata.NN;

c_eval('mu? = mu?(idxperp);',1:8);
NN = NN(idxperp);

beta1 = mu3.^2./mu2.^3;
beta2 = mu4./mu2.^2;
beta3 = mu3.*mu5./mu2.^4;
beta4 = mu6./mu2.^3;
beta5 = mu7.*mu3./mu2.^5; 
beta6 = mu8./mu2.^4;


dbeta1 = beta1.*(4*beta4 - 24*beta2 + 36 +9*beta1.*beta2 -12*beta3 + 35*beta1);
dbeta2 = beta6-4*beta2.*beta4+4*beta2.^3-beta2.^2+16*beta2.*beta1 - 8*beta3 + 16*beta1;

dbeta1 = sqrt(dbeta1)./sqrt(NN);
dbeta2 = sqrt(dbeta2)./sqrt(NN);

idx = beta1 < 0.4 & beta2 < 3.25 & beta2 > 2.78;

idx2 = beta1-dbeta1 < 0.4 & beta2-dbeta2 < 3.25 & beta2+dbeta2 > 2.78;