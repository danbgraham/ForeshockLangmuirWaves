function sgtstruct = sgtstatdata(Efac,fpemedian)

sgtstruct = struct;
sgtstruct.starttime = Efac.time.start.epochUnix;
sgtstruct.stoptime = Efac.time.stop.epochUnix;

NN = length(Efac);
sgtstruct.NN = NN;

%% Frequency calculations
Efach = Efac;
Efach.data(:,1) = Efach.data(:,1).*hamming(NN);
Efach.data(:,2) = Efach.data(:,2).*hamming(NN);
Efach.data(:,3) = Efach.data(:,3).*hamming(NN);

[freqE,powerE] = fft_mms(Efach);
freqE = freqE';
powerE = powerE(:,1)+powerE(:,2)+powerE(:,3);

% Get peak frequency
idx = freqE > fpemedian/2;
df = median(diff(freqE));

fav = sum(df*freqE(idx).*powerE(idx))/sum(df*powerE(idx));
deltafav2 = sum(df*(freqE(idx)-fav).^2.*powerE(idx))/sum(df*powerE(idx));
deltafav = sqrt(deltafav2);

sgtstruct.fav = fav;
sgtstruct.deltafav = deltafav;


%% Get Envelope function
[Eenvpar,Eenvperp] = calcenvelopepp(Efac,fpemedian/1.5);
Eenv = calcenvelope(Efac,fpemedian/1.5);

%% Basic envelope properties
Eenvmin = min(Eenv.data);
Eenvmedian = median(Eenv.data);
Eenvmax = max(Eenv.data);

Eenvparmin = min(Eenvpar.data);
Eenvparmedian = median(Eenvpar.data);
Eenvparmax = max(Eenvpar.data);

Eenvperpmin = min(Eenvperp.data);
Eenvperpmedian = median(Eenvperp.data);
Eenvperpmax = max(Eenvperp.data);

sgtstruct.Eenvmin = Eenvmin; sgtstruct.Eenvmedian = Eenvmedian; sgtstruct.Eenvmax = Eenvmax;
sgtstruct.Eenvparmin = Eenvparmin; sgtstruct.Eenvparmedian = Eenvparmedian; sgtstruct.Eenvparmax = Eenvparmax;
sgtstruct.Eenvperpmin = Eenvperpmin; sgtstruct.Eenvperpmedian = Eenvperpmedian; sgtstruct.Eenvperpmax = Eenvperpmax;

%% Get moments of the distribution

% Total
mu1 = sum(log10(Eenv.data))/NN;
mu2 = sum((log10(Eenv.data)-mu1).^2)/NN;
mu3 = sum((log10(Eenv.data)-mu1).^3)/NN;
mu4 = sum((log10(Eenv.data)-mu1).^4)/NN;
mu5 = sum((log10(Eenv.data)-mu1).^5)/NN;
mu6 = sum((log10(Eenv.data)-mu1).^6)/NN;
mu7 = sum((log10(Eenv.data)-mu1).^7)/NN;
mu8 = sum((log10(Eenv.data)-mu1).^8)/NN;

sgtstruct.mu1 = mu1; sgtstruct.mu2 = mu2; sgtstruct.mu3 = mu3; sgtstruct.mu4 = mu4; 
sgtstruct.mu5 = mu5; sgtstruct.mu6 = mu6; sgtstruct.mu7 = mu7; sgtstruct.mu8 = mu8; 

% Parallel
mu1par = sum(log10(Eenvpar.data))/NN;
mu2par = sum((log10(Eenvpar.data)-mu1par).^2)/NN;
mu3par = sum((log10(Eenvpar.data)-mu1par).^3)/NN;
mu4par = sum((log10(Eenvpar.data)-mu1par).^4)/NN;
mu5par = sum((log10(Eenvpar.data)-mu1par).^5)/NN;
mu6par = sum((log10(Eenvpar.data)-mu1par).^6)/NN;
mu7par = sum((log10(Eenvpar.data)-mu1par).^7)/NN;
mu8par = sum((log10(Eenvpar.data)-mu1par).^8)/NN;

sgtstruct.mu1par = mu1par; sgtstruct.mu2par = mu2par; sgtstruct.mu3par = mu3par; sgtstruct.mu4par = mu4par; 
sgtstruct.mu5par = mu5par; sgtstruct.mu6par = mu6par; sgtstruct.mu7par = mu7par; sgtstruct.mu8par = mu8par; 

% Perpendicular
mu1perp = sum(log10(Eenvperp.data))/NN;
mu2perp = sum((log10(Eenvperp.data)-mu1perp).^2)/NN;
mu3perp = sum((log10(Eenvperp.data)-mu1perp).^3)/NN;
mu4perp = sum((log10(Eenvperp.data)-mu1perp).^4)/NN;
mu5perp = sum((log10(Eenvperp.data)-mu1perp).^5)/NN;
mu6perp = sum((log10(Eenvperp.data)-mu1perp).^6)/NN;
mu7perp = sum((log10(Eenvperp.data)-mu1perp).^7)/NN;
mu8perp = sum((log10(Eenvperp.data)-mu1perp).^8)/NN; 

sgtstruct.mu1perp = mu1perp; sgtstruct.mu2perp = mu2perp; sgtstruct.mu3perp = mu3perp; sgtstruct.mu4perp = mu4perp; 
sgtstruct.mu5perp = mu5perp; sgtstruct.mu6perp = mu6perp; sgtstruct.mu7perp = mu7perp; sgtstruct.mu8perp = mu8perp; 

%% Make envelope distribution functions
[PlogE,logE] = hist_dg(log10(Eenv.data),'range',[log10(Eenvmin) log10(Eenvmax)],'nbins',30);
[PlogEpar,logEpar] = hist_dg(log10(Eenvpar.data),'range',[log10(Eenvparmin) log10(Eenvparmax)],'nbins',30);
[PlogEperp,logEperp] = hist_dg(log10(Eenvperp.data),'range',[log10(Eenvperpmin) log10(Eenvperpmax)],'nbins',30);

%NN = length(log10(Eenv.data(1:10:end));
diffE = median(diff(logE));
PlogE = PlogE/(NN*diffE);

diffEpar = median(diff(logEpar));
PlogEpar = PlogEpar/(NN*diffEpar);

diffEperp = median(diff(logEperp));
PlogEperp = PlogEperp/(NN*diffEperp);

sgtstruct.difflogE = diffE; sgtstruct.difflogEpar = diffEpar; sgtstruct.difflogEperp = diffEperp; 
sgtstruct.PlogE = PlogE; sgtstruct.PlogEpar = PlogEpar; sgtstruct.PlogEperp = PlogEperp; 
sgtstruct.logE = logE; sgtstruct.logEpar = logEpar; sgtstruct.logEperp = logEperp; 



%% Fit simple SGT to envelop statistics
% Define functions
SGTfun = @(mn,sigma,ampl,x) ( ampl.*exp(-(x-mn).^2./(2*sigma.^2)) ); 
%SGTNLfun = @(mn,sigma,ampl,logEc,x) ( ampl.*( exp(-(x-mn).^2./(2*sigma.^2)) - exp(-(2*logEc - x - mn).^2./(2*sigma.^2)) ) ); 

fSGT = @(x) sum( (PlogE - SGTfun(x(1),x(2),x(3),logE)).^2./SGTfun(x(1),x(2),x(3),logE),'omitnan');
fSGTpar = @(x) sum( (PlogEpar - SGTfun(x(1),x(2),x(3),logEpar)).^2./SGTfun(x(1),x(2),x(3),logEpar),'omitnan');
fSGTperp = @(x) sum( (PlogEperp - SGTfun(x(1),x(2),x(3),logEperp)).^2./SGTfun(x(1),x(2),x(3),logEperp),'omitnan');


fSGTNL = @(x) sum( (PlogE - SGTNLfun(x(1),x(2),x(3),x(4),logE)).^2./SGTNLfun(x(1),x(2),x(3),x(4),logE),'omitnan');
fSGTNLpar = @(x) sum( (PlogEpar - SGTNLfun(x(1),x(2),x(3),x(4),logEpar)).^2./SGTNLfun(x(1),x(2),x(3),x(4),logEpar),'omitnan');
fSGTNLperp = @(x) sum( (PlogEperp - SGTNLfun(x(1),x(2),x(3),x(4),logEperp)).^2./SGTNLfun(x(1),x(2),x(3),x(4),logEperp),'omitnan');

options = optimset('MaxFunEvals',1000000,'MaxIter',100000);

% Total electric field
[varSGT,fvalSGT,~,~] = fminsearch(@(x) fSGT(x),[mu1 sqrt(mu2) max(PlogE)],options);
[varSGT2,fvalSGT2,~,~] = fminsearch(@(x) fSGT(x),[mu1 sqrt(mu2) 1/(sqrt(mu2)*sqrt(2*pi))],options);

if fvalSGT2 < fvalSGT
  varSGT = varSGT2;
end


mufitESGT = varSGT(1);
sigmafitESGT = varSGT(2);
AmpfitESGT = varSGT(3);

Modseries = SGTfun(mufitESGT,sigmafitESGT,AmpfitESGT,logE);
Chiseries = (PlogE-Modseries).^2./Modseries;
Chi2ESGT = sum(Chiseries,'omitnan');
kESGT = 30-4;
Chi2ESGTr = Chi2ESGT/kESGT;

sgtstruct.mufitESGT = mufitESGT;
sgtstruct.sigmafitESGT = sigmafitESGT;
sgtstruct.AmpfitESGT = AmpfitESGT;
sgtstruct.Chi2ESGT = Chi2ESGT;
sgtstruct.Chi2ESGTr = Chi2ESGTr;

[varSGTNL,fvalSGTNL,~,~] = fminsearch(@(x) fSGTNL(x),[mufitESGT sigmafitESGT AmpfitESGT log10(Eenvmax)],options);

[varSGTNL2,fvalSGTNL2,~,~] = fminsearch(@(x) fSGTNL(x),[mu1 sqrt(mu2) 1/(sqrt(mu2)*sqrt(2*pi)) log10(Eenvmax)+0.2],options);

[varSGTNL3,fvalSGTNL3,~,~] = fminsearch(@(x) fSGTNL(x),[mufitESGT sigmafitESGT AmpfitESGT log10(Eenvmax)+2],options);

if fvalSGTNL2 < fvalSGTNL
  varSGTNL = varSGTNL2;
  fvalSGTNL = fvalSGTNL2;
end
if fvalSGTNL3 < fvalSGTNL
  varSGTNL = varSGTNL3;
end

mufitESGTNL = varSGTNL(1);
sigmafitESGTNL = varSGTNL(2);
AmpfitESGTNL = varSGTNL(3);
logEfitESGTNL = varSGTNL(4);

Modseries = SGTNLfun(mufitESGTNL,sigmafitESGTNL,AmpfitESGTNL,logEfitESGTNL,logE);
idx = Modseries < 0;
Chiseries = (PlogE-Modseries).^2./Modseries;
Chi2ESGTNL = sum(Chiseries(~idx),'omitnan');
kESGT = 30-5-sum(idx);
Chi2ESGTNLr = Chi2ESGTNL/kESGT;

sgtstruct.mufitESGTNL = mufitESGTNL;
sgtstruct.sigmafitESGTNL = sigmafitESGTNL;
sgtstruct.AmpfitESGTNL = AmpfitESGTNL;
sgtstruct.logEfitESGTNL = logEfitESGTNL;
sgtstruct.Chi2ESGTNL = Chi2ESGTNL;
sgtstruct.Chi2ESGTNLr = Chi2ESGTNLr;

% Parallel electric field
[varSGTpar,fvarSGTpar,~,~] = fminsearch(@(x) fSGTpar(x),[mu1par sqrt(mu2par) max(PlogEpar)],options);

[varSGTpar2,fvarSGTpar2,~,~] = fminsearch(@(x) fSGTpar(x),[mu1par sqrt(mu2par) 1/(sqrt(mu2par)*sqrt(2*pi))],options);

if fvarSGTpar2 < fvarSGTpar
  varSGTpar = varSGTpar2;
end

mufitEparSGT = varSGTpar(1);
sigmafitEparSGT = varSGTpar(2);
AmpfitEparSGT = varSGTpar(3);

Modseries = SGTfun(mufitEparSGT,sigmafitEparSGT,AmpfitEparSGT,logEpar);
Chiseries = (PlogEpar-Modseries).^2./Modseries;
Chi2EparSGT = sum(Chiseries,'omitnan');
kESGT = 30-4;
Chi2EparSGTr = Chi2EparSGT/kESGT;

sgtstruct.mufitEparSGT = mufitEparSGT;
sgtstruct.sigmafitEparSGT = sigmafitEparSGT;
sgtstruct.AmpfitEparSGT = AmpfitEparSGT;
sgtstruct.Chi2EparSGT = Chi2EparSGT;
sgtstruct.Chi2EparSGTr = Chi2EparSGTr;

[varSGTNLpar,fvalSGTNLpar,~,~] = fminsearch(@(x) fSGTNLpar(x),[mufitEparSGT sigmafitEparSGT AmpfitEparSGT log10(Eenvparmax)],options);

[varSGTNLpar2,fvalSGTNLpar2,~,~] = fminsearch(@(x) fSGTNLpar(x),[mu1par sqrt(mu2par) 1/(sqrt(mu2par)*sqrt(2*pi)) log10(Eenvparmax)+0.2],options);

[varSGTNLpar3,fvalSGTNLpar3,~,~] = fminsearch(@(x) fSGTNLpar(x),[mufitEparSGT sigmafitEparSGT AmpfitEparSGT log10(Eenvparmax)+2],options);

if fvalSGTNLpar2 < fvalSGTNLpar
  varSGTNLpar = varSGTNLpar2;
  fvalSGTNLpar = fvalSGTNLpar2;
end
if fvalSGTNLpar3 < fvalSGTNLpar
  varSGTNLpar = varSGTNLpar3;
end

mufitEparSGTNL = varSGTNLpar(1);
sigmafitEparSGTNL = varSGTNLpar(2);
AmpfitEparSGTNL = varSGTNLpar(3);
logEfitEparSGTNL = varSGTNLpar(4);

Modseries = SGTNLfun(mufitEparSGTNL,sigmafitEparSGTNL,AmpfitEparSGTNL,logEfitEparSGTNL,logEpar);
idx = Modseries < 0;
Chiseries = (PlogEpar-Modseries).^2./Modseries;
Chi2EparSGTNL = sum(Chiseries(~idx),'omitnan');
kESGT = 30-5-sum(idx);
Chi2EparSGTNLr = Chi2EparSGTNL/kESGT;

sgtstruct.mufitEparSGTNL = mufitEparSGTNL;
sgtstruct.sigmafitEparSGTNL = sigmafitEparSGTNL;
sgtstruct.AmpfitEparSGTNL = AmpfitEparSGTNL;
sgtstruct.logEfitEparSGTNL = logEfitEparSGTNL;
sgtstruct.Chi2EparSGTNL = Chi2EparSGTNL;
sgtstruct.Chi2EparSGTNLr = Chi2EparSGTNLr;

% Perpendicular electric field
[varSGTperp,fvarSGTperp,~,~] = fminsearch(@(x) fSGTperp(x),[mu1perp sqrt(mu2perp) max(PlogEperp)],options);

[varSGTperp2,fvarSGTperp2,~,~] = fminsearch(@(x) fSGTperp(x),[mu1perp sqrt(mu2perp) 1/(sqrt(mu2perp)*sqrt(2*pi))],options);

if fvarSGTperp2 < fvarSGTperp
  varSGTperp = varSGTperp2;
end

mufitEperpSGT = varSGTperp(1);
sigmafitEperpSGT = varSGTperp(2);
AmpfitEperpSGT = varSGTperp(3);

Modseries = SGTfun(mufitEperpSGT,sigmafitEperpSGT,AmpfitEperpSGT,logEperp);
Chiseries = (PlogEperp-Modseries).^2./Modseries;
Chi2EperpSGT = sum(Chiseries,'omitnan');
kESGT = 30-4;
Chi2EperpSGTr = Chi2EperpSGT/kESGT;

sgtstruct.mufitEperpSGT = mufitEperpSGT;
sgtstruct.sigmafitEperpSGT = sigmafitEperpSGT;
sgtstruct.AmpfitEperpSGT = AmpfitEperpSGT;
sgtstruct.Chi2EperpSGT = Chi2EperpSGT;
sgtstruct.Chi2EperpSGTr = Chi2EperpSGTr;

[varSGTNLperp,fvalSGTNLperp,~,~] = fminsearch(@(x) fSGTNLperp(x),[mufitEperpSGT sigmafitEperpSGT AmpfitEperpSGT log10(Eenvperpmax)],options);

[varSGTNLperp2,fvalSGTNLperp2,~,~] = fminsearch(@(x) fSGTNLperp(x),[mu1perp sqrt(mu2perp) 1/(sqrt(mu2perp)*sqrt(2*pi)) log10(Eenvperpmax)+0.2],options);

[varSGTNLperp3,fvalSGTNLperp3,~,~] = fminsearch(@(x) fSGTNLperp(x),[mufitEperpSGT sigmafitEperpSGT AmpfitEperpSGT log10(Eenvperpmax)+2],options);

if fvalSGTNLperp2 < fvalSGTNLperp
  varSGTNLperp = varSGTNLperp2;
  fvalSGTNLperp = fvalSGTNLperp2;
end
if fvalSGTNLperp3 < fvalSGTNLperp
  varSGTNLperp = varSGTNLperp3;
end

mufitEperpSGTNL = varSGTNLperp(1);
sigmafitEperpSGTNL = varSGTNLperp(2);
AmpfitEperpSGTNL = varSGTNLperp(3);
logEfitEperpSGTNL = varSGTNLperp(4);

Modseries = SGTNLfun(mufitEperpSGTNL,sigmafitEperpSGTNL,AmpfitEperpSGTNL,logEfitEperpSGTNL,logEperp);
idx = Modseries < 0;
Chiseries = (PlogEperp-Modseries).^2./Modseries;
Chi2EperpSGTNL = sum(Chiseries(~idx),'omitnan');
kESGT = 30-5-sum(idx);
Chi2EperpSGTNLr = Chi2EperpSGTNL/kESGT;

sgtstruct.mufitEperpSGTNL = mufitEperpSGTNL;
sgtstruct.sigmafitEperpSGTNL = sigmafitEperpSGTNL;
sgtstruct.AmpfitEperpSGTNL = AmpfitEperpSGTNL;
sgtstruct.logEfitEperpSGTNL = logEfitEperpSGTNL;
sgtstruct.Chi2EperpSGTNL = Chi2EperpSGTNL;
sgtstruct.Chi2EperpSGTNLr = Chi2EperpSGTNLr;

end


function val = SGTNLfun(mn,sigma,ampl,logEc,x)
  val = ampl.*( exp(-(x-mn).^2./(2*sigma.^2)) - exp(-(2*logEc - x - mn).^2./(2*sigma.^2)) );
  idx = val < 0;
  val(idx) = 0;
end
