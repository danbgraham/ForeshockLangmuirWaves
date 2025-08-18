function Eenv = calcenvelope(Exyz,fmin)


dfE = 1/median(diff(Exyz.time.epochUnix));
Exyzf = Exyz.filt(fmin,0,dfE,5);

E1 = irf.ts_scalar(Exyzf.time,Exyzf.data(:,1));
E2 = irf.ts_scalar(Exyzf.time,Exyzf.data(:,2));
E3 = irf.ts_scalar(Exyzf.time,Exyzf.data(:,3));
E1h = imag(hilbert(E1.data));
E2h = imag(hilbert(E2.data));
E3h = imag(hilbert(E3.data));
Eenv1 = sqrt(E1.data.^2+E1h.^2);
Eenv2 = sqrt(E2.data.^2+E2h.^2);
Eenv3 = sqrt(E3.data.^2+E3h.^2);
Eenv1 = smooth(Eenv1,15);
Eenv2 = smooth(Eenv2,15);
Eenv3 = smooth(Eenv3,15);
Eenvtot = sqrt(Eenv1.^2+Eenv2.^2+Eenv3.^2);
Eenv = irf.ts_scalar(Exyz.time,Eenvtot);

end

