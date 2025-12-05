function hfdata = loaddata(ic)
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here
scid = ['MMS' num2str(ic)];

d=dir(['hfdata/hfwaves_' scid '_*.mat']);

starttimes = [];

for ii=1:length(d)
  load(['hfdata/' d(ii).name]);
  if isempty(starttimes)
      starttimes = data.starttimes.epochUnix;
      endtimes = data.endtimes.epochUnix;
      F = data.F;
      lmax = data.lmax;
      lint = data.lint;
      lmin = data.lmin;
      thmax = data.thmax;
      thint = data.thint;
      thmin = data.thmin;
      fpemedian = data.fpemedian;
      fpemax = data.fpemax;
      fcemedian = data.fcemedian;
      fcemax = data.fcemax;
      lambdaDmedian = data.lambdaDmedian;
      lambdaDmax = data.lambdaDmax;
      rhoemedian = data.rhoemedian;
      rhoemax = data.rhoemax;
      Wmax = data.Wmax;
      Emax = data.Emax;
      fpeak = data.fpeak;
      nemedian = data.nemedian;
      Temedian = data.Temedian;
      nemax = data.nemax;
      Temax = data.Temax;
      nimedian = data.nimedian;
      Vimedian = data.Vimedian;
      Timedian = data.Timedian;
      Posgse = data.Posgse;
      Posgsm = data.Posgsm;
      Bgsemedian = data.Bgsemedian;
      Bgsmmedian = data.Bgsmmedian;
      Bgsemax = data.Bgsemax;
      Bgsmmax = data.Bgsmmax;
      Bswgse = data.Bswgse;
      Bswgsm = data.Bswgsm;
      Vsw = data.Vsw;
      nsw = data.nsw;
      Psw = data.Psw;
      EDIflag = data.EDIflag;
      ASPOCflag = data.ASPOCflag;
      regiontype = data.regiontype;
      timefromSW = data.timefromSW;
  else
      starttimes = [starttimes; data.starttimes.epochUnix];
      endtimes = [endtimes; data.endtimes.epochUnix];
      F = [F; data.F];
      lmax = [lmax; data.lmax];
      lint = [lint; data.lint];
      lmin = [lmin; data.lmin];
      thmax = [thmax; data.thmax];
      thint = [thint; data.thint];
      thmin = [thmin; data.thmin];
      fpemedian = [fpemedian; data.fpemedian];
      fpemax = [fpemax; data.fpemax];
      fcemedian = [fcemedian; data.fcemedian];
      fcemax = [fcemax; data.fcemax];
      lambdaDmedian = [lambdaDmedian; data.lambdaDmedian];
      lambdaDmax = [lambdaDmax; data.lambdaDmax];
      rhoemedian = [rhoemedian; data.rhoemedian];
      rhoemax = [rhoemax; data.rhoemax];
      Wmax = [Wmax; data.Wmax];
      Emax = [Emax; data.Emax];
      fpeak = [fpeak; data.fpeak];
      nemedian = [nemedian; data.nemedian];
      Temedian = [Temedian; data.Temedian];
      nemax = [nemax; data.nemax];
      Temax = [Temax; data.Temax];
      nimedian = [nimedian; data.nimedian];
      Vimedian = [Vimedian; data.Vimedian];
      Timedian = [Timedian; data.Timedian];
      Posgse = [Posgse; data.Posgse];
      Posgsm = [Posgsm; data.Posgsm];
      Bgsemedian = [Bgsemedian; data.Bgsemedian];
      Bgsmmedian = [Bgsmmedian; data.Bgsmmedian];
      Bgsemax = [Bgsemax; data.Bgsemax];
      Bgsmmax = [Bgsmmax; data.Bgsmmax];
      Bswgse = [Bswgse; data.Bswgse];
      Bswgsm = [Bswgsm; data.Bswgsm];
      Vsw = [Vsw; data.Vsw];
      nsw = [nsw; data.nsw];
      Psw = [Psw; data.Psw];
      EDIflag = [EDIflag; data.EDIflag];
      ASPOCflag = [ASPOCflag; data.ASPOCflag];
      regiontype = [regiontype; data.regiontype];
      timefromSW = [timefromSW; data.timefromSW];
  end
  clear data
  ii
end

hfdata = struct('starttimes',starttimes,'endtimes',endtimes,'F',F,...
    'lmax',lmax,'lint',lint,'lmin',lmin,'thmax',thmax,'thint',thint,'thmin',thmin,...
    'fpemedian',fpemedian,'fpemax',fpemax,'fcemedian',fcemedian,'fcemax',fcemax,...
    'lambdaDmedian',lambdaDmedian,'lambdaDmax',lambdaDmax,'rhoemedian',rhoemedian,'rhoemax',rhoemax,...
    'Wmax',Wmax,'Emax',Emax,'fpeak',fpeak,'nemedian',nemedian,'Temedian',Temedian,...
    'nemax',nemax,'Temax',Temax,'nimedian',nimedian,'Vimedian',Vimedian,'Timedian',Timedian,...
    'Posgse',Posgse,'Posgsm',Posgsm,'Bgsemedian',Bgsemedian,'Bgsmmedian',Bgsmmedian,...
    'Bgsemax',Bgsemax,'Bgsmmax',Bgsmmax,'Bswgse',Bswgse,'Bswgsm',Bswgsm,...
    'Vsw',Vsw,'nsw',nsw,'Psw',Psw,'EDIflag',EDIflag,'ASPOCflag',ASPOCflag,...
    'regiontype',regiontype,'timefromSW',timefromSW);

end

