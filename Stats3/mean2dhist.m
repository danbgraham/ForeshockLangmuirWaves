function histstruct = mean2dhist(xparam,yparam,epsilon,xrange,yrange,numpnts)

xvec = linspace(xrange(1),xrange(2),numpnts);
yvec = linspace(yrange(1),yrange(2),numpnts);

dx = median(diff(xvec))/2;
dy = median(diff(yvec))/2;

xedges = [xvec(1)-dx xvec+dx];
yedges = [yvec(1)-dy yvec+dy];

countshist = zeros(numpnts,numpnts);
meanhist = zeros(numpnts,numpnts);
medianhist = zeros(numpnts,numpnts);
stdhist = zeros(numpnts,numpnts);

tic
for ii=1:numpnts
  for jj=1:numpnts
    idxy = yparam > yedges(ii) & yparam < yedges(ii+1);
    idxx = xparam > xedges(jj) & xparam < xedges(jj+1);
    idx = idxy & idxx;
    countshist(ii,jj) = sum(idx);
    meanhist(ii,jj) = mean(epsilon(idx),'omitnan');
    medianhist(ii,jj) = median(epsilon(idx),'omitnan');
    stdhist(ii,jj) = std(epsilon(idx),'omitnan');
  end
end
toc

histstruct = struct('xvec',xvec,'yvec',yvec,'xedges',xedges,'yedges',yedges,'countshist',countshist,...
  'meanhist',meanhist,'medianhist',medianhist,'stdhist',stdhist);

end

