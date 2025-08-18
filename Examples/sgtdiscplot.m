dx = 0.001;
x = [-4:dx:5];

mu = 0.5;
sigma = 0.5;

y = exp(-(x-mu).^2/(2*sigma^2));

int0 = sum(y)*dx;

int1 = sum(x.*y)*dx/int0;
int2 = sum((x-mu).^2.*y)*dx/int0;
int3 = sum((x-mu).^3.*y)*dx/int0;
int4 = sum((x-mu).^4.*y)*dx/int0;

sk1 = int3/int2^(3/2);
sk2 = int3^2/int2^3;
kurt = int4/int2^2;

sigmamod = sqrt(int2);

%%
logEc = mu+sigma*1;
x2 = [-8:dx:logEc];
A = erf((logEc - mu)/(sigma*sqrt(2)));
y2 = 1/(sqrt(2*pi)*sigma*A)*(exp(-(x2 - mu).^2/(2*sigma^2)) - exp(-(2*logEc - x2 - mu).^2/(2*sigma^2)));

int02 = sum(y2)*dx;

int12 = sum(x2.*y2)*dx/int02;
int22 = sum((x2-int12).^2.*y2)*dx/int02;
int32 = sum((x2-int12).^3.*y2)*dx/int02;
int42 = sum((x2-int12).^4.*y2)*dx/int02;

sk12 = int32/int22.^(3/2);
sk22 = int32^2/int22^3
kurt2 = int42/int22^2


mup = mu + (logEc-mu)*(1 - A)/(A);

%%
mu = 0.0;
sigma = 0.5;
dx = 0.001;

xxx = [0.02:0.0103:6];

logEcvec = xxx*sigma + mu;

kurt = zeros(size(xxx));
skew1 = zeros(size(xxx));
skew2 = zeros(size(xxx));
mup = zeros(size(xxx));
sigmap = zeros(size(xxx));

for ii = 1:length(xxx)
  xvec = [-10:dx:logEcvec(ii)];
  A = erf((logEcvec(ii) - mu)/(sigma*sqrt(2)));
  yvec = 1/(sqrt(2*pi)*sigma*A)*(exp(-(xvec - mu).^2/(2*sigma^2)) - exp(-(2*logEcvec(ii) - xvec - mu).^2/(2*sigma^2)));
  
  int12 = sum(xvec.*yvec)*dx;
  int22 = sum((xvec-int12).^2.*yvec)*dx;
  int32 = sum((xvec-int12).^3.*yvec)*dx;
  int42 = sum((xvec-int12).^4.*yvec)*dx;
  
  skew2(ii) = int32^2/int22^3;
  skew1(ii) = int32/int22^(3/2);
  kurt(ii) = int42/int22^2;
  mup(ii) = int12;
  sigmap(ii) = sqrt(int22);
  
end

xxx2 = (logEcvec-mup)./sigmap;


%%
h=irf_plot(2,'newfigure');
xSize=700; ySize=280;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.41;
ywidth = 0.80;
set(h(1),'position',[0.08 0.17 xwidth ywidth]);
set(h(2),'position',[0.58 0.17 xwidth ywidth]);

plot(h(1),xxx,skew2,'k','linewidth',2)
hold(h(1),'on')
plot(h(1),xxx2,skew2,'r','linewidth',2)
hold(h(1),'off')
xlabel(h(1),'(log E_c - \mu)/\sigma','interpreter','tex')
ylabel(h(1),'\beta_1')
irf_legend(h(1),'(a)',[0.98 0.97],'color','k')

plot(h(2),xxx,kurt,'k','linewidth',2)
hold(h(2),'on')
plot(h(2),xxx2,kurt,'r','linewidth',2)
hold(h(2),'off')
xlabel(h(2),'(log E_c - \mu)/\sigma','interpreter','tex')
ylabel(h(2),'\beta_2')
irf_legend(h(2),'(b)',[0.98 0.97],'color','k')