Bz = 2;
Dp = 2;
Vs = 400;
VA = 50;

as = 14.2*Dp^(-1/6);
bs = as*(27 + 84*VA^2/Vs^2)^(-2)*Dp^(-1/3);


Y = -50:0.01:50;
X = as - bs*Y.^2;


alpha = (0.58-0.01*Bz)*(1 + 0.01*Dp);
r0 = (11.4 + 0.14*Bz)*Dp^(-1/6.6);
theta = -160:0.01:160;

rr = r0*(2./(1 + cosd(theta))).^alpha;
Xmp = rr.*cosd(theta);
Ymp = rr.*sind(theta);


%%

fn=figure;

set(fn,'Position',[10 10 300 500])
    h(1)=axes('position',[0.15 0.1 0.80 0.87]); % [x y dx dy]
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',3);
    set(fn,'defaultAxesFontSize',16)
    
plot(h(1),X,Y,'k')
hold(h(1),'on')
plot(h(1),Xmp,Ymp,'k')
fill(h(1),[X flip(Xmp)],[Y flip(Ymp)],[0.8 0.8 0.8])
quiver(h(1),20,-20,00,40,'linewidth',1.0,'color',[0 0.5 0])
quiver(h(1),27,-20,00,40,'linewidth',1.0,'color',[0 0.5 0])
quiver(h(1),as,-20,00,40,'linewidth',1.0,'color',[0 0.5 0])
plot(h(1),[27 27],[-50 50],'linewidth',1,'color',[0 0.5 0])
plot(h(1),[20 20],[-50 50],'linewidth',1,'color',[0 0.5 0])
plot(h(1),[as as],[-50 50],'linewidth',1,'color',[0 0.5 0])
quiver(h(1),as,0,00,35,'linewidth',2.0,'color','r')
quiver(h(1),as,32,-as+7.5-0.5,0,'linewidth',2.0,'color','r')
quiver(h(1),5,-35,-2.2,8,'linewidth',0.5,'color','k')
annotation('textarrow',[0.3 0.4],[0.52 0.52],'String','V_{sw}','Linewidth',2,'color','b','fontsize',16)
set(h(1), 'XDir','reverse')
plot(h(1),7,32,'color','k','Marker','h','MarkerFaceColor','k')
plot(h(1),as,32,'color','r','Marker','o','MarkerFaceColor','r')
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(h(1),-xEarth,yEarth,'k','edgecolor','none')
patch(h(1),xEarth,yEarth,'w','edgecolor','k')
axis(h(1),[-10 30 -40 40])
hold(h(1),'off')
irf_legend(h(1),'B',[0.12 0.98],'color',[0 0.5 0],'fontsize',16)
irf_legend(h(1),'MSH',[0.9 0.85],'color',[0 0 0],'fontsize',16)
irf_legend(h(1),'R',[0.37 0.75],'color','r','fontsize',16)
irf_legend(h(1),'D_f',[0.53 0.96],'color','r','fontsize',16)
irf_legend(h(1),'Bowshock',[0.8 0.03],'color','k','fontsize',16)
xlabel(h(1),'x (R_E)','fontsize',16)
ylabel(h(1),'y (R_E)','fontsize',16)
%irf_legend(h(1),'V_{sw}',[0.05 0.5],'color','b','fontsize',14

set(gcf,'color','w')