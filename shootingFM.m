clear;
close all;
addpath('MyFM/');
tic
load('Random2');

dt0=phi1;
dt1=phi2;

shift=0.00;
dt1=dt1+shift;

d=dt0;

M=0.1*ones(Nx,Ny);

map=[linspace(1,0.5,100);linspace(1,0,100);linspace(1,0,100)]';

nick='ResultShootingRandom2';

hFig = figure(1);
set(hFig, 'Position', [0 0 620 512],'Color',[1 1 1]);
axes1 = axes('Parent',hFig);
nn=0;

Mold=zeros(Nx,Ny);
crit=1;

while(crit>1e-3)
    
    t=0;
    nn=nn+1;
    % M(find(abs(dt0)>=5*dx))=0;
    k=0
    
    
    [M,d,~]=initBand(M,d,zeros(Nx,Ny),dx,dy);
    [d,M]=FastMarchingEikonalMobilityNewton(d,M,dx,dy);
    
    crit=sum(abs(M(:)-Mold(:)))/sum(abs(M(:)));
    
    Mold=M;
    
    Mvisu=M;
    Mvisu(find(d<-1))=nan;
    Mvisu(find(d>0))=nan;
    figure(1);
    colormap(map);
    imagesc(X(1:Nx,1),Y(1,1:Ny),Mvisu); hold on;
    caxis([min(Mvisu(:))-1e-6 max(Mvisu(:))]);
        
    contour(X,Y,d','LevelList',linspace(min(d(:)),-1,30),'LineColor','black','LineWidth',1);
    contour(X,Y,dt0','LevelList',[0],'LineColor','blue','LineWidth',2);
    contour(X,Y,d','LevelList',[-1],'LineColor','green','LineWidth',2);
    contour(X,Y,dt1','LevelList',[0],'LineColor','magenta','LineWidth',2);
    %contour(X,Y,phi3','LevelList',[0],'LineColor','red','LineWidth',2);
    set(axes1,'position',[0.05 0.05 0.9 0.9],'units','normalized','FontSize',15,'FontName','Arial');
    h = colorbar;
    ylabel(h, 'Mobility','FontSize',30,'FontName','Arial','Color',[0.5 0 0])
    axis image;
    hold off;
    drawnow();
    
    
    C=-dt1;
    
    [~,d2,C]=initBand(M,d+1,C,dx,dy);
    
    [C,dd]=FastMarchingExtend(d2,M,C,dx,dy);
        
    Cvisu=C;
    Cvisu(find(phi1>2*dx))=0;
    Cvisu(find(phi1<-2*dx))=0;
    
    
    %crit=0.1*max(abs(Cvisu(:)));%/max(abs(M));
    
    figure(2);
    imagesc(Cvisu); hold on;
    contour(dd,'LevelList',linspace(min(dd(:)),max(dd(:)),30),'Color','black','LineWidth',1.5);
    contour(dd,'LevelList',[0],'Color','red','LineWidth',1.5);
    contour(dd,'LevelList',[1],'Color','blue','LineWidth',1.5);
    caxis([min(C(:))-1e-10 max(C(:))]);
    axis equal;


    M=max(M-0.1*C,1e-6);
    assert(min(M(:))>0);
    
    disp(crit);
    toc();
end

save([nick '.mat']);
