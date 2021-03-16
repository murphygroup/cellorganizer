function feats=tz_calcspfeats(dists,axlns)

%function feats=tz_calcspfeats(dists,axlns)

nobj=length(dists);

dsparas=[];
lnparas=[];
lens=[];

for i=1:nobj
    sel=i;
    y=dists{sel};
    lens(i)=length(y);
    x=(0:(lens(i)-1))/(lens(i)-1); 
    
    %figure
    sp=spap2(2,4,x,y);
    fnplt(sp);
    dsparas(i,:)=[sp.knots(5),sp.coefs];
    
    hold on
    plot(x,y,'rx');
    hold off
    drawnow
    %figure
%     pause(1)
    y=axlns{sel}(:,2)';
    y1=y(x<0.33 | x>=0.67);
    y2=y(x>=0.33 & x<0.67);
    
    if mean(y1)>mean(y2)
        y=-y;
    end
    y=y-min(y);
    
    sp=spap2(2,4,x,y)
    fnplt(sp);
    lnparas(i,:)=[sp.knots(5),sp.coefs];
   
    hold on
    plot(x,y,'rx');
    hold off
    drawnow
%     pause(1)
end

feats=[dsparas,lnparas,lens'];